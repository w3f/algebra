use crate::hashing::map_to_curve_hasher::*;
use crate::hashing::*;
use ark_ff::{Field, PrimeField};
use ark_std::vec::Vec;
use digest::{Digest, FixedOutput};

// This function computes the length in bytes that a hash function should output
// for hashing a field element.
// See section 5.1 and 5.3 of the
// [IETF hash standardization draft](https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-10)
fn get_len_per_elem<F: Field>(security_parameter: usize) -> usize {
    // ceil(log(p))
    let base_field_size_in_bits = F::BasePrimeField::size_in_bits();
    // ceil(log(p)) + security_parameter
    let base_field_size_with_security_padding_in_bits =
        base_field_size_in_bits + security_parameter;
    // ceil( (ceil(log(p)) + security_parameter) / 8)
    let bytes_per_base_field_elem =
        ((base_field_size_with_security_padding_in_bits + 7) / 8) as u64;
    bytes_per_base_field_elem as usize
}

//@Skalman: not sure why do we need to that for general bytes that only need to be done by the digest.
fn map_bytes_to_field_elem<F: Field>(bz: &[u8]) -> Option<F> {
    let d: usize = F::extension_degree() as usize;
    let len_per_elem = bz.len() / d;
    let mut base_field_elems = Vec::new();
    for i in 0..d {
        let val = F::BasePrimeField::from_be_bytes_mod_order(
            &bz[i * len_per_elem..(i + 1) * len_per_elem],
        );
        base_field_elems.push(val);
    }
    F::from_base_prime_field_elems(&base_field_elems)
}

/// This field hasher constructs a Hash to Field from a variable output hash function.
/// It handles domains by hashing the input domain into 256 bits, and prefixes
/// these bits to every message it computes the hash of.
/// The state after prefixing the domain is cached.
///
/// # Examples
///
/// ```
/// use ark_test_curves::bls12_381::Fq;
/// use ark_ec::hashing::field_hashers::DefaultFieldHasher;
/// use sha2::Sha256;
/// use crate::ark_ec::hashing::map_to_curve_hasher::HashToField;
///
/// let hasher = <DefaultFieldHasher<Sha256> as HashToField<Fq>>::new_hash_to_field(&[1, 2, 3], 2).unwrap();
/// let field_elements: Vec<Fq> = hasher.hash_to_field(b"Hello, World!").unwrap();
///
/// assert_eq!(field_elements.len(), 2);
/// ```
pub struct DefaultFieldHasher<H: FixedOutput + Digest + Sized + Clone> {
    // This hasher should already have the domain applied to it.
    hasher: H,
    count: usize,
    domain: Vec<u8>,
}

// Implement HashToField from F and a variable output hash
impl<F: PrimeField, H: FixedOutput + Digest + Sized + Clone> HashToField<F>
    for DefaultFieldHasher<H>
{
    fn new_hash_to_field(domain: &[u8], count: usize) -> Result<Self, HashToCurveError> {
        // TODO check whether this is the same as
        // hasher.update(self.domain).update(self.domain.len() as u8);
        let mut dst_prime: Vec<u8> = domain.to_vec();
        dst_prime.push(domain.len() as u8);

        // Create hasher and map the error type
        let hasher = H::new();

        Ok(DefaultFieldHasher {
            hasher,
            count,
            domain: dst_prime,
        })
    }

    fn hash_to_field(&self, message: &[u8]) -> Result<Vec<F>, HashToCurveError> {
        // output size of the hash function, e.g. 32 bytes = 256 bits for sha2::Sha256
        let b_in_bytes: usize = H::output_size();
        // Hardcoded security parameter, defined as `k` in:
        // <https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-10.html#name-security-considerations-3>.
        // TODO this could be parametrised
        let security_parameter = 128;

        // The final output of `hash_to_field` will be an array of field
        // elements from F::BaseField, each of size `len_per_elem`.
        let len_per_elem = get_len_per_elem::<F>(security_parameter);
        // There are a total of `count` elements of F_p^m,
        // each comprising `m` base field elements.
        let len_in_bytes = self.count * F::extension_degree() as usize * len_per_elem;

        // TODO ensure ell is ceil(expression)
        let ell: usize = len_in_bytes / b_in_bytes;

        // Input block size of sha256
        const S_IN_BYTES: usize = 64;
        // TODO figure this out
        const l_i_b_str: usize = 128;

        let mut uniform_bytes: Vec<Vec<u8>> = Vec::with_capacity(ell);
        let mut output: Vec<F> = Vec::with_capacity(ell);

        let mut b_0_hasher = self.hasher.clone();
        b_0_hasher.update(&[0; S_IN_BYTES]);
        b_0_hasher.update(message);
        b_0_hasher.update(&[0, (l_i_b_str) as u8]);
        b_0_hasher.update(&[0]);
        b_0_hasher.update(&self.domain);
        let b_0 = b_0_hasher.finalize().to_vec();

        let mut b_1_hasher = self.hasher.clone();
        b_1_hasher.update(&b_0);
        b_1_hasher.update(&[1]);
        b_1_hasher.update(&self.domain);
        let b_1 = b_1_hasher.finalize().to_vec();

        uniform_bytes.push(b_1.clone());

        let mut b_i = b_1.clone();
        for i in 2..=ell {
            let mut b_i_hasher = self.hasher.clone();
            // update the hasher with xor of b_0 and b_i elements
            for (l, r) in b_0.iter().zip(b_i.iter()) {
                b_i_hasher.update(&[*l ^ *r]);
            }
            b_i_hasher.update(&[i as u8]);
            b_i_hasher.update(&self.domain);
            // redefine latest b_i
            b_i = b_i_hasher.finalize().to_vec();
            uniform_bytes.push(b_i.clone());
        }

        for bytes in uniform_bytes.iter() {
            let f = <F as PrimeField>::from_be_bytes_mod_order(bytes);
            output.push(f);
        }

        Ok(output)
    }
}

#[cfg(test)]
mod default_hasher_tests {
    use super::*;
    use ark_test_curves::bls12_381::{Fq, Fq2};

    #[test]
    fn test_get_len_per_elem() {
        let fq_len = get_len_per_elem::<Fq>(128);
        let fq2_len = get_len_per_elem::<Fq2>(128);
        assert_eq!(fq_len, fq2_len);
        assert_eq!(fq_len, 64);
    }
}
