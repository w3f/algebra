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
/// # use ark_test_curves::bls12_381::Fq;
/// # use ark_ec::hashing::field_hashers::DefaultFieldHasher;
/// # use blake2::VarBlake2b;
/// # use crate::ark_ec::hashing::map_to_curve_hasher::HashToField;
///
/// let hasher = <DefaultFieldHasher<VarBlake2b> as HashToField<Fq>>::new_hash_to_field(&[1, 2, 3], 2).unwrap();
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
        // Assume that the field size is 32 bytes and k is 256, where k is defined in
        // <https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-10.html#name-security-considerations-3>.
        const CHUNKLEN: usize = 64;

        let b_in_bytes: usize = H::output_size();
        // Hardcode security parameter, k
        let security_parameter = 128;
        let len_per_elem = get_len_per_elem::<F>(security_parameter);
        let len_in_bytes = self.count * F::extension_degree() as usize * len_per_elem;

        // TODO ensure ell is ceil(...)
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
        let mut b_0: Vec<u8> = Vec::with_capacity(b_in_bytes);
        b_0.copy_from_slice(&b_0_hasher.finalize());
        b_0_hasher.update(&self.domain);

        let mut b_1_hasher = self.hasher.clone();
        b_1_hasher.update(&b_0);
        b_1_hasher.update(&[1]);
        let mut b_1: Vec<u8> = Vec::with_capacity(b_in_bytes);
        b_1.copy_from_slice(&b_1_hasher.finalize());
        b_1_hasher.update(&self.domain);

        uniform_bytes.push(b_0);
        uniform_bytes.push(b_1);

        // let mut b_i = b_1.clone();
        for i in 2..=ell {
            let mut b_i_hasher = self.hasher.clone();
            // zip b_0 and b_i
            for (l, r) in uniform_bytes[0].iter().zip(uniform_bytes[i - 1].iter()) {
                b_i_hasher.update(&[*l ^ *r]);
            }
            b_i_hasher.update(&[i as u8]);
            let mut b_i = Vec::with_capacity(b_in_bytes);
            b_i.copy_from_slice(&b_i_hasher.finalize());
            // let b_i = b_i_hasher.finalize();
            uniform_bytes.push(b_i);
            b_i_hasher.update(&self.domain);
        }

        for (big, buf) in uniform_bytes.iter().zip(output.iter_mut()) {
            // let mut little = [0u8; CHUNKLEN];
            // little.copy_from_slice(big);
            // little.reverse();
            // *buf = F::from_be_bytes_mod_order(&little);
            *buf = <F as PrimeField>::from_be_bytes_mod_order(big);
        }
        Ok(output)
    }
}
