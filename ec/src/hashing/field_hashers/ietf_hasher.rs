use crate::hashing::map_to_curve_hasher::*;
use crate::hashing::*;
use ark_ff::{Field, PrimeField};
use ark_std::vec::Vec;
use digest::{Digest, FixedOutput};

use super::get_len_per_elem;

/// This field hasher constructs a Hash-To-Field based on a fixed-output hash function,
/// like SHA2, SHA3 or Blake2.
/// The implementation aims to follow the specification in [Hashing to Elliptic Curves (draft)](https://tools.ietf.org/pdf/draft-irtf-cfrg-hash-to-curve-13.pdf).
///
/// # Examples
///
/// ```
/// use ark_test_curves::bls12_381::Fq;
/// use ark_ec::hashing::field_hashers::IETFHasher;
/// use sha2::Sha256;
/// use crate::ark_ec::hashing::map_to_curve_hasher::HashToField;
///
/// let hasher = <IETFHasher<Sha256> as HashToField<Fq>>::new_hash_to_field(&[1, 2, 3], 2).unwrap();
/// let field_elements: Vec<Fq> = hasher.hash_to_field(b"Hello, World!").unwrap();
///
/// assert_eq!(field_elements.len(), 2);
/// ```
pub struct IETFHasher<H: FixedOutput + Digest + Sized + Clone> {
    hasher: H,
    count: usize,
    domain: Vec<u8>,
    len_per_base_elem: usize,
}

impl<F: Field, H: FixedOutput + Digest + Sized + Clone> HashToField<F> for IETFHasher<H> {
    fn new_hash_to_field(domain: &[u8], count: usize) -> Result<Self, HashToCurveError> {
        // TODO check whether this is the same as
        // hasher.update(self.domain).update(self.domain.len() as u8);
        let mut dst_prime: Vec<u8> = domain.to_vec();
        dst_prime.push(domain.len() as u8);

        // Create hasher and map the error type
        let hasher = H::new();

        // Hardcoded security parameter, defined as `k` in:
        // <https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-10.html#name-security-considerations-3>.
        // TODO this could be parametrised
        const SECURITY_PARAMETER: usize = 128;

        // The final output of `hash_to_field` will be an array of field
        // elements from F::BaseField, each of size `len_per_elem`.
        let len_per_base_elem = get_len_per_elem::<F, SECURITY_PARAMETER>();

        Ok(IETFHasher {
            hasher,
            count,
            domain: dst_prime,
            len_per_base_elem,
        })
    }

    fn hash_to_field(&self, message: &[u8]) -> Result<Vec<F>, HashToCurveError> {
        // output size of the hash function, e.g. 32 bytes = 256 bits for sha2::Sha256
        let b_in_bytes: usize = H::output_size();

        let m = F::extension_degree() as usize;

        // The user imposes a `count` of elements of F_p^m to output per input msg,
        // each field element comprising `m` BasePrimeField elements.
        let len_in_bytes = self.count * m * self.len_per_base_elem;

        // total length of output divided by output size of the hash function
        let ell: usize = (len_in_bytes + b_in_bytes - 1) / b_in_bytes;

        // Input block size of sha256
        const S_IN_BYTES: usize = 64;
        // Represent `len_in_bytes` as a 2-byte array.
        // As per I2OSP method outlined in https://tools.ietf.org/pdf/rfc8017.pdf,
        // The program should abort if integer that we're trying to convert is too large.
        assert!(len_in_bytes < 1<<16);
        let l_i_b_str: [u8; 2] = (len_in_bytes as u16).to_be_bytes();

        let mut uniform_bytes: Vec<u8> = Vec::with_capacity(len_in_bytes);

        let mut b_0_hasher = self.hasher.clone();
        b_0_hasher.update(&[0; S_IN_BYTES]);
        b_0_hasher.update(message);
        b_0_hasher.update(&l_i_b_str);
        b_0_hasher.update(&[0]);
        b_0_hasher.update(&self.domain);
        let b_0 = b_0_hasher.finalize().to_vec();

        let mut b_1_hasher = self.hasher.clone();
        b_1_hasher.update(&b_0);
        b_1_hasher.update(&[1]);
        b_1_hasher.update(&self.domain);
        let b_1 = b_1_hasher.finalize().to_vec();

        uniform_bytes.extend(b_1.clone());

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
            uniform_bytes.extend(b_i.clone());
        }

        let mut output: Vec<F> = Vec::with_capacity(self.count);
        for i in 0..self.count {
            let mut base_prime_field_elems: Vec<F::BasePrimeField> = Vec::new();
            for j in 0..m {
                let elm_offset = self.len_per_base_elem * (j + i * m);
                let val: F::BasePrimeField = F::BasePrimeField::from_be_bytes_mod_order(
                    &uniform_bytes[elm_offset..elm_offset + self.len_per_base_elem],
                );
                base_prime_field_elems.push(val);
            }
            let f: F = F::from_base_prime_field_elems(&base_prime_field_elems).unwrap();
            output.push(f);
        }

        Ok(output)
    }
}
