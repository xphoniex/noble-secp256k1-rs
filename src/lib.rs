#![no_std]
#![allow(non_snake_case)]
#![allow(clippy::upper_case_acronyms)]
//#![feature(const_trait_impl)]
//#![feature(const_mut_refs)]
//#![feature(const_option)]

mod point;
pub use point::Point;
//
//mod architecture;
//use architecture::{WORD, WORD_COUNT};
//
//pub type BigNum = Bn<WORD, WORD_COUNT>;

use awint::{inlawi_ty, Bits, InlAwi};

pub use awint;

pub type BigNum = inlawi_ty!(512);

mod curve;
//pub use curve::SECP256K1;
pub use curve::Curve;

//mod jacobian;
//
//#[cfg(feature = "std")]
//pub fn as_hex(num: Bn<WORD, WORD_COUNT>) -> [u8; 128] {
//    let mut buf = [0u8; 128];
//
//    num.to_hex_str(&mut buf).unwrap();
//
//    buf
//}
//
//#[cfg(not(feature = "std"))]
//pub fn as_hex(num: Bn<WORD, WORD_COUNT>, buf: &mut [u8; 128]) {
//    num.to_hex_str(buf).unwrap();
//}
