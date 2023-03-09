use awint::{inlawi, Bits, InlAwi};

use crate::{BigNum, Point};

pub struct Curve {
    pub G: Point,
    P: BigNum,
    //n: BigNum,
}

impl Curve {
    pub fn secp256k1() -> Curve {
        let G = Point {
            x: inlawi!(0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798_u512),
            y: inlawi!(0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8_u512),
        };
        let P = inlawi!(0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f_u512);
        //let n = inlawi!(0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141_u512);

        //Curve { G, P, n }
        Curve { G, P }
    }

    /// Returns `G * n mod P` using simple multiplication.
    ///
    /// This method is slow, use jacobian multiplcation [`multiply`](Self::multiply).
    #[inline]
    pub fn multiply_simple(&self, n: &mut BigNum) -> Point {
        self.G.multiply_DA(n, &self.P)
    }
}

/*
pub struct Curve {
    pub G: Point,
    P: BigNum,
    n: BigNum,
}

impl Curve {
    /// Creates a secp256k1 curve with `G`, `P`, and `n` values.
    #[inline]
    pub fn secp256k1() -> Curve {
        /*
        use fixed_bigint::num_traits::Num;
        let Gx = BigNum::from_str_radix(
            "79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798",
            16,
        )
        .unwrap();
        */

        // prefix = 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,

        let Gx = BigNum::from_le_bytes(&[
            0x98, 0x17, 0xf8, 0x16, 0x5b, 0x81, 0xf2, 0x59, 0xd9, 0x28, 0xce, 0x2d, 0xdb, 0xfc,
            0x9b, 0x02, 0x07, 0x0b, 0x87, 0xce, 0x95, 0x62, 0xa0, 0x55, 0xac, 0xbb, 0xdc, 0xf9,
            0x7e, 0x66, 0xbe, 0x79,
        ]);

        /*
         */

        /*
        let Gy = BigNum::from_str_radix(
            "483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8",
            16,
        )
        .unwrap();
        */

        let Gy = BigNum::from_le_bytes(&[
            0xb8, 0xd4, 0x10, 0xfb, 0x8f, 0xd0, 0x47, 0x9c, 0x19, 0x54, 0x85, 0xa6, 0x48, 0xb4,
            0x17, 0xfd, 0xa8, 0x08, 0x11, 0x0e, 0xfc, 0xfb, 0xa4, 0x5d, 0x65, 0xc4, 0xa3, 0x26,
            0x77, 0xda, 0x3a, 0x48,
        ]);
        /*
         */

        /*
        let P = BigNum::from_str_radix(
            "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f",
            16,
        )
        .unwrap();
        */

        let P = BigNum::from_le_bytes(&[
            0x2f, 0xfc, 0xff, 0xff, 0xfe, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
            0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
            0xff, 0xff, 0xff, 0xff,
        ]);
        /*
         */

        /*
        let n = BigNum::from_str_radix(
            "fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141",
            16,
        )
        .unwrap();
        */

        let n = BigNum::from_le_bytes(&[
            0x41, 0x41, 0x36, 0xd0, 0x8c, 0x5e, 0xd2, 0xbf, 0x3b, 0xa0, 0x48, 0xaf, 0xe6, 0xdc,
            0xae, 0xba, 0xfe, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
            0xff, 0xff, 0xff, 0xff,
        ]);
        /*
         */

        let G = Point { x: Gx, y: Gy };

        Curve { G, P, n }
    }

    /// Returns `G * n mod P` using jacobian multiplicatioin.
    #[inline]
    pub fn multiply(&self, n: &mut BigNum) -> Point {
        #[cfg(feature = "jacobian")]
        {
            self.G.multiply_JDA(n, self.P)
        }
        #[cfg(not(feature = "jacobian"))]
        {
            self.G.multiply_DA(n, &self.P)
        }
    }

    /// Returns `G * n mod P` using simple multiplication.
    ///
    /// This method is slow, use jacobian multiplcation [`multiply`](Self::multiply).
    #[inline]
    pub fn multiply_simple(&mut self, n: &mut BigNum) -> &Point {
        self.G.multiply_DA(n, &self.P);

        &self.G
    }

    /// Signs a message `m` using a random nonce `k` for private key `d`.
    //#[inline]
    pub fn sign(&self, m: BigNum, d: BigNum, k: &mut BigNum) -> Point {
        let P = self.P;
        let n = &self.n;

        let mut kG = {
            #[cfg(feature = "jacobian")]
            {
                self.G.multiply_JDA(k, P)
            }
            #[cfg(not(feature = "jacobian"))]
            {
                self.G.multiply_DA(k, &P)
            }
        };

        // r
        kG.x %= n;

        // s
        kG.y = (Point::invert(k, n) % n) * ((m + d * kG.x) % n) % n;

        kG
    }

    /// Signs a message `m` using a random `k` nonce for private key `d`.
    ///
    /// This method returns the canonical view of `s` in signature, meaning
    /// if `s` is greater than half of `n`, we subtract `n` from it.
    //#[inline]
    pub fn sign_canonical(&self, m: BigNum, d: BigNum, k: &mut BigNum) -> Point {
        let P = self.P;
        let n = &self.n;

        let mut kG = {
            #[cfg(feature = "jacobian")]
            {
                self.G.multiply_JDA(k, P)
            }
            #[cfg(not(feature = "jacobian"))]
            {
                self.G.multiply_DA(k, &P)
            }
        };

        // r
        kG.x %= n;

        // s
        kG.y = (Point::invert(k, n) % n) * ((m + d * kG.x) % n) % n;

        // canonical
        if kG.y > *n / BigNum::from_u8(2).unwrap() {
            kG.y = *n - kG.y;
        }

        kG
    }

    ///
    pub fn verify(&self) -> bool {
        todo!();
    }
}
*/
