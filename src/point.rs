use awint::{cc, inlawi, Bits, InlAwi};

use crate::BigNum;

#[derive(Clone)]
pub struct Point {
    pub x: BigNum,
    pub y: BigNum,
}

#[allow(dead_code)]
impl Point {
    #[inline]
    fn zero() -> Point {
        Point {
            x: BigNum::zero(),
            y: BigNum::zero(),
        }
    }

    /*

    fn double(&self, P: &BigNum) -> Point {
        let X1 = self.x;
        let Y1 = self.y;

        let two = BigNum::from_u8(2).unwrap();
        let three = BigNum::from_u8(3).unwrap();

        //println!("> lam");

        //println!("3 x X1 x X1 = {}", as_hex(three * X1 * X1));
        //print!("mod_inverse = ");
        //print!("{}\n", as_hex(Self::invert(two * Y1, P)));
        //let tmp_inv = (two * Y1) ^ Bn::from_i8(-1).unwrap();
        //print!("{}\n", as_hex(((two * Y1) ^ Bn::from_i8(-1).unwrap()) % P));
        //print!("{}\n", as_hex(two * Y1));

        //println!("lam[..] = {}", as_hex(three * X1 % P * X1 % P));
        let lam = ((three * ((X1 * X1) % P)) % P) * Self::invert(&(two * Y1), P) % P;
        //println!("lam = {}", as_hex(lam));
        //println!("< lam");

        //println!("> x3");
        let X3 = (lam * lam - two * X1) % P;
        //println!("< x3");
        //println!("> y3");
        //print!("X1 {} - X3 {} = ", as_hex(X1), as_hex(X3));
        //println!("{}", as_hex(X1 + P - X3));
        let X1 = if X1 < X3 { X1 + P } else { X1 };
        let Y3 = (lam * (X1 - X3) - Y1) % P;
        //println!("< y3");

        Point { x: X3, y: Y3 }
    }
    */

    // assumes `x1`, `y1`, and `p` are all 512 bits
    fn double_assign_mod(&mut self, p: &BigNum) {
        let x1 = self.x;
        let y1 = self.y;
        //let p = &SECP256K1.P;

        // in this first part, we are calculating
        // let lam = ((3 * ((x1 * x1) % p)) % p) * Self::invert(&(two * Y1), P) % P;

        // int0 = (x1 * x1) % p
        let mut pad = inlawi!(0u512);
        let mut x1_copy = inlawi!(x1; ..512).unwrap();
        x1_copy.mul_assign(&x1, &mut pad).unwrap();
        // just renaming, this is zero cost if we don't use x1_copy later
        let x1_squared = x1_copy;
        let mut int0 = inlawi!(0u512);
        Bits::udivide(&mut pad, &mut int0, &x1_squared, &p).unwrap();
        // x1_squared, pad, int0 are live (192 bytes)

        // int1 = (3 * into) % p
        int0.short_cin_mul(0, 3);
        let mut int1 = x1_squared;
        Bits::udivide(&mut pad, &mut int1, &int0, &p).unwrap();
        // int1, pad, int0 are live (192 bytes)

        let mut y1_copy = inlawi!(y1; ..512).unwrap();
        // equivalent to multiplying by 2
        y1_copy.shl_assign(1).unwrap();
        let mut y1_mul2 = y1_copy;
        // I don't know what the `Self::invert` function is, suppose we calculated it from `y1_mul2`
        // and `p` and stored it in `invert_result`, here let me just use `y1_mul2`.
        let mut invert_result = Self::invert(&mut y1_mul2, p);
        // int1, pad, int0, invert_result are live (192 bytes)

        // int2 = (int1 * invert_result) % p;
        invert_result.mul_assign(&int1, &mut pad).unwrap();
        let mut int2 = int0;
        Bits::udivide(&mut pad, &mut int2, &invert_result, &p).unwrap();
        let mut lam = int2;
        // int1, pad, lam, invert_result are live (192 bytes)

        // in this next part we are calculating
        /*let X3 = (lam * lam - two * X1) % P;
        let Y3 = if *X1 < X3 {
            (lam * (*P + X1 - X3) - Y1) % P
        } else {
            (lam * (*X1 - X3) - Y1) % P
        };

        self.x = X3;
        self.y = Y3;*/

        // (lam * lam - two * x1) % p
        let mut lam_copy = invert_result;
        cc!(lam; lam_copy).unwrap();
        lam_copy.mul_assign(&lam, &mut pad).unwrap();
        let mut lam_squared = lam_copy;
        lam_squared.sub_assign(&x1).unwrap();
        // we would have to copy anyway, just subtract twice
        lam_squared.sub_assign(&x1).unwrap();
        let mut x3 = int1;
        Bits::udivide(&mut pad, &mut self.x, &lam_squared, &p).unwrap();
        // x3, pad, lam, lam_squared are live (192 bytes)

        // (
        //   ((x1 - x3 + (`p` if x1 < x3)) * lam)
        //   - y1
        // ) % p
        let mut tmp = lam_squared;
        cc!(x1; tmp).unwrap();
        tmp.sub_assign(&self.x).unwrap();
        if x1.ult(&self.x).unwrap() {
            tmp.add_assign(&p).unwrap();
        }
        lam.mul_assign(&tmp, &mut pad).unwrap();
        lam.sub_assign(&y1).unwrap();
        let mut y3 = tmp;
        Bits::udivide(&mut pad, &mut self.y, &lam, &p).unwrap();
        //Bits::udivide(&mut pad, &mut y3, &lam, &p).unwrap();

        //self.x = x3;
        //self.y = y3;

        // x3 and y3 are our results
        //dbg!(x3, y3);
    }

    /*
        #[inline]
        fn double_assign_mod(&mut self, P: &BigNum, one: &BigNum) {
            //let X1 = self.x.clone();
            let X1 = &self.x;
            //let Y1 = self.y.clone();
            let Y1 = &self.y;

            let two = BigNum::from_u8(2).unwrap();
            let three = BigNum::from_u8(3).unwrap();
            //let three = BigNum::from_le_bytes(&[
            //    3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            //    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            //    0, 0, 0, 0, 0, 0,
            //]);

            //let mut three = BigNum::from_le_bytes(&[
            //    3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            //    0, 0, 0,
            //]);

            let lam = ((three * ((*X1 * *X1) % P)) % P) * Self::invert(&(two * Y1), P) % P;

            let X3 = (lam * lam - two * X1) % P;
            let Y3 = if *X1 < X3 {
                (lam * (*P + X1 - X3) - Y1) % P
            } else {
                (lam * (*X1 - X3) - Y1) % P
            };

            self.x = X3;
            self.y = Y3;
        }

        fn add(&self, other: Point, P: &BigNum) -> Point {
            let X1 = self.x;
            let Y1 = self.y;
            let X2 = other.x;
            let Y2 = other.y;

            if X1.is_zero() || Y1.is_zero() {
                return other;
            }

            if X2.is_zero() || Y2.is_zero() {
                return self.clone();
            }

            if X1 == X2 && Y1 == Y2 {
                return self.double(P);
            }

            if X1 == X2 && Y1 == BigNum::from_i8(-1).unwrap() * Y2 {
                return Self::zero();
            }

            let lam = ((Y2 - Y1) - Self::invert(&(X2 - X1), P)) % P;

            let X3 = (lam * lam - X1 - X2) % P;
            let Y3 = (lam * (X1 - X3) - Y1) % P;

            Point { x: X3, y: Y3 }
        }

        #[inline]
        fn add_assign_mod(&mut self, other: &Point, P: &BigNum, one: &BigNum) {
            let X1 = &self.x;
            let Y1 = &self.y;
            let X2 = &other.x;
            let Y2 = &other.y;

            if X1.is_zero() || Y1.is_zero() {
                self.x = *X2;
                self.y = *Y2;
                return;
            }

            if X2.is_zero() || Y2.is_zero() {
                return;
            }

            if X1 == X2 && Y1 == Y2 {
                //let double = self.double(*P);
                //self.x = double.x;
                //self.y = double.y;
                self.double_assign_mod(P, one);
                return;
            }

            if X1 == X2 && *Y1 == BigNum::from_i8(-1).unwrap() * Y2 {
                self.x = BigNum::zero();
                self.y = BigNum::zero();
                return;
            }

            //let lam = ((*Y2 - Y1) - Self::invert_one(&mut (*X2 - X1), P, one)) % P;
            let lam = ((*Y2 - Y1) - Self::invert(&mut (*X2 - X1), P)) % P;

            let X3 = (lam * lam - X1 - X2) % P;
            let Y3 = (lam * (*X1 - X3) - Y1) % P;

            self.x = X3;
            self.y = Y3;
        }
    */

    #[inline]
    fn add_assign_mod(&mut self, other: &mut Point, p: &BigNum) {
        let X1 = &self.x;
        let Y1 = &self.y;
        let X2 = &other.x;
        let Y2 = &mut other.y;
        //let p = &SECP256K1.P;

        if X1.is_zero() || Y1.is_zero() {
            self.x = *X2;
            self.y = *Y2;
            return;
        }

        if X2.is_zero() || Y2.is_zero() {
            return;
        }

        if X1 == X2 && Y1 == Y2 {
            self.double_assign_mod(p);
            return;
        }

        Y2.neg_assign(true);
        if X1 == X2 && Y1 == Y2 {
            self.x = BigNum::zero();
            self.y = BigNum::zero();
            return;
        }
        Y2.neg_assign(false);

        // lam = ((Y2 - Y1) - Self::invert((X2 - X1), P)) % P;
        let mut pad = inlawi!(0u512); // alloc
        cc!(Y2; pad).unwrap();
        pad.sub_assign(Y1).unwrap();
        let mut y2_y1 = pad.clone(); // alloc

        cc!(X2; pad).unwrap();
        pad.sub_assign(X1).unwrap(); // pad = X2 - X1

        let mut inverted_x2_x1 = Self::invert(&mut pad, p);

        y2_y1.sub_assign(&inverted_x2_x1).unwrap();

        Bits::udivide(&mut pad, &mut inverted_x2_x1, &y2_y1, &p).unwrap();
        let lam = inverted_x2_x1;
        // pad, y2_y1, lam

        //
        // X3 = (lam * lam - X1 - X2) % P;
        let mut lam_copy = y2_y1; // aliasing
        cc!(lam; lam_copy).unwrap();

        let mut lam_squared = lam;
        lam_squared.mul_assign(&lam_copy, &mut pad).unwrap();
        lam_squared.sub_assign(&X1).unwrap();
        lam_squared.sub_assign(&X2).unwrap();
        let lam_squared_x1_x2 = lam_squared;
        let mut rem = lam_copy; // alias

        Bits::udivide(&mut pad, &mut rem, &lam_squared_x1_x2, &p).unwrap();
        // pad, rem, lam_squared_x1_x2

        //
        // let Y3 = (lam * (*X1 - X3) - Y1) % P;
        let mut y3 = lam_squared_x1_x2;
        cc!(X1; y3).unwrap();
        y3.sub_assign(&self.x).unwrap(); // y3 = X1 - X3
        y3.mul_assign(&lam_copy, &mut pad).unwrap(); // y3 = lam * (X1 - X3)
        y3.sub_assign(Y1).unwrap(); // y3 = (lam * (X1 - X3) - Y1)
        Bits::udivide(&mut pad, &mut self.y, &y3, &p).unwrap();
    }

    #[inline]
    pub fn invert(num: &mut BigNum, P: &BigNum) -> BigNum {
        /*
        let two = BigNum::from_u8(2).unwrap();
        power_modulo(num, &mut (*modulo - two), modulo) % modulo
        */
        power_modulo(num, P)
    }

    #[inline]
    pub fn multiply_DA(&self, n: &mut BigNum, P: &BigNum) -> Point {
        let mut p = Self::zero();
        let mut d = self.clone();

        while !n.is_zero() {
            if n.lsb() {
                p.add_assign_mod(&mut d, P);
            }

            d.double_assign_mod(P);
            n.ashr_assign(1).unwrap();
        }

        p
    }

    /*
    pub fn multiply_JDA(&self, mut n: BigNum, P: &BigNum) -> Point {
        use crate::jacobian::JacobianPoint;

        let mut p = JacobianPoint::zero();
        let mut d = JacobianPoint::from_affine(self.clone());

        while !n.is_zero() {
            if n & BigNum::one() > BigNum::zero() {
                p = p.add(d.clone(), P);
            }
            let d_d = d.double(P);
            d = d_d;
            n >>= 1_usize;
        }

        //println!("to affine");

        p.to_affine(P)
    }

    pub fn as_hex(&self, buf: &mut [u8; 128]) {
        self.x.to_hex_str(&mut buf[0..64]).unwrap();
        self.y.to_hex_str(&mut buf[64..128]).unwrap();
    }
    */
}

#[cfg(feature = "std")]
impl core::fmt::Display for Point {
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        let _ = write!(f, "{{\n");
        let _ = write!(f, "\tx: {},\n", self.x);
        let _ = write!(f, "\ty: {}\n ", self.y);
        let _ = write!(f, "}}");

        Ok(())
    }
}
/*

// calculates num ^ power % modulo, faster than naive way
#[inline]
fn power_modulo(num_: &BigNum, power: &mut BigNum, modulo: &BigNum) -> BigNum {
    let one = BigNum::one();

    let mut result = one;
    let mut num = *num_;

    loop {
        // NOTE: is_even() shaves 0.3sec from 1.75s vs. is_multiple_of(&two)
        while power.is_even() {
            *power >>= 1_usize;
            num = (num % modulo) * (num % modulo);
        }

        if !power.is_zero() {
            result = num % modulo * result % modulo;
            *power -= one;
        }

        if power.is_zero() {
            break;
        }
    }

    result
}

*/

// TODO(check if you can only use num instad of `num_` AND `num`
// calculates num ^ power % modulo, faster than naive way
// num = num
// power = P - 2
// modulo = P
//
// return should be %= P
#[inline]
fn power_modulo(num: &mut BigNum, P: &BigNum) -> BigNum {
    //dbg!(&num_);

    let mut result = inlawi!(1_u512);

    let mut power = inlawi!(2_u512);
    power.rsb_assign(P).unwrap();
    /*
     */

    /*
    let mut power = inlawi!(0_u512);
    cc!(SECP256K1.P; power).unwrap();
    power.sub_assign(&TWO).unwrap();
     */

    //dbg!(&power);

    let mut pad = inlawi!(0_u512);
    let mut clone = inlawi!(0_u512);

    //let mut num = inlawi!(0_u512);
    //cc!(num_; num);

    loop {
        // power is even
        while !power.lsb() {
            power.ashr_assign(1).unwrap();
            //println!(">> pwr = {}", power);

            // (num % modulo) * (num % modulo)
            //println!(">> num = {}", num);
            cc!(num; clone).unwrap();
            //println!(">> cln = {}", clone);
            Bits::udivide(&mut pad, num, &clone, P).unwrap();
            cc!(num; clone).unwrap();
            num.mul_assign(&mut clone, &mut pad).unwrap();

            //println!(">> num = {}", num);
        }

        if !power.is_zero() {
            // result = num % modulo * result % modulo
            //
            // OR
            //
            // result %= modulo
            // result *= (num % modulo)
            //
            // 1. result = result % modulo

            cc!(result; clone).unwrap();
            Bits::udivide(&mut pad, &mut result, &clone, P).unwrap();

            // 2. clone =  num % modulo

            cc!(num; clone).unwrap();
            Bits::udivide(&mut pad, &mut clone, &num, P).unwrap();

            // 3. result *= clone
            result.mul_assign(&mut clone, &mut pad).unwrap();

            // power -= 1
            //power.sub_assign(&ONE).unwrap();
            power.dec_assign(false);
        }

        if power.is_zero() {
            break;
        }
    }

    cc!(result; clone).unwrap();
    Bits::udivide(&mut pad, &mut result, &clone, P).unwrap();
    //dbg!(&result);

    result
}
