use num_integer::Integer; 
 
pub fn egcd<T: Copy + Integer>(a: T, b: T) -> (T, T, T) { 
    if a == T::zero() { 
        (b, T::zero(), T::one()) 
    } 
    else { 
        let (g, x, y) = egcd(b % a, a); 
        (g, y - (b / a) * x, x) 
    } 
} 
 
pub fn modinverse<T: Copy + Integer>(a: T, m: T) -> Option<T> { 
    let (g, x, _) = egcd(a, m); 
    if g != T::one() { 
        None 
    } 
    else { 
        Some((x % m + m) % m) 
    } 
} 
#[derive(Debug)]
struct ECpoint {
    curve: String, // por ejemplo
    x: i64,
    y: i64,
    z: i64,
}

impl ECpoint {
    fn new(curve: String, x: i64, y: i64, z: i64) -> Self {
        Self {
            curve,
            x,
            y,
            z,
        }
    }
}

#[allow(dead_code)]
struct ECcurve { 
    p: i64, 
    a: i64, 
    b: i64, 
} 
 
impl ECcurve { 
    fn field_mul(&self, a: i64, b: i64) -> i64 { 
        (a * b) % self.p 
    } 
     
    fn field_div(&self, num: i64, den: i64) -> i64 { 
        let inverse_den = modinverse(den % self.p, self.p); 
        match inverse_den { 
            Some(inv) => self.field_mul(num % self.p, inv), 
            None => -1, 
        } 
    }
    
    fn field_exp(&self, num: i64, power: i64) -> i64 {
        num.pow(power as u32) % self.p
    }

    fn identity(&self) -> ECpoint {
        ECpoint { curve: String::new(), x: self.p, y: 0, z: 0 }
    }
    fn touches(&self, q: &ECpoint) -> bool {
        let x = q.x;
        let y = q.y;
        let y2 = (y * y) % self.p;
        let x3ab = (self.field_mul((x * x) % self.p + self.a, x) + self.b) % self.p;
        y2 == (x3ab) % self.p
    }

    fn tangent(&self, q: &ECpoint) -> i64 {
        self.field_div(q.x * q.x * 3 + self.a, q.y * 2)
    }
    fn double(&self, q: ECpoint) -> ECpoint {
        if q.x == self.p {
            return q;
        } 
    
        let s = (4*q.x*q.y.pow(2))%self.p;
        let z2 = q.z.pow(2);
        let z4 = (z2.pow(2))%self.p;
        let m = 3*q.x.pow(2) + self.a*z4; 
        let x = (m.pow(2) - 2*s)%self.p;
        let y2 = q.y.pow(2);
        let y = (m*(s-x) - 8*y2.pow(2))%self.p; 
        let z = (2*q.y*q.z)%self.p;
        ECpoint{ curve: String::new(), x, y, z }
    }
    fn add(&self, q1: ECpoint, q2: ECpoint) -> ECpoint {
        // Identity special cases
        if q1.x == self.p {
            return q2;
        }
        if q2.x == self.p {
            return q1;
        }
        let q1z2 = q1.z * q1.z;
        let q2z2 = q2.z * q2.z;
        let xs1 = self.field_mul(q1.x, q2z2);
        let xs2 = self.field_mul(q2.x, q1z2);
        let ys1 = self.field_mul(self.field_mul(q1.y, q2z2), q2.z);
        let ys2 = self.field_mul(self.field_mul(q2.y, q1z2), q1.z);
    
        // Equality special cases
        if xs1 == xs2 {
            if ys1 == ys2 {
                return self.double(q1);
            } else {
                return self.identity();
            }
        }
    
        // Ordinary case
        let xd = (xs2 - xs1) % self.p;
        let yd = (ys2 - ys1) % self.p;
        let xd2 = (xd * xd) % self.p;
        let xd3 = (xd2 * xd) % self.p;
        let x = (self.field_mul(yd, yd) - xd3 - self.field_mul(2 * xs1, xd2)) % self.p;
        let y = (self.field_mul(yd, self.field_mul(xs1, xd2) - x) - self.field_mul(ys1, xd3)) % self.p;
        let z = self.field_mul(xd, self.field_mul(q1.z, q2.z)) % self.p;
        ECpoint { curve: String::new(), x, y, z }
    }
    // Multiply this elliptic curve point Q by the scalar (integer) m 
    //    Often the point Q will be the generator G 
    fn mul(&self, m: u32, mut q: (f64, f64)) -> ECpoint {
        let mut r = self.identity(); // return point
    
        let mut m = m;
    
        while m != 0 { // binary multiply loop
            if m & 1 == 1 { // bit is set
                println!("  mul: adding q to r = {:?}", r);
                r = self.add(r, ECpoint { curve: String::new(), x: q.0 as i64, y: q.1 as i64, z: 1 });
                println!("  mul: added q to r = {:?}", r);
            }
            m >>= 1;
            if m != 0 {
                println!("  mul: doubling q = {:?}", q);
                q = (self.double(ECpoint { curve: String::new(), x: q.0 as i64, y: q.1 as i64, z: 1 }).x as f64, self.double(ECpoint { curve: String::new(), x: q.0 as i64, y: q.1 as i64, z: 1 }).y as f64);
            }
        }
    
        ECpoint { curve: String::new(), x: r.x, y: r.y, z: r.z }
    }
}
fn main() {
    let curve = ECcurve { p: 17, a: 2, b: 2 };

    let a = 3;
    let b = 5;
    let c = curve.field_mul(a, b);
    println!("{} * {} (mod {}) = {}", a, b, curve.p, c);

    let d = 8;
    let e = 3;
    let f = curve.field_div(d, e);
    println!("{} / {} (mod {}) = {}", d, e, curve.p, f);

    let g = 4;
    let h = 3;
    let i = curve.field_exp(g, h);
    println!("{} ^ {} (mod {}) = {}", g, h, curve.p, i);

    let point = curve.identity();
    println!("Identity point: ({}, {})", point.x, point.y);

    let q = ECpoint { curve: String::new(), x: 1, y: 4, z: 1 };
    let touches_curve = curve.touches(&q);
    println!("Q touches curve: {}", touches_curve);

    let slope = curve.tangent(&q);
    println!("Slope of tangent at Q: {}", slope);

    let r = curve.double(q);
    println!("2Q = ({}, {})", r.x, r.y);
    
    let curve = ECcurve { p: 17, a: 2, b: 2 };
 
    // Crea dos puntos de la curva
    let q1 = ECpoint { curve: String::new(), x: 1, y: 4, z: 1 };
    let q2 = ECpoint { curve: String::new(), x: 3, y: 5, z: 1 };
 
    // Agrega los puntos
    let result = curve.add(q1, q2);
 
    // Imprime el resultado
    println!("Resultado: ({}, {})", result.x, result.y);

    let curve = ECcurve { p: 17, a: 2, b: 2 }; // define the curve to use
    let q = (5.0, 1.0); // define the point to multiply
    let m = 3; // define the scalar to multiply by
    let result = curve.mul(m, q); // call the mul function
    println!("Result: {:?}", result); // print the result
}
