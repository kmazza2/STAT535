fn main() {
    println!("Hello, world!");
}

fn trapezoid(f: fn(f64) -> f64, a: f64, b: f64, eps: f64, max_steps: u64) -> Result<f64, String> {
    assert!(b > a);
    let mut old_val: f64 = 0.0;
    let mut val: f64 = 0.0;
    let mut converged: bool = false;
    for i in 1..=max_steps {
        if i == 1 {
            val = 0.5 * (b - a) * (f(a) + f(b));
        } else {
            let num_points = 2_u64.pow((i - 2) as u32);
            let h = (b - a) / (num_points as f64);
            val = 0.5 * (val + (b - a) * (0..num_points).map(|j| f(a + ((j as f64) + 0.5) * h)).sum::<f64>() / (num_points as f64))
        }
        if i > 5 &&
                (
                    (val - old_val).abs() < eps * old_val.abs() ||
                    (old_val == 0.0 && val == 0.0)
                ) {
            converged = true;
            break;
        }
        old_val = val;
    }
    match converged {
        true => Ok(val),
        false => Err(String::from("failed to converge"))
    }
}

#[test]
fn test_trapezoid() {
    match trapezoid(|_| 1.0, 2.0, 21.0, 1.0E-13_f64, 10000) {
        Ok(val) => assert!((19.0 - val).abs() < 0.00001),
        Err(s) => panic!("{}", s)
    };
    match trapezoid(|x| x.cos(), 0.35, 11.1, 1.0E-13_f64, 10000) {
        Ok(val) => assert!((-1.337450395659441 - val).abs() < 0.00001),
        Err(s) => panic!("{}", s)
    };
    match trapezoid(|x| x.exp(), 2.39, 5.28, 1.0E-13_f64, 10000) {
        Ok(val) => assert!((185.4563803458403 - val).abs() < 0.00001),
        Err(s) => panic!("{}", s)
    };
}

fn log_gamma(x: f64) -> f64 {
    assert!(x > 0.0);
    const COEF: [f64; 14] = [
        57.1562356658629235_f64, -59.5979603554754912_f64, 14.1360979747417471_f64,
        -0.491913816097620199_f64, 0.339946499848118887E-4_f64, 0.465236289270485756E-4_f64,
        -0.983744753048795646E-4_f64, 0.158088703224912494E-3_f64, -0.210264441724104883E-3_f64,
        0.217439618115212643E-3_f64, -0.164318106536763890E-3_f64, 0.844182239838527433E-4_f64,
        -0.261908384015814087E-4_f64, 0.368991826595316234E-5_f64];
    const STP: f64 = 2.5066282746310005_f64;
    let mut y = x;
    let tmp = x + 5.24218750000000000_f64;
    let tmp = (x + 0.5) * tmp.ln() - tmp;
    let mut ser = 0.999999999999997092_f64;
    for j in 0..14 {
        y += 1.0;       
        ser += COEF[j] / y; 
    }                       
    tmp + (STP * ser / x).ln()
}       

fn nchoosek(n: usize, k: usize) -> usize {
    (log_gamma((n + 1) as f64) - log_gamma((k + 1) as f64) - log_gamma((n - k + 1) as f64))
        .exp()              
        .round() as usize       
}         

fn beta(x1: f64, x2: f64) -> f64 {
    (log_gamma(x1) + log_gamma(x2) - log_gamma(x1 + x2)).exp()
}

#[test]
fn test_beta() {
    assert!((0.1158881575662717E-3_f64 - beta(5.1_f64, 9.4_f64)).abs() < 1.0E-9_f64);
}
