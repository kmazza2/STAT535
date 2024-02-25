fn main() {
    let prior1 = PriorParams {
        alpha: 0.3,
        beta: 0.7,
    };
    let prior2 = PriorParams {
        alpha: 0.3,
        beta: 0.7,
    };
    let delta = 0.2;
    let target = 0.9;
    println!("{:<7}{:<7}{:8}{:8}", "y1", "y2", "post", "pred");
    println!("------------------------------");
    for y1 in 5..=15 {
        for y2 in 5..=y1 {
            let data1 = Data {
                y: y1,
                N: 50,
                n: 25,
            };
            let data2 = Data {
                y: y2,
                N: 50,
                n: 25,
            };
            let prob = post_prob(delta, &data1, &data2, &prior1, &prior2);
            let pred = pred_prob(target, &data1, &data2, &prior1, &prior2);
            println!("{:<7}{:<7}{:<8.4}{:<8.4}", y1, y2, prob, pred);
        }
    }
}

// This will fail if the parameters of the beta distribution
// yield singularities at the endpoints of its support; beware!
fn trapezoid(
    f: impl Fn(f64) -> f64,
    a: f64,
    b: f64,
    eps: f64,
    max_steps: u64,
) -> Result<f64, String> {
    assert!(b >= a);
    let mut old_val: f64 = 0.0;
    let mut val: f64 = 0.0;
    let mut converged: bool = false;
    for i in 1..=max_steps {
        if i == 1 {
            val = 0.5 * (b - a) * (f(a) + f(b));
        } else {
            let num_points = 2_u64.pow((i - 2) as u32);
            let h = (b - a) / (num_points as f64);
            val = 0.5
                * (val
                    + (b - a)
                        * (0..num_points)
                            .map(|j| f(a + ((j as f64) + 0.5) * h))
                            .sum::<f64>()
                        / (num_points as f64));
        }
        if i > 5 && ((val - old_val).abs() < eps * old_val.abs() || (old_val == 0.0 && val == 0.0))
        {
            converged = true;
            break;
        }
        old_val = val;
    }
    match converged {
        true => Ok(val),
        false => Err(String::from("failed to converge")),
    }
}

#[test]
fn test_trapezoid() {
    match trapezoid(|_| 1.0, 2.0, 21.0, 1.0E-13_f64, 10000) {
        Ok(val) => assert!((19.0 - val).abs() < 0.00001),
        Err(s) => panic!("{}", s),
    };
    match trapezoid(|x| x.cos(), 0.35, 11.1, 1.0E-13_f64, 10000) {
        Ok(val) => assert!((-1.337450395659441 - val).abs() < 0.00001),
        Err(s) => panic!("{}", s),
    };
    match trapezoid(|x| x.exp(), 2.39, 5.28, 1.0E-13_f64, 10000) {
        Ok(val) => assert!((185.4563803458403 - val).abs() < 0.00001),
        Err(s) => panic!("{}", s),
    };
    let inner: &dyn Fn(f64) -> f64 = &|x2| x2;
    let outer: &dyn Fn(f64) -> f64 = &|x1| {
        x1.cos()
            * trapezoid(inner, 0.0, 2.0, 1.0E-5, 100).expect("inner integral failed to converge")
    };
    match trapezoid(outer, 9.0, 13.0, 1.0E-5, 100) {
        Ok(val) => assert!((0.01609710316976865_f64 - val).abs() < 0.00001),
        Err(s) => panic!("{}", s),
    }
}

fn log_gamma(x: f64) -> f64 {
    assert!(x > 0.0);
    const COEF: [f64; 14] = [
        57.1562356658629235_f64,
        -59.5979603554754912_f64,
        14.1360979747417471_f64,
        -0.491913816097620199_f64,
        0.339946499848118887E-4_f64,
        0.465236289270485756E-4_f64,
        -0.983744753048795646E-4_f64,
        0.158088703224912494E-3_f64,
        -0.210264441724104883E-3_f64,
        0.217439618115212643E-3_f64,
        -0.164318106536763890E-3_f64,
        0.844182239838527433E-4_f64,
        -0.261908384015814087E-4_f64,
        0.368991826595316234E-5_f64,
    ];
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

fn beta(x1: f64, x2: f64) -> f64 {
    (log_gamma(x1) + log_gamma(x2) - log_gamma(x1 + x2)).exp()
}

#[test]
fn test_beta() {
    assert!((0.1158881575662717E-3_f64 - beta(5.1_f64, 9.4_f64)).abs() < 1.0E-9_f64);
}

fn nchoosek(n: u64, k: u64) -> u64 {
    (log_gamma((n + 1) as f64) - log_gamma((k + 1) as f64) - log_gamma((n - k + 1) as f64))
        .exp()
        .round() as u64
}

struct PriorParams {
    alpha: f64,
    beta: f64,
}

struct Data {
    y: u64,
    N: u64,
    n: u64,
}

fn beta_dens(x: f64, a: f64, b: f64) -> f64 {
    x.powf(a - 1.0) * (1.0 - x).powf(b - 1.0) / beta(a, b)
}

fn betabin_dens(x: u64, n: u64, a: f64, b: f64) -> f64 {
    (nchoosek(n, x) as f64) * beta((x as f64) + a, (n as f64) - (x as f64) + b) / beta(a, b)
}

fn post_prob(
    delta: f64,
    data1: &Data,
    data2: &Data,
    prior1: &PriorParams,
    prior2: &PriorParams,
) -> f64 {
    assert!(0.0_f64 < delta && delta < 1.0_f64);
    let inner: &dyn Fn(f64) -> f64 = &|p2| {
        beta_dens(
            p2,
            prior2.alpha + (data2.y as f64),
            prior2.beta + (data2.n as f64) - (data2.y as f64),
        )
    };
    let outer: &dyn Fn(f64) -> f64 = &|p1| {
        beta_dens(
            p1,
            prior1.alpha + (data1.y as f64),
            prior1.beta + (data1.n as f64) - (data1.y as f64),
        ) * trapezoid(inner, 0.0, p1 - delta, 1.0E-7, 1000)
            .expect("inner integral failed to converge")
    };
    trapezoid(outer, delta, 1.0, 1.0E-7, 1000).expect("posterior probability failed to converge")
}
fn pred_prob(
    target: f64,
    data1: &Data,
    data2: &Data,
    prior1: &PriorParams,
    prior2: &PriorParams,
) -> f64 {
    let mut sum = 0.0;
    for x1 in 0..=(data1.N - data1.n) {
        for x2 in 0..=(data2.N - data2.n) {
            let inner_alpha: f64 = prior1.alpha + (data1.y as f64) + (x1 as f64);
            let inner_beta: f64 = prior1.beta + (data1.N as f64) - (data1.y as f64) - (x1 as f64);
            let inner: &dyn Fn(f64) -> f64 = &|p1| beta_dens(p1, inner_alpha, inner_beta);
            let outer_alpha: f64 = prior2.alpha + (data2.y as f64) + (x2 as f64);
            let outer_beta: f64 = prior2.beta + (data2.N as f64) - (data2.y as f64) - (x2 as f64);
            let outer: &dyn Fn(f64) -> f64 = &|p2| {
                beta_dens(p2, outer_alpha, outer_beta)
                    * trapezoid(inner, p2, 1.0, 1.0E-2, 1000)
                        .expect("inner integral failed to converge")
            };
            let prob = trapezoid(outer, 0.0, 1.0, 1.0E-3, 1000)
                .expect("posterior probability failed to converge");
            if prob >= target {
                sum += betabin_dens(
                    x1,
                    data1.N - data1.n,
                    prior1.alpha + (data1.y as f64),
                    prior1.beta + (data1.n as f64) - (data1.y as f64),
                ) * betabin_dens(
                    x2,
                    data2.N - data2.n,
                    prior2.alpha + (data2.y as f64),
                    prior2.beta + (data2.n as f64) - (data2.y as f64),
                );
            }
        }
    }
    sum
}
