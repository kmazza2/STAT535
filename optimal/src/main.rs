use std::cmp;

fn main() {
    const MAX_SAMPLE: usize = 100;
    let given_params: [GivenParam; 2] = [
        GivenParam {
            alpha: 0.05,
            beta: 0.1,
            p0: 0.1,
            p1: 0.3,
        },
        GivenParam {
            alpha: 0.1,
            beta: 0.1,
            p0: 0.2,
            p1: 0.4,
        },
    ];
    let mut optimal_params: [OptimalParam; 2] = [
        OptimalParam {
            n1: 0,
            n2: 0,
            r1: 0,
            r: 0,
            expectation: f64::INFINITY,
        },
        OptimalParam {
            n1: 0,
            n2: 0,
            r1: 0,
            r: 0,
            expectation: f64::INFINITY,
        },
    ];
    for i in 0..given_params.len() {
        for n in 1..=MAX_SAMPLE {
            // expectation is at least n1, so don't check values of n1 greater
            // than lowest expectation found yet
            for n1 in 1..=cmp::min(n - 1, optimal_params[i].expectation.floor() as usize) {
                let n2 = n - n1;
                for r1 in 0..=n1 {
                    for r in r1..=r1 + n2 {
                        let param = OptimalParam {
                            n1,
                            n2,
                            r1,
                            r,
                            expectation: f64::NAN,
                        };
                        // power decreases as r increases, so can stop
                        // iterating when there isn't enough power
                        if param.power(&given_params[i]) < 1.0 - given_params[i].beta {
                            break;
                        }
                        if param.size(&given_params[i]) < given_params[i].alpha {
                            let param = OptimalParam {
                                expectation: param.expectation(&given_params[i]),
                                ..param
                            };
                            if param.expectation < optimal_params[i].expectation {
                                optimal_params[i] = param;
                                // increasing r increases expectation, so
                                // next value of r will yield larger expectation
                                // than current one
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    println!("First result: {:?}", optimal_params[0]);
    println!("Second result: {:?}", optimal_params[1]);
}

struct GivenParam {
    alpha: f64,
    beta: f64,
    p0: f64,
    p1: f64,
}

#[derive(Debug)]
struct OptimalParam {
    n1: usize,
    n2: usize,
    r1: usize,
    r: usize,
    expectation: f64,
}

impl OptimalParam {
    fn size(&self, given: &GivenParam) -> f64 {
        let n1 = self.n1;
        let n2 = self.n2;
        let r1 = self.r1;
        let r = self.r;
        let p0 = given.p0;
        let mut result: f64 = 0.0;
        for y1 in (r1 + 1)..=n1 {
            for y2 in ((r as i32) - (y1 as i32) + 1)..=(n2 as i32) {
                if y2 >= 0 {
                    let y2 = y2 as usize;
                    result += (nchoosek(n1, y1) as f64)
                        * p0.powf(y1 as f64)
                        * (1.0 - p0).powf((n1 - y1) as f64)
                        * (nchoosek(n2, y2) as f64)
                        * p0.powf(y2 as f64)
                        * (1.0 - p0).powf((n2 - y2) as f64);
                }
            }
        }
        result
    }
    fn power(&self, given: &GivenParam) -> f64 {
        let n1 = self.n1;
        let n2 = self.n2;
        let r1 = self.r1;
        let r = self.r;
        let p1 = given.p1;
        let mut result: f64 = 0.0;
        for y1 in (r1 + 1)..=n1 {
            for y2 in ((r as i32) - (y1 as i32) + 1)..=(n2 as i32) {
                if y2 >= 0 {
                    let y2 = y2 as usize;
                    result += (nchoosek(n1, y1) as f64)
                        * p1.powf(y1 as f64)
                        * (1.0 - p1).powf((n1 - y1) as f64)
                        * (nchoosek(n2, y2) as f64)
                        * p1.powf(y2 as f64)
                        * (1.0 - p1).powf((n2 - y2) as f64);
                }
            }
        }
        result
    }
    fn expectation(&self, given: &GivenParam) -> f64 {
        let n1 = self.n1;
        let n2 = self.n2;
        let r1 = self.r1;
        let r = self.r;
        let p0 = given.p0;
        let mut factor = 0.0;
        for y1 in (r1 + 1)..=r {
            if y1 <= n1 {
                factor += (nchoosek(n1, y1) as f64)
                    * p0.powf(y1 as f64)
                    * (1.0 - p0).powf((n1 - y1) as f64)
            }
        }
        (n1 as f64) + (n2 as f64) * factor
    }
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

#[test]
fn test_log_gamma() {
    assert!((log_gamma(1.0) - 0.0).abs() < 0.001);
    assert!((log_gamma(10.0) - 12.80182748008).abs() < 0.001);
    assert!((log_gamma(100.0) - 359.13420536958).abs() < 0.001);
    assert!((log_gamma(1000.0) - 5905.22042320918).abs() < 0.001);
    assert!((log_gamma(10000.0) - 82099.71749644238).abs() < 0.001);
}

fn nchoosek(n: usize, k: usize) -> usize {
    (log_gamma((n + 1) as f64) - log_gamma((k + 1) as f64) - log_gamma((n - k + 1) as f64))
        .exp()
        .round() as usize
}

#[test]
fn test_nchoosek() {
    assert_eq!(4845_usize, nchoosek(20, 4));
    assert!(494705855605882455_f64 - (nchoosek(87, 17) as f64) < 10_000_000_f64);
    assert_eq!(84_usize, nchoosek(9, 3));
}
