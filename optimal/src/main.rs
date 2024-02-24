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
    println!("First result: {:?}", optimal_params[1]);
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
    fn power(&self, given: &GivenParam) -> f64 {
        unimplemented!()
    }
    fn size(&self, given: &GivenParam) -> f64 {
        unimplemented!()
    }
    fn expectation(&self, given: &GivenParam) -> f64 {
        unimplemented!()
    }
}

fn log_gamma(x: f64) -> f64 {
    const COEF: [f64; 6] = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179E-2,
        -0.5395239384953E-5];
    const STP: f64 = 2.5066282746310005;
    let mut y = x;
    let tmp = x + 5.5;
    let tmp = (x + 0.5) * tmp.ln() - tmp;
    let mut ser = 1.000000000190015;
    for j in 0..6 {
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
