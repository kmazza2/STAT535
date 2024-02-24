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
            size: f64::NAN,
            power: f64::NAN,
            expectation: f64::NAN,
        },
        OptimalParam {
            n1: 0,
            n2: 0,
            r1: 0,
            r: 0,
            size: f64::NAN,
            power: f64::NAN,
            expectation: f64::NAN,
        },
    ];
    for i in 0..given_params.len() {
        let mut found: bool = false;
        'given_param_loop: for n in 1..=MAX_SAMPLE {
            for n1 in 1..=(n - 1) {
                let n2 = n - n1;
                for r1 in 0..=n1 {
                    for r in r1..=(r1 + n2) {
                        let param = OptimalParam {
                            n1,
                            n2,
                            r1,
                            r,
                            size: f64::NAN,
                            power: f64::NAN,
                            expectation: f64::NAN,
                        };
                        let param = OptimalParam {
                            size: param.size(&given_params[i]),
                            power: param.power(&given_params[i]),
                            ..param
                        };
                        if param.power >= 1.0 - given_params[i].beta
                            && param.size < given_params[i].alpha
                        {
                            let param = OptimalParam {
                                expectation: param.expectation(&given_params[i]),
                                ..param
                            };
                            optimal_params[i] = param;
                            found = true;
                        }
                        if found {
                            break 'given_param_loop;
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
    size: f64,
    power: f64,
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
    const COEF: [f64; 6] = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179E-2,
        -0.5395239384953E-5,
    ];
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
