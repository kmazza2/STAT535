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
