use rand::Rng;
fn main() {
    let mut rng = rand::thread_rng();
    const SCENARIOS: usize = 4;
    const REPLICATIONS: usize = 10_000;
    const LEVELS: usize = 6;
    const PROPERTIES: usize = 3;
    const SCENARIO_SIZE: usize = REPLICATIONS * LEVELS * PROPERTIES;
    const REPLICATION_SIZE: usize = LEVELS * PROPERTIES;
    const LEVEL_SIZE: usize = PROPERTIES;
    const PATIENTS_PROP: usize = 0;
    const DLTS_PROP: usize = 1;
    const MTD_PROP: usize = 2;
    // idx is a hack I used to store structured usize data in a single-dimensional
    // array (called "data" below); don't do this
    fn idx(scenario: usize, replication: usize, level: usize, property: usize) -> usize {
        return scenario * SCENARIO_SIZE
            + replication * REPLICATION_SIZE
            + level * LEVEL_SIZE
            + property;
    }
    let mut test_cohort = |tox_prob| {
        [
            (rng.gen::<f64>() < tox_prob) as usize,
            (rng.gen::<f64>() < tox_prob) as usize,
            (rng.gen::<f64>() < tox_prob) as usize,
        ]
        .iter()
        .sum::<usize>()
    };
    let mut data = [0; SCENARIOS * REPLICATIONS * LEVELS * PROPERTIES];
    let tox_prob = [
        0.00, 0.30, 0.40, 0.55, 0.60, 0.65, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.00, 0.02, 0.06,
        0.10, 0.20, 0.30, 0.00, 0.20, 0.30, 0.60, 0.70, 0.75,
    ];
    for scenario in 0..SCENARIOS {
        for replication in 0..REPLICATIONS {
            let mut level: usize = 1;
            loop {
                if level == 0 {
                    data[idx(scenario, replication, level, MTD_PROP)] = 1;
                    break;
                }
                match data[idx(scenario, replication, level, PATIENTS_PROP)] {
                    0 => {
                        data[idx(scenario, replication, level, DLTS_PROP)] +=
                            test_cohort(tox_prob[scenario * 6 + level]);
                        data[idx(scenario, replication, level, PATIENTS_PROP)] += 3;
                        match data[idx(scenario, replication, level, DLTS_PROP)] {
                            0 => {
                                if level < 5 {
                                    level += 1;
                                }
                            }
                            1 => {}
                            2 | 3 => {
                                level -= 1;
                            }
                            _ => panic!("impossible DLT size"),
                        }
                    }
                    3 => {
                        data[idx(scenario, replication, level, DLTS_PROP)] +=
                            test_cohort(tox_prob[scenario * 6 + level]);
                        data[idx(scenario, replication, level, PATIENTS_PROP)] += 3;
                        match data[idx(scenario, replication, level, DLTS_PROP)] {
                            0 | 1 => {
                                if level < 5
                                    && data[idx(scenario, replication, level + 1, DLTS_PROP)] < 2
                                {
                                    level += 1;
                                }
                            }
                            2 => {
                                data[idx(scenario, replication, level - 1, MTD_PROP)] = 1;
                                break;
                            }
                            3 | 4 | 5 | 6 => {
                                level -= 1;
                            }
                            _ => panic!("impossible DLT size"),
                        }
                    }
                    6 => match data[idx(scenario, replication, level, DLTS_PROP)] {
                        0 | 1 => {
                            if level == 5
                                || data[idx(scenario, replication, level + 1, DLTS_PROP)] > 1
                            {
                                data[idx(scenario, replication, level, MTD_PROP)] = 1;
                                break;
                            } else {
                                level += 1;
                            }
                        }
                        2 | 3 | 4 | 5 | 6 => {
                            level -= 1;
                        }
                        _ => panic!("impossible DLT size"),
                    },
                    _ => panic!("impossible patient size"),
                }
            }
        }
    }
    for scenario in 0..SCENARIOS {
        println!("Scenario: {}", scenario + 1);
        let mut total_patients = 0;
        let mut total_dlts = 0;
        for replication in 0..REPLICATIONS {
            for level in 0..LEVELS {
                total_patients += data[idx(scenario, replication, level, PATIENTS_PROP)];
                total_dlts += data[idx(scenario, replication, level, DLTS_PROP)];
            }
        }
        let avg_patients = (total_patients as f64) / (REPLICATIONS as f64);
        let avg_dlts = (total_dlts as f64) / (REPLICATIONS as f64);
        println!(" Average patients: {}", avg_patients);
        println!(" Average DLTs: {}", avg_dlts);
        for level in 0..LEVELS {
            println!(" Level {}:", level);
            let mut total_level_mtds = 0;
            let mut total_level_patients = 0;
            for replication in 0..REPLICATIONS {
                total_level_mtds += data[idx(scenario, replication, level, MTD_PROP)];
                total_level_patients += data[idx(scenario, replication, level, PATIENTS_PROP)];
            }
            let percent_level_mtd = 100.0 * (total_level_mtds as f64) / (REPLICATIONS as f64);
            let avg_level_patients = (total_level_patients as f64) / (REPLICATIONS as f64);
            println!(" Percent MTD: {}%", percent_level_mtd);
            println!(" Average patients: {}", avg_level_patients);
        }
    }
}
