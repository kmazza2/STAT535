tic;
parpool;
rng('default');
scenarios = 4;
num_levels = 5;
replications = 10000;
max_patients = 30;
% tox_prob(scenario, level)
tox_prob = [ 0.3 0.4 0.55 0.6 0.65 ;
    0.1 0.2 0.3 0.4 0.5 ;
    0.02 0.06 0.1 0.2 0.3 ;
    0.2 0.3 0.6 0.7 0.75 ];
skeleton = [0.1 0.2 0.3 0.4 0.5];
target = 0.3;
cohort_size = 3;
for scenario = 1:scenarios
    parfor replication = 1:replications
        patients = zeros(num_levels, 1);
        dlts = zeros(num_levels, 1);
        level = 1;
        trial_failed = false;
        while sum(patients) < max_patients
            dlts(level) = dlts(level) + binornd(cohort_size, tox_prob(scenario, level));
            patients(level) = patients(level) + cohort_size;
            means = post_means(num_levels, skeleton, patients, dlts);
            estimate = level_closest_to_target(num_levels, means, target);
            level = next_level(estimate, level);
            if lowest_level_too_toxic(num_levels, skeleton, patients, dlts, target)
                results(scenario, replication).conclusive = false;
                results(scenario, replication).patients = patients;
                results(scenario, replication).dlts = dlts;
                results(scenario, replication).estimate = -1;
                trial_failed = true;
                break;
            end
        end
        if ~trial_failed
            means = post_means(num_levels, skeleton, patients, dlts);
            estimate = level_closest_to_target(num_levels, means, target);
            results(scenario, replication).conclusive = true;
            results(scenario, replication).patients = patients;
            results(scenario, replication).dlts = dlts;
            results(scenario, replication).estimate = estimate;
        end
    end
end
toc;
for scenario = 1:scenarios
    fprintf('Scenario %u:\n', scenario);
    fprintf(' Average total toxicity: %g\n', ...
        avg_total_tox(scenario, replications, num_levels, results));
    fprintf(' Average number of patients: %g\n', ...
        avg_num_patients(scenario, replications, num_levels, results));
    for level = 1:num_levels
        fprintf(' Level %u:\n', level);
        fprintf(' MTD Percent: %g\n', ...
            percent_mtd(scenario, replications, level, results));
        fprintf(' Avg Patients: %g\n', ...
            avg_patients(scenario, replications, level, results));
    end
end
delete(gcp);
function result = cond_likeli(num_levels, skeleton, patients, dlts, a)
result = 1;
for level = 1:num_levels
    result = result .* ...
        ((skeleton(level).^exp(a)).^dlts(level)) .* ...
        ((1 - skeleton(level).^exp(a)).^(patients(level) - dlts(level)));
end
end
function result = marginal_likeli(num_levels, skeleton, patients, dlts)
result = integral( ...
    @(a) ...
    cond_likeli(num_levels, skeleton, patients, dlts, a) .* ...
    normpdf(a, 0, 1.34), ...
    -Inf, Inf ...
    );
end
function result = joint_dens(num_levels, skeleton, patients, dlts, a)
result = cond_likeli(num_levels, skeleton, patients, dlts, a) .* ...
    normpdf(a, 0, 1.34);
end
function result = expec_integrand(num_levels, skeleton, patients, dlts, level, a)
result = (skeleton(level).^exp(a)) .* ...
    joint_dens(num_levels, skeleton, patients, dlts, a) ./ ...
    marginal_likeli(num_levels, skeleton, patients, dlts);
end
function mean = post_mean(num_levels, skeleton, patients, dlts, level)
mean = integral( ...
    @(a) expec_integrand(num_levels, skeleton, patients, dlts, level, a), ...
    -Inf, Inf);
end
function means = post_means(num_levels, skeleton, patients, dlts)
means = zeros(num_levels, 1);
for level = 1:num_levels
    means(level) = post_mean(num_levels, skeleton, patients, dlts, level);
end
end
function result = level_closest_to_target(num_levels, means, target)
result = 1;
for level = 2:num_levels
    if abs(means(level) - target) < abs(means(result) - target)
        result = level;
    end
end
end
function new_level = next_level(estimate, level)
if estimate > level && level < 5
    new_level = level + 1;
elseif estimate < level && level > 1
    new_level = level - 1;
else
    new_level = level;
end
end
function result = post_a(num_levels, skeleton, patients, dlts, a)
result = joint_dens(num_levels, skeleton, patients, dlts, a) ./ ...
    marginal_likeli(num_levels, skeleton, patients, dlts);
end
function result = lowest_level_too_toxic(num_levels, skeleton, patients, dlts, target)
if integral(@(a) post_a(num_levels, skeleton, patients, dlts, a), ...
        -Inf, log(log(target)/log(skeleton(1)))) > 0.9
    result = true;
else
    result = false;
end
end
function result = avg_total_tox(scenario, replications, num_levels, results)
total_tox = 0;
for replication = 1:replications
    for level = 1:num_levels
        total_tox = ...
            total_tox + results(scenario, replication).dlts(level);
    end
end
result = total_tox / replications;
end
function result = avg_num_patients(scenario, replications, num_levels, results)
total_patients = 0;
for replication = 1:replications
    for level = 1:num_levels
        total_patients = ...
            total_patients + results(scenario, replication).patients(level);
    end
end
result = total_patients / replications;
end
function result = percent_mtd(scenario, replications, level, results)
total_mtd = 0;
for replication = 1:replications
    if results(scenario, replication).estimate == level
        total_mtd = total_mtd + 1;
    end
end
result = 100.0 * total_mtd / replications;
end
function result = avg_patients(scenario, replications, level, results)
total_patients = 0;
for replication = 1:replications
    total_patients = total_patients + ...
        results(scenario, replication).patients(level);
end
result = total_patients / replications;
end