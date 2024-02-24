warning('off', 'all');

MAX_SAMPLE = 100;

given_params = struct( ...
    'alpha', {0.05, 0.1}, ...
    'beta' , {0.1 , 0.1}, ...
    'p0'   , {0.1 , 0.2}, ...
    'p1'   , {0.3 , 0.4}  ...
);

optimal_params = struct( ...
    'n1',          {NaN, NaN}, ...
    'n2',          {NaN, NaN}, ...
    'r1',          {NaN, NaN}, ...
    'r',           {NaN, NaN}, ...
    'expectation', {Inf, Inf}  ...
);

for i = 1 : length(given_params)
    bar = waitbar(0, sprintf('Progress on param %d', i));
    for n = 1 : MAX_SAMPLE
        waitbar(n / MAX_SAMPLE, bar, sprintf('Progress on param %d: %.1f%%', i, 100 * n / MAX_SAMPLE));
        % expectation is at least n1, so don't check values of n1 greater
        % than lowest expectation found yet
        for n1 = 1 : min(n - 1, floor(optimal_params(i).expectation))
            n2 = n - n1;
            for r1 = 0 : n1
                for r = r1 : r1 + n2
                    params = struct( ...
                        'n1'         , n1, ...
                        'n2'         , n2, ...
                        'r1'         , r1, ...
                        'r'          , r , ...
                        'expectation', NaN ...
                    );
                    % power decreases as r increases, so can stop iterating 
                    % when there isn't enough power
                    if power(params, given_params(i)) < 1 - given_params(i).beta
                        break;
                    end
                    if size(params, given_params(i)) < given_params(i).alpha
                        params.expectation = expectation(params, given_params(i));
                        if params.expectation < optimal_params(i).expectation
                            optimal_params(i) = params;
                            % increasing r increases expectation, so next
                            % value of r will yield larger expectation than
                            % current one
                            break;
                        end
                    end
                end
            end
        end
    end
    close(bar);
end

optimal_params(1)
optimal_params(2)

% MATLAB uses a slow implementation of nchoosek, so I use a different one
function result = mynchoosek(n, k)
    result = round(exp(log(gamma(n+1))-log(gamma(k+1))-log(gamma(n-k+1))));
end

function result = size(params, given_params)
    n1 = params.n1;
    n2 = params.n2;
    r1 = params.r1;
    r = params.r;
    p0 = given_params.p0;
    result = 0;
    for y1 = r1 + 1 : n1
        for y2 = r - y1 + 1 : n2
            if y1 >= 0 && y2 >= 0
                result = result + ...
                    mynchoosek(n1, y1) * p0^y1 * (1 - p0)^(n1 - y1) * ...
                    mynchoosek(n2, y2) * p0^y2 * (1 - p0)^(n2 - y2);
            end
        end
    end
end

function result = power(params, given_params)
    n1 = params.n1;
    n2 = params.n2;
    r1 = params.r1;
    r = params.r;
    p1 = given_params.p1;
    result = 0;
    for y1 = r1 + 1 : n1
        for y2 = r - y1 + 1 : n2
            if y1 >= 0 && y2 >= 0
                result = result + ...
                    mynchoosek(n1, y1) * p1^y1 * (1 - p1)^(n1 - y1) * ...
                    mynchoosek(n2, y2) * p1^y2 * (1 - p1)^(n2 - y2);
            end
        end
    end
end

function result = expectation(params, given_params)
    n1 = params.n1;
    n2 = params.n2;
    r1 = params.r1;
    r = params.r;
    p0 = given_params.p0;
    factor = 0;
    for y1 = r1 + 1 : r
        if y1 <= n1
            factor = factor + mynchoosek(n1, y1) * p0^y1 * (1 - p0)^(n1 - y1);
        end
    end
    result = n1 + n2 * factor;
end