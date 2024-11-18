function [corrected_p_tscore, tscore_orig] = permuttest2_unbalanced_tmax(x, nPerms, labels)

n_groups = numel(x);
pairwise_mean_diffs = nan(n_groups, n_groups);
tscore_orig = nan(n_groups, n_groups);

perm_tscore_max = nan(nPerms, 1);  % To store max t-scores for each permutation

% Compute actual pairwise mean differences (observed data)
for i = 1:n_groups
    for j = i+1:n_groups
        % Calculate mean difference between group i and group j
        pairwise_mean_diffs(i, j) = mean(x{i}) - mean(x{j});
        [~, ~, ~, tscore_orig_stat] = ttest2(x{i}, x{j}); % w.r.t baseline test
        tscore_orig(i, j) = tscore_orig_stat.tstat;
    end
end

% Permutation test
for perm = 1:nPerms
    perm_mean_diffs = nan(n_groups, n_groups);  % Reset permuted mean differences
    perm_tscore_diffs = nan(n_groups, n_groups);  % Reset permuted t-scores

    % Permute the group labels and compute permuted mean differences
    for i = 1:n_groups
        for j = i+1:n_groups
            % Concatenate the two groups
            combined = [x{i}; x{j}];
            n1 = length(x{i});
            n2 = length(x{j});

            % Shuffle combined data and split into new permuted groups
            perm_combined = combined(randperm(length(combined)));
            perm_group1 = perm_combined(1:n1);
            perm_group2 = perm_combined(n1+1:end);

            % Compute permuted mean difference
            perm_mean_diffs(i, j) = mean(perm_group1) - mean(perm_group2);
            [~, ~, ~, tscore_perm_stat] = ttest2(perm_group1, perm_group2); % w.r.t baseline test
            perm_tscore_diffs(i, j) = tscore_perm_stat.tstat;
        end
    end

    % Store the maximum absolute mean difference and t-score from this permutation
    perm_tscore_max(perm) = max(abs(perm_tscore_diffs(:)));
end

fprintf('-------------------------------------------------- \n')

% Print results
for i = 1:n_groups
    for j = i+1:n_groups
        % P-value for T-max based on t-scores
        observed_diff_tscore = tscore_orig(i, j);
        pdabs_tscore = abs(perm_tscore_max);
        corrected_p_tscore = (sum(abs(observed_diff_tscore) <= pdabs_tscore) + 1) / (nPerms + 1);
        fprintf("P-val (T-max) %s vs %s: obs-T = %1.4f [t], p = %1.3f \n", labels{i}, labels{j}, observed_diff_tscore, corrected_p_tscore);
    end
end
fprintf('-------------------------------------------------- \n')
