function [r,delay, pval, f] = compute_delay_withperm(freqs,phase, nbtsp)
freqs_ = unique(freqs);
phase_mean = arrayfun(@(x) circ_mean(phase(freqs  == x)), freqs_);
if isempty(nbtsp)
    [r,delay, f, pval] = compute_delay(freqs_,phase_mean);
else
    if nargout == 3
        [r,delay] = compute_delay(freqs_,phase_mean);
    else
        [r,delay, f] = compute_delay(freqs_,phase_mean);
    end
    % checvk permutation test
    %phase_groups = arrayfun(@(x) phase(freqs  == x), freqs_, 'uni',false);
    r_random_shuffle = nan(1, nbtsp);

    for btsp = 1 : nbtsp
        %phase_random = cellfun(@(x) x(randi(numel(x))), phase_groups);
        phase_unwrap = unwrap(phase_mean);
        phase_random_shuffle = phase_unwrap(randperm(numel(phase_unwrap)));
        r_random_shuffle(btsp) = corr(freqs_,phase_random_shuffle);
    end
    % (1 + sum(abs(plv_perm1) >= abs(cplv_true)))/(iter_num + 1)
    pval = (1 + sum(abs(r_random_shuffle) >= abs(r)))/(nbtsp + 1);
end



