function store = comp_clusterPLV(origPLV, PLV_perm, MPC)

if size(origPLV,1) > 1 % if we calculated it for several frequencies
    NUM_PERMS = size(PLV_perm,3);

    zscores_orig = (origPLV - mean(PLV_perm,3)) ./ std(PLV_perm,[],3);
    zscores_perm = bsxfun(@rdivide, bsxfun(@minus, PLV_perm, mean(PLV_perm,3)), std(PLV_perm,[],3));

    p_perm =  2 * (1 - normcdf(abs(zscores_perm), 0, 1)); % get p-values from the zscore, abs to make it 2-tailed
    meanPerm = repmat(mean(zscores_perm,3), 1, 1, NUM_PERMS);

    %p_orig = 2 * (1 - normcdf(abs(zscores_orig), 0, 1)); % get p-values from the zscore, abs to make it 2-tailed

    %         % Alternative way of computing the p-value
    p_orig = (sum(abs(zscores_perm - meanPerm) >= abs(zscores_orig - meanPerm), 3)+1) / (NUM_PERMS+1); %V2 from Ernst, Permutation Methods: A Basis for Exact Inference, 2004

else
    % if I only calculated it for a single frequency
    NUM_PERMS = size(PLV_perm,2);

    PLV_perm = squeeze(PLV_perm);
    zscores_orig = [(origPLV - mean(PLV_perm,2)') ./ std(PLV_perm,[],2)']';
    zscores_perm = bsxfun(@rdivide, bsxfun(@minus, PLV_perm, mean(PLV_perm)), std(PLV_perm));
    meanPerm = repmat(mean(zscores_perm,2), 1,NUM_PERMS);
    p_perm =  2 * (1 - normcdf(abs(zscores_perm), 0, 1)); % get p-values from the zscore, abs to make it 2-tailed

    %p_orig = 2 * (1 - normcdf(abs(zscores_orig), 0, 1)); % get p-values from the zscore, abs to make it 2-tailed

    %         % Alternative way of computing the p-value
    p_orig = (sum(abs(zscores_perm - meanPerm) >= abs(zscores_orig - meanPerm), 2)+1) / (NUM_PERMS+1); %V2 from Ernst, Permutation Methods: A Basis for Exact Inference, 2004
    
    p_orig = p_orig';
     zscores_orig = zscores_orig';
   
    

end


% calculate clusters
store = getSignifClusters(p_orig, zscores_orig, p_perm, zscores_perm, MPC);



end