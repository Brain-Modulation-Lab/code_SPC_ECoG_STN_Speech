function store = run_Cluster2DMaps(data,xQ,yQ, cfg, nperms)


%baselne = cfg.baseline;
stat_test = cfg.stat_test; % can be 'paired-test-conds', 'baseline-test'

%data_perm = nan(size(data));
%temp = cat(3,PPC_mat_all_adj{:});
tscore_orig = nan(numel(yQ), numel(xQ));
tscore_perm = nan(nperms, numel(yQ), numel(xQ));

disp(' Prepare baseline for true PPCz map')

switch stat_test
    case  'paired-test-conds'
        assert(iscell(data),'Data must be a cell array!')
        disp(' Calculating true tscore...')

        contrast_data = data{1} - data{2};

        for fi = 1 : numel(yQ)
            fprintf(' Computing true t-score in %d - %d freq bin \n',fi, numel(yQ))
            for ti = 1 : numel(xQ)
                [~,~,~,tscore_orig_stat] = ttest(squeeze(contrast_data(fi,ti,:))); % one-sample t test vs 0
                tscore_orig(fi,ti) = tscore_orig_stat.tstat;
            end
        end
        disp(' Done computation true t-score...')

        disp(' Starting permuted t-score...')



        for perm_i = 1 : nperms
            fprintf(' Computing permuted t-score in %d - %d permutation \n',perm_i, nperms)
            
            arrangement = 2 * randi([0, 1], size(contrast_data,3),1 ) - 1; % random +1 and -1 trick
            
            %data_perm = data{1}(:,:,randperm(size(contrast_data,3),size(contrast_data,3))) - data{2}(:,:,randperm(size(contrast_data,3),size(contrast_data,3)));
            %splitIndx = randi(numel(xQ),1,size(data,3));

            %baseline_dataperm = squeeze(mean(data_perm(:,xQ <= baseline(2) & xQ >= baseline(1),:),2,'omitnan'));
            for fi = 1 : numel(yQ)
                for ti = 1 : numel(xQ)
                    [~,~,~,tscore_perm_stat] = ttest(squeeze(contrast_data(fi,ti,:)).*arrangement);
                    tscore_perm(perm_i,fi,ti) = tscore_perm_stat.tstat;
                end
            end

        end

        disp(' Done computation permuted t-score...')

        disp(' Calculating p-values:true and perm')

        tscore_perm = permute(tscore_perm, [2 3 1]);
        p_perm =  2 * (1 - tcdf(abs(tscore_perm), size(contrast_data,3)-1)); % get p-values from the zscore, abs to make it 2-tailed
        meanPerm = repmat(mean(tscore_perm,3,'omitnan'), 1, 1, nperms);
        p_orig = (sum(abs(tscore_perm - meanPerm) >= abs(tscore_orig - meanPerm), 3)+1) / (nperms+1); %V2 from Ernst, Permutation Methods: A Basis for Exact Inference, 2004



        disp('Calculating clusters')
        store = getSignifClusters(p_orig, tscore_orig, p_perm, tscore_perm, 'zstat',0.05,0.01);
        store.TscoreMap = tscore_orig;

    case 'baseline-test'
        baseline = cfg.Baseline;
        baseline_data = squeeze(mean(data(:,xQ <= baseline(2) & xQ >= baseline(1),:),2,'omitnan'));


        for fi = 1 : numel(yQ)
            fprintf(' Computing true t-score in %d - %d freq bin \n',fi, numel(yQ))
            for ti = 1 : numel(xQ)
                [~,~,~,tscore_orig_stat] = ttest(baseline_data(fi,:)',squeeze(data(fi,ti,:))); % w.r.t baseline test
                tscore_orig(fi,ti) = tscore_orig_stat.tstat;
            end
        end
        disp(' Done computation true t-score...')

        disp(' Starting permuted t-score...')


        for perm_i = 1 : nperms
            fprintf(' Computing permuted t-score in %d - %d permutation \n',perm_i, nperms)

            data_perm = nan(numel(yQ),numel(xQ),size(data,3));
            splitIndx = randi(numel(xQ),1,size(data,3));
            for pair_i = 1 : size(data,3)
                data_perm(:,:,pair_i) = data(:, [(splitIndx(pair_i)+1) : numel(xQ)  1 : splitIndx(pair_i)],pair_i );
            end
            baseline_dataperm = squeeze(mean(data_perm(:,xQ <= baseline(2) & xQ >= baseline(1),:),2,'omitnan'));
            for fi = 1 : numel(yQ)
                for ti = 1 : numel(xQ)
                    [~,~,~,tscore_perm_stat] = ttest(baseline_dataperm(fi,:)',squeeze(data_perm(fi,ti,:)));
                    tscore_perm(perm_i,fi,ti) = tscore_perm_stat.tstat;
                end
            end

        end

        tscore_orig = -tscore_orig;
        tscore_perm = -tscore_perm;
        disp(' Done computation permuted t-score...')

        disp(' Calculating p-values:true and perm')

        tscore_perm = permute(tscore_perm, [2 3 1]);
        p_perm =  2 * (1 - tcdf(abs(tscore_perm), size(data,3)-1)); % get p-values from the zscore, abs to make it 2-tailed
        meanPerm = repmat(mean(tscore_perm,3,'omitnan'), 1, 1, nperms);
        p_orig = (sum(abs(tscore_perm - meanPerm) >= abs(tscore_orig - meanPerm), 3)+1) / (nperms+1); %V2 from Ernst, Permutation Methods: A Basis for Exact Inference, 2004



        disp('Calculating clusters')
        store = getSignifClusters(p_orig, tscore_orig, p_perm, tscore_perm, 'zstat');
        store.TscoreMap = tscore_orig;

end


% if baseline == 0 || isempty(baseline)
%     baseline_data = 0;
% else
% end




