function fh = compare_ClusterOccurrence(DB_acc,cfg)

T = cfg.T;
bands = cfg.bands;
nbands = numel(bands);
fh = cell(1,nbands);

for fi = 1 : nbands
    cfg.which_band = bands{fi};
    [cluster_prob_accurate_ci, cluster_prob_error_ci, store] = calculate_ClusterOccurrence(DB_acc, cfg);


    fh{fi} = figure('Position',[300 300 300 300]);
    nexttile
    plotShaded(T,cluster_prob_accurate_ci, [ 1 0 0 ], '-',.7)
    hold on
    plotShaded(T,cluster_prob_error_ci, [0 0 1 ], '-',.7)
    xlabel('Time [s]')
    ylabel('Cluster [a.u.]')
    box off
    ylim([0 20])
    for s = 1 : store.nSign
        plot(T(store.LabelMatrix  == s),20*ones(1,sum((store.LabelMatrix == s))),'color','k','linewidth',4);
    end
end
end




function [cluster_prob_accurate_ci, cluster_prob_error_ci, store] = calculate_ClusterOccurrence(DB_acc, cfg)
T = cfg.T;

nperms = cfg.nperms;
which_band = cfg.which_band;
[ClusterCentroidBand_error,clust_prob_error] = get_SPCClustersProperties(cfg, DB_acc.Clusters_error);
[ClusterCentroidBand_accurate,clust_prob_accurate] = get_SPCClustersProperties(cfg, DB_acc.Clusters_accurate);

cluster_error = (clust_prob_error(ismember(ClusterCentroidBand_error,which_band),:));
cluster_accurate = (clust_prob_accurate(ismember(ClusterCentroidBand_accurate,which_band),:));


cluster_prob = {};
cluster_prob{1} = cluster_error;
cluster_prob{2} = cluster_accurate;

norm_f = size(cluster_error,1) + size(cluster_accurate,1);

clust_prob_perm = cell(1,2);
for fi = 1:2

    clust_prob_perm{fi} = nan(numel(T),nperms);

    %clust_prob_ROI_perms = nan(nperms,numel(T));
    for perm_i = 1 : nperms
        tmp = [];
        cuts = randi(numel(T),1, size(cluster_prob{fi},1));
        for cut_i = 1 : numel(cuts)
            tmp = [tmp; cluster_prob{fi}(cut_i,[cuts(cut_i) : numel(T) 1 : (cuts(cut_i)-1)])];
        end
        clust_prob_perm{fi}(:,perm_i) = 100*sum(tmp)/norm_f;
    end

end


% get bootstrap and permutation significance

clust_prb_func = @(x) sum(x)/norm_f*100;
cluster_prob_error_ci = mybootci(cluster_prob{1},nperms,clust_prb_func);
cluster_prob_accurate_ci = mybootci(cluster_prob{2},nperms,clust_prb_func);
nerror = size(cluster_prob{1},1);
naccurate = size(cluster_prob{2},1);
cluster_prob_all = [cluster_prob{1}; cluster_prob{2}];

t_orig = nan(1,numel(T));

for ti = 1 : numel(T)
    [~,~,~,stats]  = ttest2(cluster_prob{2}(:,ti),cluster_prob{1}(:,ti));
    t_orig(ti) = stats.tstat;
end

clust_prob_perm = cell(nperms,2);
nperms = 100;
t_perm = nan(nperms, numel(T));
for perm_i = 1 : nperms
    idx_error_perm = randi(nerror+ naccurate,1,nerror);
    idx_accurate_perm = randi(nerror+ naccurate,1,naccurate);
    clust_prob_perm{perm_i,1} = cluster_prob_all(idx_accurate_perm,:);
    clust_prob_perm{perm_i,2} = cluster_prob_all(idx_error_perm,:);

    %[~,~,~,stats]  = ttest2(cluster_prob_all(idx_accurate_perm,:),cluster_prob_all(idx_error_perm,:));
    for ti = 1 : numel(T)
        [~,~,~,stats] = ttest2(clust_prob_perm{perm_i,1}(:,ti),clust_prob_perm{perm_i,2}(:,ti));
        t_perm(perm_i,ti) = stats.tstat;
    end
end

meanPerm = repmat(mean(t_perm,1,'omitnan'), nperms,1);
p_perm =  2 * (1 - tcdf(abs(t_perm), nperms-1)); % get p-values from the zscore, abs to make it 2-tailed

p_orig = (sum(abs(t_perm - meanPerm) >= abs(t_orig - meanPerm), 1)+1) / (nperms+1); %V2 from Ernst, Permutation Methods: A Basis for Exact Inference, 2004
store = getSignifClusters(p_orig, t_orig, p_perm, t_perm, 'zstat');
end