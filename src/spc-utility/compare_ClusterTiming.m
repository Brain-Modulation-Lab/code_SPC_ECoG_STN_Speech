function [fh, varargout] = compare_ClusterTiming(DB,cfg)

bands = cfg.bands;
[ClusterCentroidBand_error,ClusterOnOff_error] = get_SPCClustersProperties(cfg, DB.Clusters_error);
[ClusterCentroidBand_accurate,ClusterOnOff_accurate] = get_SPCClustersProperties(cfg, DB.Clusters_accurate);

cluster_error = (ClusterOnOff_error(ismember(ClusterCentroidBand_error,bands),:));
cluster_accurate = (ClusterOnOff_accurate(ismember(ClusterCentroidBand_accurate,bands),:));
%ClusterSubjects_error = ClusterSubjects_error(ismember(ClusterCentroidBand_error,bands));
%ClusterSubjects_accurate = ClusterSubjects_accurate(ismember(ClusterCentroidBand_accurate,bands));


[unique_error,idx_error] = getUniquePairs(DB.Clusters_error);
idx_error = idx_error(ismember(ClusterCentroidBand_error,[2 3]));

[unique_accurate,idx_accurate] = getUniquePairs(DB.Clusters_accurate);
idx_accurate = idx_accurate(ismember(ClusterCentroidBand_accurate,[2 3]));

unique_all = unique([unique_error; unique_accurate],'rows');

OnOff_all = nan(size(unique_all,1),2);
Off_all = nan(size(unique_all,1),2);
Duration_all = nan(size(unique_all,1),2);
Duration_all_pre = nan(size(unique_all,1),2);
Duration_all_post = nan(size(unique_all,1),2);

%OnOff_all_std = nan(size(unique_all,1),2);

duration_cluster_error = diff(cluster_error,[],2);
duration_cluster_accurate = diff(cluster_accurate,[],2);



for u = 1 : size(unique_all,1)
    [~,idx] = ismember(unique_all(u,:),unique_error,'rows');
    OnOff_all(u,1) = mean(cluster_error(idx_error == idx ,1),'omitnan');
    Off_all(u,1) = mean(cluster_error(idx_error == idx ,2),'omitnan');
    
    %OnOff_all_std(u,1) = std(meanpoint_cluster_error(idx_error == idx),'omitnan');
    Duration_all(u,1) = mean(duration_cluster_error(idx_error == idx ),'omitnan');
    Duration_all_pre(u,1) = mean(duration_cluster_error(idx_error == idx &  cluster_error(:,1) < 0),'omitnan');
    Duration_all_post(u,1) = mean(duration_cluster_error(idx_error == idx &  cluster_error(:,1) >= 0),'omitnan');
    


    
    [~,idx] = ismember(unique_all(u,:),unique_accurate,'rows');
  
    OnOff_all(u,2) = mean(cluster_accurate(idx_accurate == idx ,1),'omitnan');
    Off_all(u,2) = mean(cluster_accurate(idx_accurate == idx ,2),'omitnan');
    
    Duration_all(u,2) = mean(duration_cluster_accurate(idx_accurate == idx ),'omitnan');
    %OnOff_all_std(u,2) = std(meanpoint_cluster_accurate(idx_accurate == idx),'omitnan');
    Duration_all_pre(u,2) = mean(duration_cluster_accurate(idx_accurate == idx &  cluster_accurate(:,1) < 0),'omitnan');
    Duration_all_post(u,2) = mean(duration_cluster_accurate(idx_accurate == idx &  cluster_accurate(:,1) >= 0),'omitnan');
    
end

store = struct();
store.OnOff_all = OnOff_all;
store.Off_all = Off_all;
store.Pair_info = unique_all;
varargout{1} = store;



fh{1} = figure('Position',[200 200 300 300]);
scatter_jitter_paired({OnOff_all(:,1),OnOff_all(:,2)},{'error','accurate'},.03,1:2,[0 0 1;1 0 0])
[Tstat,pTstat] = compare_stat_scatterplot(OnOff_all(:,1),OnOff_all(:,2),'ttest_dep');
title(sprintf("(t = %1.2f, p = %1.3f)",Tstat,pTstat))
ylim([-2.75 2.5])
ylabel('Time [s]')

fh{2} = figure('Position',[200 200 600 300]);
nexttile(1)
scatter_jitter_paired({Duration_all_pre(:,1),Duration_all_pre(:,2)},{'error','accurate'},.03,1:2,[0 0 1;1 0 0])
[Tstat,pTstat] = compare_stat_scatterplot(Duration_all_pre(:,1),Duration_all_pre(:,2),'ttest_dep');
title(sprintf("(t = %1.2f, p = %1.3f)",Tstat,pTstat))
ylim([-2.75 2.5])
ylabel('Duration [s]')
ylim([0 3.5])
nexttile(2)
scatter_jitter_paired({Duration_all_post(:,1),Duration_all_post(:,2)},{'error','accurate'},.03,1:2,[0 0 1;1 0 0])
[Tstat,pTstat] = compare_stat_scatterplot(Duration_all_post(:,1),Duration_all_post(:,2),'ttest_dep');
title(sprintf("(t = %1.2f, p = %1.3f)",Tstat,pTstat))
ylim([-2.75 2.5])
ylabel('Duration [s]')
ylim([0 3.5])