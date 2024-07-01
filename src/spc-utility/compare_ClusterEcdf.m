function fh = compare_ClusterEcdf(DB,cfg)

bands = cfg.bands;
[ClusterCentroidBand_error,ClusterOnOff_error] = get_SPCClustersProperties(cfg, DB.Clusters_error);
[ClusterCentroidBand_accurate,ClusterOnOff_accurate] = get_SPCClustersProperties(cfg, DB.Clusters_accurate);

cluster_error = (ClusterOnOff_error(ismember(ClusterCentroidBand_error,bands),:));
cluster_accurate = (ClusterOnOff_accurate(ismember(ClusterCentroidBand_accurate,bands),:));


[fOn_error, xOn_error] = ecdf(cluster_error(:,1));
[fOff_error, xOff_error] = ecdf(cluster_error(:,2));
[fOn_accurate, xOn_accurate] = ecdf(cluster_accurate(:,1));
[fOff_accurate, xOff_accurate] = ecdf(cluster_accurate(:,2));

[~,idx_sortON] = sort(cluster_error(:,1),'ascend');
cluster_error_sortedbyON = cluster_error(idx_sortON,:);
[~,idx_sortON] = sort(cluster_accurate(:,1),'ascend');
cluster_accurate_sortedbyON = cluster_accurate(idx_sortON,:);


fh = figure('Position',[200 200 400 300]);
plot(xOn_error,fOn_error)
hold on
plot(xOff_error,fOff_error)
plot(xOn_accurate,fOn_accurate)
plot(xOff_accurate,fOff_accurate)
legend({'ON error','OFF error','ON accurate','OFF accurate'},'location','best')
xlim([-2.75 2.5])
box off
xlabel('Time [s]')
ylabel(' cdf')

