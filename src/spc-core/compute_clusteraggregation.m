function [cluster_aggr, clust_prob_perm] = compute_clusteraggregation(cfg, clust_prob)

nperms = cfg.nperms;
T = cfg.T;
% find temporal cluster trndency
clust_prob_perm = nan(numel(T),nperms);
%clust_prob_ROI_perms = nan(nperms,numel(T));
for perm_i = 1 : nperms
    tmp = [];
    cuts = randi(numel(T),1,size(clust_prob,1));
    for cut_i = 1 : numel(cuts)
        tmp = [tmp; clust_prob(cut_i,[cuts(cut_i) : numel(T) 1 : (cuts(cut_i)-1)])];
    end
    clust_prob_perm(:,perm_i) = mean(tmp);
end

%
clust_prob_lthr = prctile(clust_prob_perm,5,2)';
clust_prob_uthr = prctile(clust_prob_perm,95,2)';

cluster_aggr = zeros(1,size(clust_prob,2));
cluster_aggr(bsxfun(@lt,mean(clust_prob),clust_prob_lthr) ) = -1;
cluster_aggr(bsxfun(@gt,mean(clust_prob),clust_prob_uthr)) = 1;
