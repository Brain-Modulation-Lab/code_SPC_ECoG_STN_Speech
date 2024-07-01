function store = getSignifClusters(p_sig, zscores, p_perm, zscores_perm, MPC, varargin)
% GETSIGNIFCLUSTERS Find significant clusters (in original data) based on
% clusters obtained by a permutation procedure
% (based on Maris & Oostenveld (2007) J Neurosci Meth 164:177-190)
%
% INPUTS
% psig         - significance matrix (FREQ x TIME)
% zscores      - zscores to divide into positive and negative sig. clusters
% p_perm       - matrix of p-values derived from permutations
% zscores_perm - zscores derived from permutation to divide into positive and
%                negative sig. clusters (NPERM x FREQ x TIME) where NPERM
%                is the number of permutation numPerms
% MPC = 'zstat' or 'csize' or 'noMCP'

% VARARGIN (only if MPC != noMPC)
% THRESH_SUPRACLUSTER - significance threshold to determine supra-threshold
%                          clusters (default 0.05)
% ALPHA - significance threshold, i.e. percentage of the largest
%         sums of suprathreshold cluster z-scores (or cluster sizes) from
%         the permutation distribution that may be larger than the sum from
%         the original unpermuted data
%
% OUTPUTS

% clusLabelSign_Z_Stat = cluslabels (matrix form)
% clusPval_Z_Stat   - P-values of all possible clusters based on summed z-scores
% clusZval_Z_Stat   - z-values of all possible clusters based on summed z-scores

if nargin < 6
    THRESH_SUPRACLUSTER = .05;
    ALPHA = .05;
elseif nargin < 7
    ALPHA = .05;
    THRESH_SUPRACLUSTER = varargin{1}'
elseif nargin == 7
    THRESH_SUPRACLUSTER = varargin{1};
    ALPHA = varargin{2};
end

numPerms = size(p_perm, 1);


store = [];

% get all pre-cluster thresholds in the orignal sample
[clusLabel, numClus] = getPosAndNegClusters(p_sig, zscores, THRESH_SUPRACLUSTER);

switch MPC
    case 'zstat'
        
        clus  = zeros(1, numClus);
        % for each supra-threshold cluster, sum up the z-scores or the number
        % of pixels in this cluster
        for c = 1:numClus
            clus(c)   = sum(abs(zscores(clusLabel == c))); % abs() for two-tailed testing
        end
        permDist_maxSum   = zeros(1, numPerms);
        %% Get cluster sums (for z-scores and pixels) for all permutations
        
        for i = 1:size(p_perm, 1)
            if mod(i, 500) == 0; fprintf(['   ' num2str(i)]);  end
            [clusLabel_perm, numClus_perm] = getPosAndNegClusters(squeeze(p_perm(i,:,:)), squeeze(zscores_perm(i,:,:)), THRESH_SUPRACLUSTER);
            
            permClus   = zeros(1, numClus_perm);
            % get the summed cluster stats
            for c = 1:numClus_perm
                permClus(c)   = sum(abs(zscores_perm(i, clusLabel_perm == c))); % abs() for two-tailed testing
            end
            
            % store only the sum of the largest cluster (based on z-scores or size) for this iteration
            if numClus_perm>0 % if significant clusters were present, otherwise leave the field at 0
                permDist_maxSum(i)   = max(permClus);  % minimum because z-values are all negative as norminv from p
            end
        end
        
        %fprintf('\n');
        % now compare the obtained clusters to the permutation clusters
        clusPval   = nan(numClus,1);
        
        % For each original supra-threshold cluster check how many of the
        % maximum sums of the permutation distribution exceeds the sum of
        % z-scores/pixels of the current cluster
        for c = 1:numClus
            valExtreme         = sum(permDist_maxSum >= clus(c));
            clusPval(c) = valExtreme / numPerms;
        end
        
        nSign = sum(clusPval < ALPHA);
        
        if nSign > 0
            clusPos = ismember(clusLabel, find(clusPval < ALPHA));
            LabelMatrix = clusPos.*clusLabel;
            clust_id = unique(LabelMatrix(LabelMatrix~= 0));
            for clus_i = 1 : nSign
                [LabelMatrix(LabelMatrix == clust_id(clus_i))] = deal(clus_i);
            end
            store.LabelMatrix = LabelMatrix;
            store.Zstat = clus(clusPval < ALPHA);
            store.Pval = clusPval(clusPval < ALPHA);
            
        else
            store.LabelMatrix = nan;
            store.Zstat = nan;
            store.Pval = nan;
        end
        
        % Tag all the significant clusters derived from the "sum of pixels"
        % (i.e. cluster size) method as 1 in a boolean array
        % clusPos_Z_Stat = (FREQ x TIME) or (1 x TIME/FREQ) in 1-d case
        %     clusPval_Z_Stat % print the p-values to the console
        %     clusPval_clusSize
        
        store.nSign = nSign;
        
    case 'csize'
        
        clus= zeros(1, numClus);
        
        for c = 1:numClus
            clus(c) = sum(clusLabel(:) == c);
        end
        
        permDist_maxSum = zeros(1, numPerms);
        
        for i = 1:size(p_perm, 1)
            if mod(i, 500) == 0; fprintf(['   ' num2str(i)]);  end
            [clusLabel_perm, numClus_perm] = getPosAndNegClusters(squeeze(p_perm(i,:,:)), squeeze(zscores_perm(i,:,:)), THRESH_SUPRACLUSTER);
            
            permClus = zeros(1, numClus_perm);
            % get the summed cluster stats
            for c = 1:numClus_perm
                permClus(c) = sum(clusLabel_perm(:) == c);
            end
            
            % store only the sum of the largest cluster (based on z-scores or size) for this iteration
            if numClus_perm>0 % if significant clusters were present, otherwise leave the field at 0
                permDist_maxSum(i) = max(permClus);
            end
        end
        
        %fprintf('\n');
        % now compare the obtained clusters to the permutation clusters
        clusPval = nan(numClus,1);
        
        % For each original supra-threshold cluster check how many of the
        % maximum sums of the permutation distribution exceeds the sum of
        % z-scores/pixels of the current cluster
        for c = 1:numClus
            valExtreme           = sum(permDist_maxSum>= clusc);
            clusPval(c) = valExtreme / numPerms;
        end
        
        % Tag all the significant clusters derived from the "sum of z-scores"
        % method as 1 in a boolean array
        % clusPos_Z_Stat = (FREQ x TIME) or (1 x TIME/FREQ) in 1-d case
        % Tag all the significant clusters derived from the "sum of pixels"
        % (i.e. cluster size) method as 1 in a boolean array
        % clusPos_Z_Stat = (FREQ x TIME) or (1 x TIME/FREQ) in 1-d case
        nSign = sum(clusPval < ALPHA);
        
        if nSign > 0
            clusPos = ismember(clusLabel, find(clusPval < ALPHA));
            LabelMatrix = clusPos.*clusLabel;
            clust_id = unique(LabelMatrix(LabelMatrix~= 0));
            for clus_i = 1 : nsign
                [LabelMatrix(LabelMatrix == clust_id(clus_i))] = deal(clus_i);
            end
            store.LabelMatrix = LabelMatrix;
            store.Zstat = clus(clusPval < ALPHA);
            store.Pval = clusPval(clusPval < ALPHA);
        else
            store.LabelMatrix = nan;
            store.Zstat = nan;
            store.Pval = nan;
        end
        
        store.nSign = nSign;
        
        % Tag all the significant clusters derived from the "sum of pixels"
        % (i.e. cluster size) method as 1 in a boolean array
        % clusPos_Z_Stat = (FREQ x TIME) or (1 x TIME/FREQ) in 1-d case
        %     clusPval_Z_Stat % print the p-values to the console
        %     clusPval_clusSize
        
        
    case 'noMPC'
        store.LabelMatrix = p_orig < 0.05;
        store.Zstat = z_scores;
        store.Pval = p_orig;
        store.nSign = nan;
end
store.MPC = MPC;
end

function [clusLabel, numClus] = getPosAndNegClusters(p_sig, zscores, THRESH_SUPRACLUSTER)
% GETPOSANDNEGCLUSTERS Find significant clusters but separate them into
%                      positive and negative clusters
% INPUTS
% psig         - significance matrix (FREQ x TIME)
% zscores      - zscores to divide into positive and negative sig. clusters
% THRESH_SUPRACLUSTER - significance threshold to determine supra-threshold
%                          clusters (default 0.05)
%
% OUTPUTS
% clusLabel - Cluster labels (each cluster is labeled as int number
%             starting at 1)
% numClus   - Number of clusters

% threshold the data
p_subThreshold = p_sig < THRESH_SUPRACLUSTER;
zscores_thresh = zscores.*p_subThreshold;
% find the negative and positive clusters separately
[clusNegative, numClus1] = bwlabeln(zscores_thresh < 0);
[clusPositive, numClus2] = bwlabeln(zscores_thresh > 0);
clus_tmp = clusPositive + numClus1; % increase labels by the num of neg clusters
clus_tmp(clusPositive == 0) = 0;   % and make sure to set old 0s back to 0s
clusLabel = clusNegative + clus_tmp;
numClus = numClus1 + numClus2;
end

