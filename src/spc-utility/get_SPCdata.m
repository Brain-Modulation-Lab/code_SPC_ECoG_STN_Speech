function DB = get_SPCdata(SUBJECTS,toDo, PATHS,FILENAME)

% toDo can be:
% 'all' : all patients 
% [idx] : patients idx

nSUBJECTS = numel(SUBJECTS);

if strcmpi(toDo,'all')
    toDo = 1 : nSUBJECTS;
end

if ~exist('FILENAME','var')
    FILENAME = 'ClustersPLV.mat';
end

nSign_pairs_all = [];
nPairs_all = [];
nSign_pairs_perm_all = {};
Clusters_all = [];
PPC_mat_all = [];
PPCz_mat_all = [];
PPCmedianperm_mat_all = [];
PPCmeanperm_mat_all = [];
PPCstdperm_mat_all = [];
PPCz_condmat_all = [];
PPC_condmat_all = [];
PPCz_condmat_values_all = [];

ES_mat_all = [];
Phase_mat_all = [];
PairsLocation_MNI_all = [];
PLVTimeEvts_all = [];
list_subj = [];
for subj_i = toDo %
    SUBJECT = SUBJECTS{subj_i};
    fprintf('Pooling Clusters  for Subject %s %d - %d \n', SUBJECT, subj_i, nSUBJECTS)

    PATH_clusters = fullfile(PATHS.DataDir,SUBJECT,FILENAME);

    if isfile(PATH_clusters)
        fprintf(' Stacking clusters results in %s \n', PATH_clusters)
        load(PATH_clusters,'Clusters', 'nSign_pairs','nSign_perm_pairs','n_pairs','PairsLocation_MNI', 'PPC_mat','PPCmedianperm_mat','PPCmeanperm_mat','PPCstdperm_mat','PPCz_mat','ES_mat','PLVTimeEvts','Phase_mat','PPCz_condmat','PPC_condmat','PPCz_condmat_values');
        
        
        % put info about subjects id
        if exist('PairsLocation_MNI','var')
            PairsLocation_MNI.subj_id = repmat(SUBJECT, height(PairsLocation_MNI),1);
            PairsLocation_MNI_all = [PairsLocation_MNI_all; PairsLocation_MNI];
        end
        % stack Clusters
        if exist('Clusters','var') && nSign_pairs > 0
            [Clusters.subj_id] = deal(SUBJECT);
            Clusters_all = [Clusters_all Clusters(~isnan([Clusters.S_typeFRmod]))]; % eliminate nan firemod
        end
        if exist('PPC_mat','var')
            PPC_mat_all = [PPC_mat_all  PPC_mat];
        end
        if exist('PPCz_mat','var')
            PPCz_mat_all = [PPCz_mat_all  PPCz_mat];
        end

        if exist('PPCz_condmat','var') % handle when we have two conditions
            PPCz_condmat_all = [PPCz_condmat_all  PPCz_condmat];
        end
        if exist('PPC_condmat','var') % handle when we have two conditions
            PPC_condmat_all = [PPC_condmat_all  PPC_condmat];
        end
        if exist('PPCz_condmat_values','var') % handle when we have two conditions
            PPCz_condmat_values_all = [PPCz_condmat_values_all  PPCz_condmat_values];
        end
        if exist('PPCmedianperm_mat','var')
            PPCmedianperm_mat_all = [PPCmedianperm_mat_all  PPCmedianperm_mat];
            PPCmeanperm_mat_all = [PPCmeanperm_mat_all  PPCmeanperm_mat];
            PPCstdperm_mat_all = [PPCstdperm_mat_all  PPCstdperm_mat];
        end
        if exist('ES_mat','var')
            ES_mat_all = [ES_mat_all  ES_mat];
        end
        if exist('Phase_mat','var')
            Phase_mat_all = [Phase_mat_all Phase_mat];
        end
        if exist('PLVTimeEvts','var')
            PLVTimeEvts_all = [PLVTimeEvts_all PLVTimeEvts];
        end
        % grab information about significant pairs and n_pairs
        if exist('nSign_pairs','var')
            nSign_pairs_all = [ nSign_pairs_all nSign_pairs];
        end
        if exist('nSign_perm_pairs','var')
            nSign_pairs_perm_all{end+1} = nSign_perm_pairs;
        end        
        if exist('n_pairs','var')
            nPairs_all = [nPairs_all n_pairs];
        end
        list_subj = [list_subj subj_i];
        fprintf(' Completed Stacking clusters results in %s \n', PATH_clusters)
        % clear variables for consistency
        clear Clusters nSign_pairs nSign_perm_pairs n_pairs PairsLocation_MNI PPC_mat PPCz_mat ES_mat PPCmedianperm_mat PPCmeanperm_mat PPCstdperm_mat PLVTimeEvts Phase_mat
    else
        warning('Analysis still running: Clusters  is not available yet for Subject %s %d - %d \n', SUBJECT, subj_i, nSUBJECTS)
    end
end

nClusters = numel(Clusters_all);
%nonan = ~isnan(nSign_pairs_all); % to manage the lack of data for patient 5


%% save results in DB

DB = struct();
DB.Clusters = Clusters_all;
DB.Pairs.Location = PairsLocation_MNI_all;
DB.Pairs.PPCmap = PPC_mat_all;
DB.Pairs.PPCzmap = PPCz_mat_all;
DB.Pairs.PPC_condmap = PPC_condmat_all;
DB.Pairs.PPCz_condmap = PPCz_condmat_all;
DB.Pairs.PPCz_condmap_values = PPCz_condmat_values_all;

DB.Pairs.medianPPCperm = PPCmedianperm_mat_all;
DB.Pairs.meanPPCperm = PPCmeanperm_mat_all;
DB.Pairs.stdPPCperm = PPCstdperm_mat_all;

DB.Pairs.ESmap = ES_mat_all;
DB.Pairs.Phasemap = Phase_mat_all;
DB.Pairs.TimeEvts = PLVTimeEvts_all;
DB.Pairs.nSignpairs = nSign_pairs_all;
DB.Pairs.nPairs = nPairs_all;
DB.Pairs.nSignpairs_perm = nSign_pairs_perm_all;


DB.Subjects = list_subj;



