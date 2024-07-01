function storeSPCinformation(PLV, cfg_analysis)

PATH_OUTPUT = cfg_analysis.PATH_OUTPUT;

% Convert to PLV into a convenient struct
fprintf('Convert to PLV into a convenient struct')
PLV_struct = convert_PLV2PLVstruct(PLV,cfg_analysis);


% define variables of interest
toJoin = {'S_channel','E_channel','centerEvts','flip_polarity','S_unitType','S_unitGrade','S_depth','S_MNI_X',...
    'S_MNI_Y','S_MNI_Z','S_DISTAL_LABEL1','S_typeFRmod','E_MNI_X', ...
    'E_MNI_Y','E_MNI_Z','E_atlas_label_Destrieux','E_atlas_label_Desikan','E_HCPMMP1_label_1'};
toJoinMinimal = {'S_depth','S_MNI_X',...
    'S_MNI_Y','S_MNI_Z','S_typeFRmod','S_DISTAL_LABEL1','E_MNI_X', ...
    'E_MNI_Y','E_MNI_Z','E_atlas_label_Destrieux','E_atlas_label_Desikan','E_HCPMMP1_label_1'};


fprintf('Starting Cluster analysis')

Clusters = [];
PairsLocation_MNI = table();
PPC_mat = {};
PPCz_mat = {};
PPCmedianperm_mat = {};
PPCmeanperm_mat = {};
PPCstdperm_mat = {};

ES_mat = {};
PLVTimeEvts = struct();
Phase_mat  ={};

% find clusters
n_pairs = numel(PLV_struct);
fprintf(' Save information about location ')

fprintf(' Computing clusters acroos %d pairs... \n',n_pairs )
Clusters = arrayfun(@(x) comp_clusterPLV(x.SPC.PPC,x.SPC.PPC_perm, 'zstat'),PLV_struct);
nSign = [Clusters.nSign];

% save information about pairs location
PairsLocation_MNI.id = [1:n_pairs]';
PairsLocation_MNI.nSign = nSign';

for ff = 1 : numel(toJoinMinimal)
    PairsLocation_MNI.(toJoinMinimal{ff}) = [PLV_struct.(toJoinMinimal{ff})]';
end

sum(isnan(PairsLocation_MNI.S_typeFRmod))
for kk = 1 : n_pairs
    PPC_mat{kk} =  PLV_struct(kk).SPC.PPC;
    PPCz_mat{kk} = (PLV_struct(kk).SPC.PPC - mean(PLV_struct(kk).SPC.PPC_perm,3,'omitnan')) ./ std(PLV_struct(kk).SPC.PPC_perm,[],3,'omitnan');
    PPCmedianperm_mat{kk} = squeeze(median(PLV_struct(kk).SPC.PPC_perm,3,'omitnan'));
    PPCmeanperm_mat{kk} = mean(PLV_struct(kk).SPC.PPC_perm,3,'omitnan');
    PPCstdperm_mat{kk} = std(PLV_struct(kk).SPC.PPC_perm,[],3,'omitnan');
    ES_mat{kk} =  PLV_struct(kk).SPC.ES;
    Phase_mat{kk} = PLV_struct(kk).SPC.phase;

    % get timing
    PLVTimeEvts(kk).Cue.time = mean(PLV_struct(kk).winCenters -PLV_struct(kk).winCenters(:,PLV_struct(kk).centerEvts(2))) ;
    PLVTimeEvts(kk).Speech.time = mean(PLV_struct(kk).winCenters -PLV_struct(kk).winCenters(:,PLV_struct(kk).centerEvts(4))) ;
    for ev = 1 : 6
        PLVTimeEvts(kk).Cue.Evts(:,ev)= PLV_struct(kk).winCenters(:,PLV_struct(kk).centerEvts(ev)) -PLV_struct(kk).winCenters(:,PLV_struct(kk).centerEvts(2));
        PLVTimeEvts(kk).Speech.Evts(:,ev)= PLV_struct(kk).winCenters(:,PLV_struct(kk).centerEvts(ev)) -PLV_struct(kk).winCenters(:,PLV_struct(kk).centerEvts(4));
    end

    % get pair information
    PairsLocation_MNI.S_channel(kk) = string(PLV_struct(kk).S_channel);
    PairsLocation_MNI.E_channel(kk) = string(PLV_struct(kk).E_channel);

end
fprintf(' Adjusting clusters... \n')
% adjust clusters
idx_sign = find([Clusters.nSign] > 0);
Clusters = Clusters(idx_sign);
tmp = num2cell(idx_sign);
[Clusters.pair_i] = deal(tmp{:});
PLV_struct_toDo = PLV_struct(idx_sign);
PLVTimeEvts_toDo = PLVTimeEvts(idx_sign);
nSign_pairs = numel(PLV_struct_toDo);
fprintf(' Identified %d-%d clusters pairs... \n',nSign_pairs, n_pairs )
% define clusters properties: task location, time-frequency
fprintf(' Computing clusters properties... \n')
Clusters = arrayfun(@(x,y) comp_propTFclusters(x,y),Clusters, PLV_struct_toDo);
% merge important information from PLV

for clus_i = 1 : nSign_pairs
    for ff = 1 : numel(toJoin)
        Clusters(clus_i).(toJoin{ff}) = PLV_struct_toDo(clus_i).(toJoin{ff});
    end
    Clusters(clus_i).TimeEvts = PLVTimeEvts_toDo(clus_i);
end

fprintf(' ...Completed clusters properties \n')

PATH_clusters = fullfile(PATH_OUTPUT,'ClustersPLV.mat');
fprintf(' Saving results in %s \n', PATH_clusters)
save(PATH_clusters,'Clusters', 'PLV_struct_toDo', 'nSign_pairs','n_pairs','PairsLocation_MNI', 'PPC_mat','PPCmedianperm_mat','PPCmeanperm_mat','PPCstdperm_mat','PPCz_mat','ES_mat','PLVTimeEvts','Phase_mat');
fprintf(' Done  \n')