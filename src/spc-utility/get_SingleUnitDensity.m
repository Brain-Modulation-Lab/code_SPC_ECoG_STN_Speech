function [SPCbandsunits_Density,BandSpecifity_btsp] = get_SingleUnitDensity(DB, nperms)


% get units pairs
[Units, idxUnitPairs] = getUniquePairs(DB.Pairs.Location);
nUnits = size(Units,1);
nUniquePairs = accumarray(idxUnitPairs,1);

% get unit clusters
cfg = [];
cfg.VarFields = {'Sign_Taskid','ClusterCentroidBand','ClusterSChannel','ClusterSubjects','ClusterOnOff','ClusterCentroid'};
cfg.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg.T = [-2.5 2];

[Sign_Task,ClusterCentroidBand,Cluster_SUnit, Cluster_Subject,ClusterOnOff,ClusterCentroid] = get_SPCClustersProperties(cfg,DB.Clusters);
%Units_FRmod = arrayfun(@(x) , unique(idxUnitPairs));


%[UnitsClusters, IdxUnitClusters] = getUniquePairs(DB.Clusters);
SPCbeta = nan(4,nUnits);
UnitFRMod = nan(1,nUnits);
SPCbandsunits_Density = table();

BandSpecifity_btsp = [];

nBands_units =nan(1,nUnits);

for unit_i = 1 : nUnits
    idxU = all(contains([Cluster_SUnit string(Cluster_Subject)],Units(unit_i,:)),2);
    sum(idxU)
    UnitClusterCentroidBand = ClusterCentroidBand(idxU);
    UnitClusterCentroidFrequency = ClusterCentroid(idxU,2);
    UnitClusterOnOff = ClusterOnOff(idxU,:)';
    UnitSign_Task = Sign_Task(:,idxU);
    nBands_units(unit_i) = numel(unique(UnitClusterCentroidBand));
    %tmp = ClusterS_FRMod(idxU == 1);
    UnitFRMod(unit_i) = mean(DB.Pairs.Location.S_typeFRmod(idxUnitPairs == unit_i));
    SPCbeta(4,unit_i) = sum(UnitClusterCentroidBand == 4 & UnitSign_Task(5,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i);
    SPCbeta(3,unit_i) = sum(UnitClusterCentroidBand == 4 & UnitSign_Task(1,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i);
    SPCbeta(2,unit_i) = sum(UnitClusterCentroidBand == 3 & UnitSign_Task(4,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i);
    SPCbeta(1,unit_i) = sum(UnitClusterCentroidBand == 2 & UnitSign_Task(4,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i);

    SPCbandsunits_Density.Unit(unit_i) = Units(unit_i,2);
    SPCbandsunits_Density.Subject(unit_i)  = Units(unit_i,1);
    SPCbandsunits_Density.nPairs(unit_i)  = nUniquePairs(unit_i);
    SPCbandsunits_Density.FRmod(unit_i) = UnitFRMod(unit_i);
    SPCbandsunits_Density.Density(unit_i)  = numel(UnitClusterCentroidBand)/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.theta(unit_i)  = sum(UnitClusterCentroidBand == 2)/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.alpha(unit_i)  = sum(UnitClusterCentroidBand == 3)/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.beta(unit_i)  = sum(UnitClusterCentroidBand == 4)/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.gammaL(unit_i)  = sum(UnitClusterCentroidBand == 5)/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.gammaH(unit_i)  = sum(UnitClusterCentroidBand == 6)/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.theta_ITI(unit_i) = sum(UnitClusterCentroidBand == 2 & UnitSign_Task(1,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.alpha_ITI(unit_i) = sum(UnitClusterCentroidBand == 3 & UnitSign_Task(1,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.theta_speech(unit_i)  = sum(UnitClusterCentroidBand == 2 & UnitSign_Task(4,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.alpha_speech(unit_i)  = sum(UnitClusterCentroidBand == 3 & UnitSign_Task(4,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;

    SPCbandsunits_Density.thetaalpha_ITI(unit_i)  = sum(UnitClusterCentroidBand < 4 & UnitSign_Task(1,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;

    SPCbandsunits_Density.beta_ITI(unit_i)  = sum(UnitClusterCentroidBand == 4 & UnitSign_Task(1,:) == 1 & UnitClusterOnOff(1,:) >= cfg.T(1) & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;
    SPCbandsunits_Density.beta_speech(unit_i)  = sum(UnitClusterCentroidBand == 4 & UnitClusterOnOff(1,:) >= -0.2 & UnitClusterOnOff(2,:) <= 0.5)/nUniquePairs(unit_i)*100;

    SPCbandsunits_Density.beta_postspeech(unit_i)  = sum(UnitClusterCentroidBand == 4  & UnitClusterOnOff(1,:) >= 0.51 & UnitClusterOnOff(2,:) <= cfg.T(end))/nUniquePairs(unit_i)*100;

    SPCbandsunits_Density.meanfreq(unit_i)  = mean(UnitClusterCentroidFrequency);
    SPCbandsunits_Density.medianfreq(unit_i)  = median(UnitClusterCentroidFrequency);

    freq_clus = [SPCbandsunits_Density.theta(unit_i) SPCbandsunits_Density.alpha(unit_i) SPCbandsunits_Density.beta(unit_i) ];





    prob_clus = freq_clus/sum(freq_clus);
    Entropyband = sum(-prob_clus.*log2(prob_clus + eps));
    if Entropyband < 0
        Entropyband = 0;
    end
    SPCbandsunits_Density.band_spec(unit_i)  =  1 - Entropyband/log2(3);

    for btsp_i = 1 : nperms

         UnitClusterCentroidBand_perm = randi(6,1,nUniquePairs(unit_i));UnitClusterCentroidBand_perm = UnitClusterCentroidBand_perm + 1;
     
        theta_btsp  = sum(UnitClusterCentroidBand_perm == 2)/nUniquePairs(unit_i)*100;
        alpha_btsp  = sum(UnitClusterCentroidBand_perm == 3)/nUniquePairs(unit_i)*100;
        beta_btsp  = sum(UnitClusterCentroidBand_perm == 4)/nUniquePairs(unit_i)*100;
        freq_clus_btsp = [theta_btsp alpha_btsp beta_btsp ];


        prob_clus_btsp = freq_clus_btsp/sum(freq_clus_btsp);
        Entropyband_btsp = sum(-prob_clus_btsp.*log2(prob_clus_btsp + eps));
        BandSpecifity_btsp = [BandSpecifity_btsp 1 - Entropyband_btsp/log2(3)];

    end


end