function varargout = get_SPCClustersProperties(cfg, Clusters_all)

VarFields = cfg.VarFields;
nVarField = numel(VarFields);
nClusters = numel(Clusters_all);


ClusternSign = [Clusters_all.nSign];
nClusters_all = sum(ClusternSign);

bands = bml_get_canonical_bands();

for var_i = 1 : nVarField

    VarField = VarFields{var_i};

    switch VarField
        case'PairId'
             varargout{var_i} = [Clusters_all.pair_i];           
        case 'PPC'
            varargout{var_i} = [Clusters_all.PPC];
        case 'PPCz'
            varargout{var_i} = [Clusters_all.PPCz];
        case 'ES'
            varargout{var_i} = [Clusters_all.ES];
        case 'ClusterDur'
            varargout{var_i} = [Clusters_all.TimeSpan];
        case 'ClusterCycles'
            varargout{var_i} = [Clusters_all.nCycles];
        case 'Zstat'
            varargout{var_i} = [Clusters_all.Zstat];
        case 'ClusterSignEffect'
            varargout{var_i} = sign([Clusters_all.Zstat_rel]);
        case 'FreqSpread'
            varargout{var_i} = [Clusters_all.FreqSpread];
        case 'nSign'
            varargout{var_i} = [Clusters_all.nSign];
        case 'PhaseFlip'
            varargout{var_i} = [Clusters_all.flip_polarity];
        case 'PhaseinCluster'
            ClusterPhaseinCluster = {Clusters_all.PhaseInCluster};
            varargout{var_i} = [ClusterPhaseinCluster{:}];
        case 'ClusterPhase'
            ClusterPhase = [];
            for clus_i = 1 : nClusters
                ClusterPhase =[ ClusterPhase; Clusters_all(clus_i).Phase];
            end
            varargout{var_i} = ClusterPhase;
        case 'ClusterCentroid'
            ClusterCentroid = [];
            for clus_i = 1 : nClusters
                ClusterCentroid =[ClusterCentroid; Clusters_all(clus_i).Centroid];
            end
            varargout{var_i} = ClusterCentroid;

        case  'ClusterOnOff'
            ClusterOnOff = [];
            for clus_i = 1 : nClusters
                ClusterOnOff = [ClusterOnOff; [Clusters_all(clus_i).On' Clusters_all(clus_i).Off']];
            end
            varargout{var_i} = ClusterOnOff;

        case 'centerEvts'
            centerEvts = [];
            for clus_i = 1 : nClusters
                centerEvts = [centerEvts; Clusters_all(clus_i).centerEvts'];
            end
            varargout{var_i} = centerEvts;

        case 'ClusterE_MNI'
            ClusterE_MNI = [];
            for clus_i = 1 : nClusters
                ClusterE_MNI = [ ClusterE_MNI; repmat([Clusters_all(clus_i).E_MNI_X,Clusters_all(clus_i).E_MNI_Y, Clusters_all(clus_i).E_MNI_Z],ClusternSign(clus_i),1)];
            end
            varargout{var_i} = ClusterE_MNI;

        case 'ClusterS_MNI'
            ClusterS_MNI = [];
            for clus_i = 1 : nClusters
                ClusterS_MNI = [ ClusterS_MNI; repmat([Clusters_all(clus_i).S_MNI_X,Clusters_all(clus_i).S_MNI_Y, Clusters_all(clus_i).S_MNI_Z],ClusternSign(clus_i),1)];
            end
            varargout{var_i} = ClusterS_MNI;
        case 'ClusterPair'
            ClusterPair = [];
            for clus_i = 1 : nClusters
                ClusterPair = [ClusterPair; repmat(Clusters_all(clus_i).pair_i,ClusternSign(clus_i),1)];
            end
            varargout{var_i} = ClusterPair;

        case 'ClusterSubjects'
            ClusterSubjects = [];
            for clus_i = 1 : nClusters
                ClusterSubjects = [ClusterSubjects; repmat(Clusters_all(clus_i).subj_id,ClusternSign(clus_i),1)];
            end
            ClusterSubjects = cellstr(ClusterSubjects);
            varargout{var_i} = ClusterSubjects;

        case     'ClusterEAtlasLabel'
            ClusterEAtlasLabel = {};
            for clus_i = 1 : nClusters
                ClusterEAtlasLabel = [ClusterEAtlasLabel; repmat(Clusters_all(clus_i).E_atlas_label_Destrieux,ClusternSign(clus_i),1)];
            end
            varargout{var_i} = ClusterEAtlasLabel;

        case 'ClusterEvtTimes'
            ClusterEvtTimes = [];
            for clus_i = 1 : nClusters
                ClusterEvtTimes = [ClusterEvtTimes; repmat(median(Clusters_all(clus_i).TimeEvts.Evts),ClusternSign(clus_i),1)];
            end
            varargout{var_i} = ClusterEvtTimes;

        case 'ClusterEAtlasLabel_HCMMP'
            ClusterEAtlasLabel_HCMMP = {};
            for clus_i = 1 : nClusters
                ClusterEAtlasLabel_HCMMP = [ClusterEAtlasLabel_HCMMP; repmat(Clusters_all(clus_i).E_HCPMMP1_label_1,ClusternSign(clus_i),1)];
            end
            varargout{var_i} = ClusterEAtlasLabel_HCMMP;

        case 'ClusterS_FRMod'
            ClusterS_FRMod = [];
            for clus_i = 1 : nClusters
                ClusterS_FRMod = [ClusterS_FRMod;  repmat(Clusters_all(clus_i).S_typeFRmod,ClusternSign(clus_i),1)];
            end
            varargout{var_i} = ClusterS_FRMod;

        case 'ClusterSChannel'
            ClusterSChannel = [];
            for clus_i = 1 : nClusters
                ClusterSChannel = [ClusterSChannel; repmat(string(Clusters_all(clus_i).S_channel),ClusternSign(clus_i),1)];
            end
            varargout{var_i} = ClusterSChannel;

        case 'ClusterEChannel'
            ClusterEChannel = [];
            for clus_i = 1 : nClusters
                ClusterEChannel = [ClusterEChannel; repmat(string(Clusters_all(clus_i).E_channel),ClusternSign(clus_i),1)];
            end
            varargout{var_i} = ClusterEChannel;

        case 'ClusterSunitGrade'
            ClusterSunitGrade = {};
            for clus_i = 1 : nClusters
                ClusterSunitGrade = [ClusterSunitGrade; repmat(string(Clusters_all(clus_i).S_unitGrade),ClusternSign(clus_i),1)];
            end
            varargout{var_i} = ClusterSunitGrade;
        case 'ClusterCentroidBand'
            ClusterCentroid = [];
            for clus_i = 1 : nClusters
                ClusterCentroid =[ClusterCentroid; Clusters_all(clus_i).Centroid];
            end


            ClusterCentroidBand = nan(1,nClusters_all);
            for clus_i = 1 : nClusters_all
                ClusterCentroidBand(clus_i) = find(ClusterCentroid(clus_i,2) < bands.fends,1);
            end
            varargout{var_i} = ClusterCentroidBand;

            
        case 'Sign_Taskid'
            Task_labels = cfg.Task_labels;

            Sign_Task = [];
            for clus_i = 1 : nClusters
                Sign_Task =  [Sign_Task ;Clusters_all(clus_i).window_task(:)];
            end

            Sign_Taskid = zeros(5,nClusters_all);
            for ii = 1 : nClusters_all
                tmp = Sign_Task{ii};
                for kk = 1:numel(tmp)
                    Sign_Taskid(find(strcmpi(Task_labels,tmp{kk})),ii) = 1;
                end
            end
            varargout{var_i} = Sign_Taskid;
      
        case 'clust_prob'
            
            T = cfg.T;
%             if ~contains(VarFields,'ClusterOnOff')
                ClusterOnOff = [];
                for clus_i = 1 : nClusters
                    ClusterOnOff = [ClusterOnOff; [Clusters_all(clus_i).On' Clusters_all(clus_i).Off']];
                end
%             end
            cluster_prob = zeros(nClusters_all,numel(T));
            for ii = 1 : nClusters_all
                %     cluster_prob(ii,T<= (ClusterOnOff(ii,2) - meanSpeechOnset) & T >= (ClusterOnOff(ii,1) - meanSpeechOnset)) = 1;
                cluster_prob(ii,T<= (ClusterOnOff(ii,2)) & T >= (ClusterOnOff(ii,1))) = 1;
            end
            varargout{var_i} = cluster_prob;

    end
end
