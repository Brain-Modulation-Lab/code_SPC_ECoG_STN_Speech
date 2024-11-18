function varargout = getUniquePairsByModality(x, modality)
switch modality
    case "STN"

        if istable(x) % pairs type
            Subjects = x.subj_id;
            Units = x.S_channel;
        elseif isstruct(x) % clusters type
            cfg.VarFields = {'ClusterSubjects','ClusterSChannel'};
            [Subjects, Units] = get_SPCClustersProperties(cfg, x);

        end
        if nargout == 2
            [UniquePairs, ~,idxUniquePairs] = unique([Subjects, Units],"rows");
            varargout{1} = UniquePairs;
            varargout{2} = idxUniquePairs;
        elseif nargout == 1
            UniquePairs  = unique([Subjects, Units],"rows");
            varargout{1} = UniquePairs;
        end
    case "ECoG"
        if istable(x) % pairs type
            Subjects = x.subj_id;
            Electrodes = x.E_channel;
        elseif isstruct(x) % clusters type
            cfg.VarFields = {'ClusterSubjects','ClusterEChannel'};
            [Subjects, Electrodes] = get_SPCClustersProperties(cfg, x);

        end
        if nargout == 2
            [UniquePairs, ~,idxUniquePairs] = unique([Subjects, Electrodes],"rows");
            varargout{1} = UniquePairs;
            varargout{2} = idxUniquePairs;
        elseif nargout == 1
            UniquePairs  = unique([Subjects, Electrodes],"rows");
            varargout{1} = UniquePairs;
        end
end