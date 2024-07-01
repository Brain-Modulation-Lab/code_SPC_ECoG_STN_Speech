function varargout = getUniquePairs(x)
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