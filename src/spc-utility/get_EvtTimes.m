function EvtTimes = get_EvtTimes(Pairs)
EvtTimes = struct();
EvtTimes.Cue = [];
EvtTimes.Speech = [];
% ClusternSign = [Clusters.nSign];
% nClusters = numel(Clusters);
% for clus_i = 1 : nClusters
%     EvtTimes = [EvtTimes; repmat(median(Clusters(clus_i).TimeEvts.Speech.Evts),ClusternSign(clus_i),1)];
% end

for pair_i = 1 : numel(Pairs.TimeEvts)
    EvtTimes.Cue = [EvtTimes.Cue; median(Pairs.TimeEvts(pair_i).Cue.Evts)];
    EvtTimes.Speech = [EvtTimes.Speech; median(Pairs.TimeEvts(pair_i).Speech.Evts)];  
end
EvtTimes.Cue = median(EvtTimes.Cue);
EvtTimes.Speech = median(EvtTimes.Speech);
