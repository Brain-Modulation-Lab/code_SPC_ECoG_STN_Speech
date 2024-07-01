function SPCpair = get_SPCpair(DB, cfg)

% to plot need:
% 1. PPCz
% 2. LabelMatrix
% 3. Event_names
% 4. Center events
% 5. Unit
% 6. Subject
% 7. Ecog
% 8. Ecog Label


Subjects = cfg.Subjects;
ClusterSubjects = {DB.Clusters.subj_id};
PairSubjects = cellstr(DB.Pairs.Location.subj_id)';

SPCpair = struct();
pair_id = [];
LabelMatrix = {};
CenterEvts = {};
Subj = {};
Unit = {};
ECoG = {};
ECoG_atlas = {};
PPCz = {};
PPC = {};
nSign = [];
for subj_i = 1 : numel(Subjects)
    % get participants
    subject = Subjects{subj_i};
    % get t-SPC and Pairs in a participant
    ClusterInSubject = DB.Clusters(strcmpi(ClusterSubjects,subject));
    LocationInSubject = DB.Pairs.Location(strcmpi(PairSubjects,subject),:);
    PPCInSubject = DB.Pairs.PPCmap(strcmpi(PairSubjects,subject));
    PPCzInSubject = DB.Pairs.PPCzmap(strcmpi(PairSubjects,subject));

    % get idx pairs
    idx_pairs = [ClusterInSubject.pair_i];
    pair_id = [pair_id idx_pairs];
    % get t-SPC properties
    LabelMatrix = cat(2,LabelMatrix, {ClusterInSubject.LabelMatrix});
    CenterEvts = cat(2,CenterEvts, {ClusterInSubject.centerEvts});
    Subj = cat(2,Subj, {ClusterInSubject.subj_id});
    Unit = cat(2,Unit, {ClusterInSubject.S_channel});
    ECoG = cat(2,ECoG, {ClusterInSubject.E_channel});
    ECoG_atlas = cat(2,ECoG_atlas, {ClusterInSubject.E_atlas_label_Destrieux});
    % use idx_pairs to connect pairs and t-SPC
    PPCz = cat(2,PPCz, PPCzInSubject(idx_pairs));
    PPC = cat(2,PPC, PPCInSubject(idx_pairs));
    nSign = [nSign LocationInSubject.nSign(idx_pairs)'];

end

% save it
SPCpair.pair_id = pair_id;
SPCpair.PPCz = PPCz;
SPCpair.PPC = PPC;
SPCpair.LabelMatrix = LabelMatrix;
SPCpair.nSign = nSign;
SPCpair.CenterEvts = CenterEvts;
SPCpair.freqQ =  logspace(log10(cfg.freqsOfInt(1,1)), log10(cfg.freqsOfInt(end,1)), size(SPCpair.PPCz{1,1},1));

SPCpair.Subj = Subj;
SPCpair.Unit = Unit;
SPCpair.ECoG = ECoG;
SPCpair.ECoG_atlas = ECoG_atlas;







