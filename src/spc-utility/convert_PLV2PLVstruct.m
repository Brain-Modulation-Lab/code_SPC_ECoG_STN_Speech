function PLV_struct = convert_PLV2PLVstruct(PLV,cfg_analysis)


% get fields from cfg_analysis
MNI = cfg_analysis.MNI;
sortnotes = cfg_analysis.sortnotes;
FR_STATS = cfg_analysis.FR_STATS;
SUBJECT = cfg_analysis.SUBJECT;
cfg = cfg_analysis.cfg;

if isstruct(FR_STATS) % just a check, ensure that is a table
    FR_STATS = struct2table(FR_STATS);
end


n_pairs = numel(PLV);
% initialize counter
cont = 1;

fprintf(" Stacking # pairs = %d in a struct: \n",n_pairs);
PLV_struct = struct();

for pair_i = 1 : n_pairs

    if mod(pair_i, 10) == 0
        fprintf(" Stacking %d - %d pairs [%1.2f %%]... \n",pair_i, n_pairs, round(pair_i/n_pairs*100,2));
    end
    PLV_struct(cont).id = cont;
    PLV_struct(cont).SUBJECT = SUBJECT;
    PLV_struct(cont).E_channel = PLV(pair_i).E_channel;
    PLV_struct(cont).S_channel = PLV(pair_i).S_channel;
    PLV_struct(cont).n_trials = PLV(pair_i).n_trials;

    PLV_struct(cont).winCenters = PLV(pair_i).winCenters;
    PLV_struct(cont).centerEvts = PLV(pair_i).centerEvts;
    PLV_struct(cont).freqsOfInt = PLV(pair_i).freqsOfInt;
    PLV_struct(cont).flip_polarity = PLV(pair_i).flip_polarity;
    PLV_struct(cont).SPC = PLV(pair_i).store;
    PLV_struct(cont).cfg  = cfg; 


    % get unit from sortnotes
    idx_unit = regexp(PLV_struct(cont).S_channel,'\d*','Match');
    idx_unit = sortnotes.id == str2double(idx_unit{1}); % modified for the demo
    % idx_unit = str2double(idx_unit{1}); % previous version
    

    % save S info
    PLV_struct(cont).S_unitType = sortnotes.unitType{idx_unit};
    PLV_struct(cont).S_unitGrade = sortnotes.unitGrade{idx_unit};
    PLV_struct(cont).S_starts = sortnotes.starts(idx_unit);
    PLV_struct(cont).S_ends = sortnotes.ends(idx_unit);
    PLV_struct(cont).S_duration = sortnotes.duration(idx_unit);
    PLV_struct(cont).S_session = sortnotes.session_id(idx_unit);
    PLV_struct(cont).S_traj = string(sortnotes.channel(idx_unit));

    PLV_struct(cont).S_depth = MNI.S{contains(MNI.S.electrode,PLV_struct(cont).S_traj) & MNI.S.session_id == PLV_struct(cont).S_session,{'depth'}};
    PLV_struct(cont).S_MNI_X = MNI.S{contains(MNI.S.electrode,PLV_struct(cont).S_traj) & MNI.S.session_id == PLV_struct(cont).S_session,{'mni_nonlinear_x'}};
    PLV_struct(cont).S_MNI_Y = MNI.S{contains(MNI.S.electrode,PLV_struct(cont).S_traj) & MNI.S.session_id == PLV_struct(cont).S_session,{'mni_nonlinear_y'}};
    PLV_struct(cont).S_MNI_Z = MNI.S{contains(MNI.S.electrode,PLV_struct(cont).S_traj) & MNI.S.session_id == PLV_struct(cont).S_session,{'mni_nonlinear_z'}};
    PLV_struct(cont).S_DISTAL_LABEL1 = MNI.S{contains(MNI.S.electrode,PLV_struct(cont).S_traj) & MNI.S.session_id == PLV_struct(cont).S_session,{'DISTAL_label_1'}};
    PLV_struct(cont).S_DISTAL_LABEL2 = MNI.S{contains(MNI.S.electrode,PLV_struct(cont).S_traj) & MNI.S.session_id == PLV_struct(cont).S_session,{'DISTAL_label_2'}};

    idx_stats = FR_STATS.unit_id == idx_unit & FR_STATS.session == PLV_struct(cont).S_session & contains(FR_STATS.SubjectID,SUBJECT);
    if sum(idx_stats)> 0
        if cellfun(@(x) ~isempty(x),FR_STATS.Excit(idx_stats)) && cellfun(@(x) ~isempty(x),FR_STATS.Inhib(idx_stats))
            PLV_struct(cont).S_typeFRmod = 2;
        elseif cellfun(@(x) ~isempty(x),FR_STATS.Excit(idx_stats)) && cellfun(@(x) isempty(x),FR_STATS.Inhib(idx_stats))
            PLV_struct(cont).S_typeFRmod = 1;
        elseif cellfun(@(x) isempty(x),FR_STATS.Excit(idx_stats)) && cellfun(@(x) ~isempty(x),FR_STATS.Inhib(idx_stats))
            PLV_struct(cont).S_typeFRmod = -1;
        elseif cellfun(@(x) isempty(x),FR_STATS.Excit(idx_stats)) && cellfun(@(x) isempty(x),FR_STATS.Inhib(idx_stats))
            PLV_struct(cont).S_typeFRmod = 0;
        end
    else
        PLV_struct(cont).S_typeFRmod = nan;
    end

    % save E info

    % PLV_struct(cont).ecogChannel = PLV.lfplabel{chan_i,1};%PPC.speech{1,unit_i}.labelcmb{chan_i,2};

    idx_electrode = contains(MNI.E.electrode,PLV_struct(cont).E_channel);
    PLV_struct(cont).E_strip = MNI.E.strip(idx_electrode);
    PLV_struct(cont).E_target = MNI.E.target(idx_electrode);
    PLV_struct(cont).E_side = MNI.E.side(idx_electrode);
    PLV_struct(cont).E_MNI_X = MNI.E.mni_nonlinear_x(idx_electrode);
    PLV_struct(cont).E_MNI_Y = MNI.E.mni_nonlinear_y(idx_electrode);
    PLV_struct(cont).E_MNI_Z = MNI.E.mni_nonlinear_z(idx_electrode);
    PLV_struct(cont).E_atlas_label_Desikan = MNI.E.atlas_label_Desikan(idx_electrode);
    PLV_struct(cont).E_atlas_label_Destrieux = MNI.E.atlas_label_Destrieux(idx_electrode);

    PLV_struct(cont).E_HCPMMP1_label_1 = MNI.E.HCPMMP1_label_1(idx_electrode);
    PLV_struct(cont).E_HCPMMP1_weight_1 = MNI.E.HCPMMP1_weight_1(idx_electrode);
    PLV_struct(cont).E_HCPMMP1_label_2 = MNI.E.HCPMMP1_label_2(idx_electrode);
    PLV_struct(cont).E_HCPMMP1_weight1 = MNI.E.HCPMMP1_weight_2(idx_electrode);


    % update cont
    cont = cont+1;
end
clear PLV data

fprintf('PLV struct has %d rows \n',numel(PLV_struct))


