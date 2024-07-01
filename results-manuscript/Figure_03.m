%% This code reproduces the Figure 3 of the manuscript. This figure describes
% the spatial distribution of the spike-phase coupling on the cortex

clc
clear all
close all

% Define paths for script, data, and output
script_folder = pwd;  % Current folder of the script
[parentFolder,~,~] = fileparts(script_folder);


% Define paths for source scripts and data files
PATH_SPC = fullfile(parentFolder,'src');  % Path to source folder
PATH_DATA = fullfile(parentFolder,'data');  % Path to source folder
PATH_DB = fullfile(PATH_DATA,'DB_main_analysis.mat');  % Path to source folder

% Add source path and its subdirectories to the search path
addpath(genpath(PATH_SPC))

% initialize libraries
bml_defaults;
ft_defaults;

% make plot prettier
prettify_plots;

% define frequency bands
bands = bml_get_canonical_bands([0,150]);
bands.fmid = sqrt(bands.fstarts .* bands.fends);
newcolor = {'#d7191c', '#fdae61','#fecc5c','#abdda4','#2b83ba'};
bands.color(2:6) = newcolor;



% Define your participant list
SUBJECTS = {'DBS3001','DBS3002','DBS3003','DBS3004','DBS3005','DBS3008', 'DBS3010', 'DBS3011', ...
            'DBS3012', 'DBS3014', 'DBS3015', 'DBS3017', 'DBS3018', 'DBS3019', 'DBS3020', 'DBS3022', ...
            'DBS3023', 'DBS3024', 'DBS3026', 'DBS3027', 'DBS3028', 'DBS3029', 'DBS3030', 'DBS3031', ...
            'DBS3032'};


% Step 1. gather all SPC (pairs & clusters) information
load(PATH_DB,'DB')
cfg = set_configs('default');

% Step 2. Get times of events
EvtTimes = get_EvtTimes(DB.Pairs);


% time resolution
clust_res = .005;
T = [-2.5 2];

% this may take a while (~20 min on personal laptop)

% Panel Figure 3B (panelsa arounds cortex: Plot t-SPC occurrence in cortical
% ROIs:
%     {'G_precentral L'         }
%     {'G_front_middle L'       }
%     {'G_pariet_inf-Supramar L'}
%     {'G_temp_sup-Lateral L'   }
%     {'G_postcentral L'        }
%     {'G_and_S_subcentral L'   }
%     {'G_front_inf-Opercular L'}
%     {'G_front_inf-Triangul L' }
%     {'G_temporal_middle L'    }
%     {'G_temp_sup-Plan_tempo L'}

cfg_sp = struct();
cfg_sp.time = 'all'; 
cfg_sp.cfg = cfg;
cfg_sp.atlas = 'E_atlas_label_Destrieux';
cfg_sp.min_subj = 7;
cfg_sp.min_pairs = 100;
cfg_sp.flag_sort = true;
cfg_sp.nperms = 500;
cfg_sp.bands = bands;
cfg_sp.EvtTimes = EvtTimes.Speech;
cfg_sp.T = T(1) : clust_res : T(2);
cfg_sp.Visible = 'off'; % to visualize them put 'on'
[fh3B2, ROI_AtlasLabels] = create_SPCROI_time(cfg_sp, DB);

% Panel Figure 3A-B and Supplementary Figures 8-9 (8 panel total)
cfg_sp = struct();
cfg_sp.bands = bands;
cfg_sp.map = fullfile(PATH_DATA,'tSPC_density_ECoG.txt');
cfg_sp.ref = fullfile(PATH_DATA,'cortex_MNI.mat');
cfg_sp.nperms = 500;
cfg_sp.atlas = 'E_atlas_label_Destrieux';
cfg_sp.min_subj = 7;
cfg_sp.min_pairs = 100;
cfg_sp.flag_sort = true;
fh3AB = create_SPCCortex_CSFS(cfg_sp, DB); %

% alternatively you can project the .node file on SurfIce.





