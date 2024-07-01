%% This code reproduces the Figure 3 of the manuscript. This figure describes
% the error analysis in the manuscript

clc
clear all
close all

% Define paths for script, data, and output
script_folder = pwd;  % Current folder of the script
[parentFolder,~,~] = fileparts(script_folder);


% Define paths for source scripts and data files
PATH_SPC = fullfile(parentFolder,'src');  % Path to source folder
PATH_DATA = fullfile(parentFolder,'data');  % Path to source folder
PATH_DB = fullfile(PATH_DATA,'DB_error_analysis.mat');  % Path to source folder

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
load(PATH_DB,'DB_acc')
cfg = set_configs("default-5_30Hz");

% Step 2. Get times of events
EvtTimes = get_EvtTimes(DB_acc.Pairs);

which_var = 'PPCz';
which_pairs = 'all';
Evt_lock = 'Speech';
T = [-2.5 2];
Tres = 0.050;
nPerms = 1000;
clust_res = .005;
Tc = T(1) : clust_res : T(2);


% Figure 5 panel A: t-score map of the difference bewteen error and accurate
% trials
cfg.stat_test = 'paired-test-conds';
cfg.Filestats = fullfile(PATH_DATA,'permutation_avgmaps',[which_var,'Stat_behavsplit-phonetic_ontarget_pairs-',which_pairs,'_evt-',Evt_lock,'_stat-',cfg.stat_test,'_nperm-',num2str(nPerms),'_new-implementation.mat']); % use the pre=saved permutaiton because they take long time!

% pre-computed stats
SPCmap = get_SPCmap(DB_acc.Pairs,'PPCz_cond', T, ...
    Tres,cfg,'Speech','All',nPerms);

% flip to change sign (flip colors)
SPCmap.map = -1*SPCmap.map;
SPCmap.Stat.TscoreMap = -1*SPCmap.Stat.TscoreMap;

cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTimes = EvtTimes.Speech;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CLim = {[0 .3],[-5 5]};
cfg_fh.CMap = redwhitegreen;
cfg_fh.Visible = 'off';
cfg_fh.NClustplot = 0;

fh5A = plotter_SPCmap(SPCmap, cfg_fh);
set(fh{2},'visible','on')


%% Figure 5 panel B: Cluster occurrence in theta-alpha and beta in error vs accurate trials ( this might take a while)
%  % fh{1}: theta-alpha
%  % fh{2}: beta

cfg_fh = [];
cfg_fh.VarFields = {'ClusterCentroidBand','clust_prob'};
cfg_fh.nperms = nPerms; 
cfg_fh.T = Tc;
cfg_fh.bands = {[2 3],4};
fh5B = compare_ClusterOccurrence(DB_acc,cfg_fh);

%% Figure 5 panel C: comapre cumulative probability of cluster onset and offset in error vs accurate trials
cfg_fh = [];
cfg_fh.VarFields = {'ClusterCentroidBand','ClusterOnOff'};
cfg_fh.T = Tc;
cfg_fh.bands = [2 3];
fh5C = compare_ClusterEcdf(DB_acc,cfg_fh);

%% Figure 5 panelS D-E: comapre cluster timing and duration in error vs accurate trials
%  % fh{1}: cluster onset
%  % fh{2}: cluster duration

cfg_fh = [];
cfg_fh.VarFields = {'ClusterCentroidBand','ClusterOnOff'};
cfg_fh.bands = [2 3];
fh = compare_ClusterTiming(DB_acc,cfg_fh);


