% This script calculates the time-resolved spike-phase coupling for 2
% exemplary ECoG channels and 2 exemplary STN neurons.

% Clear command window, workspace, and close all figures
clc
clear all
close all

% Define paths for script, data, and output
script_folder = pwd;  % Current folder of the script
[parentFolder,~,~] = fileparts(script_folder);
PATH_DATA = fullfile(script_folder,'intracranial-data-examples');  % Path to data folder
PATH_OUTPUT = fullfile(script_folder,'output');  % Path to output folder

% Create output folder if it doesn't exist
if ~isfolder(PATH_OUTPUT)
    mkdir(PATH_OUTPUT);
end

% Define paths for source scripts and data files
PATH_SPC = fullfile(parentFolder,'src');  % Path to source folder
PATH_DATASOURCE = fullfile(PATH_DATA ,'intracranial_data.mat'); 
% Path to ECoG data
% Path to spike data
% Path to spike sorting notes
% Path to anatomical localization
% Path to firing rate stats


% Add source path and its subdirectories to the search path
addpath(genpath(PATH_SPC))

% Ensure required toolboxes are added to the path
prettify_plots;  % Custom function to improve plot appearance
ft_defaults;     % Initialize FieldTrip defaults
bml_defaults;    % Initialize BML toolbox defaults

% Get configuration settings (edit set_configs.m for customization)
cfg = set_configs('default');

% Display script start time
fprintf(" Script starts at %s \n", datetime("now"))
tStart = tic;  % Start timer

% Load ECoG data (E), neuron spikes (S), spike sorting notes (sortnotes), 
% anatomical localization (MNI), firing rate stats (FR_STATS)
load(PATH_DATASOURCE)

% select participant
SUBJECT = 'DBS3008'; % exemplary data are from this patient

% Compute Phase Locking Value (PLV)
PLV = calc_spike_PLV_all(E, S, cfg);  % Main function for PLV calculation (it might take a while)

% Extract t-SPC information and pair information from the struct (data will
% be saved in a struct called ClusterPLV.mat)
cfg_s = [];
cfg_s.MNI = MNI;
cfg_s.sortnotes = sortnotes;
cfg_s.FR_STATS = FR_STATS;
cfg_s.SUBJECT = SUBJECT;
cfg_s.PATH_OUTPUT = PATH_OUTPUT;
cfg_s.cfg = cfg;
storeSPCinformation(PLV, cfg_s)

% Plot PLV results
for ii = 1 : numel(PLV)
    plot_spike_PLV(PLV(ii)) % plot each pair
end


% End of script
tEnd = toc(tStart);  % End timer
disp(" Script ends ...")
fprintf(" Time elapsed %.2f \n",tEnd)