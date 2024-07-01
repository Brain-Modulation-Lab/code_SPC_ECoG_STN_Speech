function [] = avg_numSpikes_PLV(Paths, Const)


allFiles = {'STN_20131003_H4', 'STN_20131010_H6', 'STN_20131011_H1', 'STN_20131011_H2', 'STN_20131011_H3', 'STN_20131010_H4', ...
            'STN_20131009_H3'};% , 'STN_20130130_C4'};
        
WIN_GET_SPIKENR = 0.2; % 0.2 s, can be smaller than in the patient data as we have more trials
Const.NBINS     = 6;   % used for the binning procedure, num of bins for each interval between two events


LOW_FREQ = false;
addit_4ms = '_4ms';
addit_4ms = '';



for f = 1:numel(allFiles)

    currFile = allFiles{f};
    currFileName = [Paths.sourceData, currFile, '\', currFile];
    if ~LOW_FREQ
        load([currFileName, '_downsamp', addit_4ms, '.mat']);
    elseif LOW_FREQ
        load([currFileName, '_downsamp_lowFreq', addit_4ms, '.mat']);
    end
    
    clear store
    SR = data.SR;

    %% Load spikes
    Evt_data = load([Paths.sourceData, currFile, '\', currFile, '.mat']);
    [Evt] = readEvts(Evt_data); 
    currSpikeTimes = data.spikeTimes_sec;
    spike_Vect = zeros(1, numel(data.sfree_data));
    spike_Vect(round(currSpikeTimes*SR)) = 1;


    tWidth = 2;
    tOffset = 1;    
    [spike_atMvmt_trials]  =  getTrialsMatrix(Evt.mvmt_onset, spike_Vect, 1/SR, tWidth, tOffset);
    spikeTrialSum = sum(spike_atMvmt_trials,2);
    nSpikes(f)       = round(sum(spikeTrialSum) / tWidth * WIN_GET_SPIKENR);   
    
    MTs(f) = nanmean(Evt.MT);
    min_MTs(f) = min(Evt.MT);
    max_MTs(f) = max(Evt.MT);
end
max_MTs
nSpikes