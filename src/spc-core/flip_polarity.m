function flip_flag = flip_polarity(data, FS, cfg)
% this is to correct polarity of signal to compare phases across different
% electrode and brain regions. We use the HFA approach following the
% procedure in Fischer et al 2020 Elife
% HFA flipping procedure
% 1. Take raw data without any digital filters
% 2. High pass filter at F_hp (e.g., 300 Hz)
% 3. Full-wave rectification 
% 4. Low-pass filter at F_lp (e.g., 100 Hz)
% 5. Take your band-filtered signal and subdivide it in 126 equally spaced
%    phase bins (.05pi rad)
% 6. Average HFAn in each bin and smooth (20 samples) it
% 7. Compute average HFA in 0.4pi rad window centered around 0 rad
% 8. Compute average HFA in 0.2 pi at two extremes
% 9. If HFA (7) < HFA (8) => flip the signal


flip_flag = false;

foi = cfg.foi;
ordFfoi = cfg.ordFfoi;
nPhasebins = cfg.nPhasebins;
smooth_win = cfg.smooth_win;
phase_win = cfg.phase_win;
fig_flag = cfg.fig_flag;

% ATTENTION! IF NAN VALUES ----> EVERYTHING DOES NOT WORK ANYMORE
data = data(~isnan(data)); % rid off nan values

% Calculate HFA
FN = FS/2;
HFA = calc_HFA(data, FS, cfg);


% 5. Take your band-filtered signal and subdivide it in 126 equally spaced
%    phase bins (.05pi rad)
[bFOI, aFOI] = butter(ordFfoi, foi/FN);
dataFOI = filtfilt(bFOI, aFOI, data);
phaseFOI = angle(hilbert(dataFOI));


phaseEdges = linspace(-pi, pi, nPhasebins);
[~,~,loc] = histcounts(phaseFOI,phaseEdges);
dataFOI_bin = accumarray(loc(:), dataFOI(:))./accumarray(loc(:),1);
HFA_bin = accumarray(loc(:), HFA(:))./accumarray(loc(:),1);
phaseMids = .5*(phaseEdges(1:end-1) + phaseEdges(2:end));

% 6. Average HFAn in each bin and smooth (20 samples) it
HFA_bin_smoothed = zscore(smoothdata(HFA_bin,'movmean',smooth_win));

% 7. Compute average HFA in 0.4pi rad window centered around 0 rad
phaseMids_centers = phaseMids <= phase_win/2 & phaseMids >= -phase_win/2;
phaseMids_boundaries = (phaseMids >= -pi & phaseMids <= (-pi + phase_win/2)) | (phaseMids >= (pi - phase_win/2)  & phaseMids <= (pi - phase_win/2));
HFA_center = mean(HFA_bin_smoothed(phaseMids_centers));
HFA_boundaries = mean(HFA_bin_smoothed(phaseMids_boundaries));

if HFA_center < HFA_boundaries
    flip_flag = true;
end

% 8. Sanity check plot flipping

if fig_flag
    figure("renderer","painters","position",[400 400 400 400])
    plot(phaseMids,zscore(dataFOI_bin),"linewidth",1.5)
    hold on
    plot(phaseMids,HFA_bin_smoothed,"linewidth",1.5)
    xlabel(" FOI phase [rad] ")
    ylabel(" Amplitude [a.u.] ")
    
    if flip_flag
        %HFA_flipped = calc_HFA(-data_copy, FS, cfg);

        %dataFOI = filtfilt(bFOI, aFOI, -data_copy);
        phaseFOI = calc_shiftPhase(phaseFOI,pi);
        [~,~,loc] = histcounts(phaseFOI,phaseEdges);
        HFA_bin_flipped = accumarray(loc(:), HFA(:))./accumarray(loc(:),1);
        
        % 6. Average HFAn in each bin and smooth (20 samples) it
        HFA_bin_smoothed_flipped = zscore(smooth(HFA_bin_flipped, smooth_win));
        phaseMids = .5*(phaseEdges(1:end-1) + phaseEdges(2:end));
        plot(phaseMids,zscore(HFA_bin_smoothed_flipped),"linewidth",1.5)
        legend([" FOI data " ," HFA ", "Flipped"])
        title(sprintf("Flipped: T, Delta = %2.2f ", HFA_center - HFA_boundaries))
    else
        title(sprintf("Flipped: F, Delta = %2.2f ", HFA_center - HFA_boundaries))
        legend([" FOI data " ," HFA "])
    end
    

end

