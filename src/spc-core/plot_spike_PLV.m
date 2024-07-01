function [] = plot_spike_PLV(PLV)


store = PLV.store;
cfg = PLV.cfg;
Const = cfg.plot;

E_channel = PLV.E_channel;
S_channel = PLV.S_channel;

NUM_PERMS = cfg.plv.NUM_PERMS;

DEACTIVATE_MPC = cfg.stat.DEACTIVATE_MPC; %  deactivate multiple comparison correction if true
FOCUS_SIGTEST = cfg.stat.FOCUS_SIGTEST; % default fals





if FOCUS_SIGTEST
    addit_focus_sigTest = 'focus_sigTest_';
else
    addit_focus_sigTest = '';
end



figure("renderer","painters","position",[400 400 600 500],"Visible","on");
sgtitle(strcat(E_channel,' - ', S_channel))

for sp = 1:4
    ax(sp) = subplot(2,2,sp);
    addit_BSL = '';
    if sp == 1
        REM_BASELINE = false; % no baseline normalization
        PLOT_PLV_OR_ZSCORE = 1;
    elseif sp == 2
        REM_BASELINE = false; % no baseline normalization
        PLOT_PLV_OR_ZSCORE = 2;
    elseif sp == 3
        REM_BASELINE = true; % yes baseline normalization
        PLOT_PLV_OR_ZSCORE = 1;
    elseif sp == 4
        REM_BASELINE = true; % yes baseline normalization
        PLOT_PLV_OR_ZSCORE = 2;
    end
    

    addit_medSplit = '';
    
    if size(store.PLV,1) > 1 % if we calculated it for several frequencies
        
        if strcmpi(cfg.plot.var,'PLV')
            origPLV =  store.PLV;
            PLV_perm = squeeze(store.PLV_perm);
        elseif strcmpi(cfg.plot.var,'PPC')
            origPLV =  store.PPC;
            PLV_perm = squeeze(store.PPC_perm);
        end
        
        if REM_BASELINE
            bsl_win = 1:PLV.centerEvts(2);
            baseline_PLV = mean(origPLV(:,bsl_win),2);
            origPLV = bsxfun(@minus, origPLV, baseline_PLV);
            
            baseline_PLV_perm = repmat(mean(PLV_perm(:,bsl_win,:),2),1,size(PLV_perm,2),1);
            PLV_perm = (PLV_perm - baseline_PLV_perm);
            addit_BSL = '_rem_baseline';
        end
        zscores_orig = (origPLV - mean(PLV_perm,3)) ./ std(PLV_perm,[],3);
        zscores_perm = bsxfun(@rdivide, bsxfun(@minus, PLV_perm, mean(PLV_perm,3)), std(PLV_perm,[],3));
    else
        % if I only calculated it for a single frequency
        origPLV =  store.PLV';
        PLV_perm = squeeze(store.PLV_perm)';
        zscores_orig = (origPLV' - mean(PLV_perm)) ./ std(PLV_perm);
        zscores_perm = bsxfun(@rdivide, bsxfun(@minus, PLV_perm, mean(PLV_perm)), std(PLV_perm));
    end
       p_perm =  2 * (1 - normcdf(abs(zscores_perm), 0, 1)); % get p-values from the zscore, abs to make it 2-tailed
    meanPerm = repmat(mean(zscores_perm,3), 1, 1, NUM_PERMS);
    p_orig = (sum(abs(zscores_perm - meanPerm) >= abs(zscores_orig - meanPerm), 3)+1) / (NUM_PERMS+1); %V2 from Ernst, Permutation Methods: A Basis for Exact Inference, 2004
    

    
    
    if FOCUS_SIGTEST
        freqNan = 1:size(p_orig,1);
        timeNan = 20:size(p_orig,2);
        p_orig(freqNan,timeNan)         = nan;
        zscores_orig(freqNan,timeNan)   = nan;
        zscores_orig(freqNan,timeNan,:) = nan;
        zscores_perm(freqNan,timeNan,:) = nan;
        origPLV(freqNan,timeNan,:)      = nan;
    end
    
    
    if DEACTIVATE_MPC
        clusPos_Z_Stat = p_orig < 0.05;
    else
        stats = getSignifClusters(p_orig, zscores_orig, p_perm, zscores_perm, 'zstat', 0.05, 0.05);
        clusPos_Z_Stat = stats.LabelMatrix;
    end
    
    if PLOT_PLV_OR_ZSCORE == 1
        plotMat = origPLV;
    elseif PLOT_PLV_OR_ZSCORE == 2
        plotMat = zscores_orig;
    end
    
    imagesc(plotMat); hold on
    cb = colorbar;
    addit_behav = '';

    if PLOT_PLV_OR_ZSCORE == 1 && ~REM_BASELINE

        colormap(ax(sp), Const.greyBlue_div_CM)

        cb.Label.String = [cfg.plot.var, addit_behav];
    elseif PLOT_PLV_OR_ZSCORE == 2 && ~REM_BASELINE
        colormap(ax(sp), Const.greyRed_div_CM)
        cb.Label.String = ['Z-scores', addit_behav];
    elseif REM_BASELINE
        colormap(ax(sp), Const.greyBlue_div_CM)
        if PLOT_PLV_OR_ZSCORE == 1
            cb.Label.String = [cfg.plot.var 'BSL subtracted', addit_behav];
        else
            cb.Label.String = [cfg.plot.var 'BSL subtracted Z-sc', addit_behav];
        end
    end
    
    axis xy
    if ~isnan(clusPos_Z_Stat)
        contour(clusPos_Z_Stat, 1, 'linecolor','k', 'linewidth', 1)
    end
    centerEvts = PLV.centerEvts;
    set(gca, 'XTick', centerEvts)
    xLabeling = Const.EvtTypes;
    set(gca, 'XTickLabel', xLabeling)
    
    plot(centerEvts(2)*[1,1], [0,size(plotMat,1)], '--', 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
    plot(centerEvts(3)*[1,1], [0,size(plotMat,1)], '--', 'Color', [0.5,0.5,0.5], 'LineWidth', 1.6)
    plot(centerEvts(4)*[1,1], [0,size(plotMat,1)], '--', 'Color', "k", 'LineWidth', 1)
    plot(centerEvts(5)*[1,1], [0,size(plotMat,1)], '--', 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
   
    yTicks = 1:2:size(PLV.freqsOfInt,1);
    set(gca, 'YTick', yTicks)
    newLabels    = mean(PLV.freqsOfInt(yTicks,:),2);
    set(gca, 'YTickLabel', newLabels)
    set(gca, 'FontSize', 8)
end

end
