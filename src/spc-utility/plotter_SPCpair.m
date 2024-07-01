function fh = plotter_SPCpair(SPCpair, cfg)

% unfold cfg
Position  = cfg.Position;
EvtTypes = cfg.EvtTypes;

if isfield(cfg,'CLim')
    CLim = cfg.CLim;
else
    CLim = [];
end

if isfield(cfg,'CMap')
    CMap = cfg.CMap;
end

if isfield(cfg,'Visible')
    Visible = cfg.Visible;
else
    Visible = true;
end


nPairs = numel(SPCpair.PPCz);
fh = cell(1,nPairs);
for pair_i = 1 : nPairs
    EvtTimes = SPCpair.CenterEvts{1,pair_i}';

    fh{pair_i} = figure('renderer','painters','position',Position,'visible',Visible);
    if cfg.smooth
        imagesc(1:size(SPCpair.PPCz{1,pair_i},2),SPCpair.freqQ, imgaussfilt(SPCpair.PPCz{1,pair_i},[1 1]))
    else
        imagesc(1:size(SPCpair.PPCz{1,pair_i},2),SPCpair.freqQ,SPCpair.PPCz{1,pair_i})
    end
    hold on
    contour(1:size(SPCpair.PPCz{1,pair_i},2),SPCpair.freqQ,SPCpair.LabelMatrix{1,pair_i}, 1, 'linecolor','k', 'linewidth', 1)
    for evt = [2 3 4 5]'
        plot(EvtTimes(evt)' *[1,1], [min(SPCpair.freqQ) max(SPCpair.freqQ)], '--', 'Color', 'w', 'LineWidth', 2)
    end
    set(gca, 'XTick', EvtTimes([2 3 4 5]))
    xLabeling = EvtTypes(2:end-1);
    set(gca, 'XTickLabel', xLabeling)
    set(gca,'YDir','Normal')
    xtickangle(30)
 
    %set(gca, 'YTick', 1:1:size(SPCpair.PPCz{1,pair_i},1))
    %yticklabels(round(SPCpair.freqQ(1:end),2))
    set(gca,'YMinorTick','off')
    cb = colorbar;
    ylabel(cb, 'PPC [z-score]')
    ylabel('Frequency [Hz] ')
    set(gca,'YScale','log')
    colormap(CMap)
    if ~isempty(CLim)
        clim(CLim)
    end
    axis tight
    box off
    title(sprintf('SPC sub-%s (pair-%d): %s - %s (%s)', ...
        SPCpair.Subj{pair_i},SPCpair.pair_id(pair_i),SPCpair.Unit{pair_i}, ...
        SPCpair.ECoG{pair_i}, char(SPCpair.ECoG_atlas{pair_i})),'Interpreter','none');
%     figname = sprintf('SPC_sub-%s_pair-%d_%s_%s_[%s]', ...
%         SPCpair.Subj{pair_i},SPCpair.pair_id(pair_i),SPCpair.Unit{pair_i}, ...
%         SPCpair.ECoG{pair_i}, char(SPCpair.ECoG_atlas{pair_i}));
    
    %saveFigures(fh,fullfile(PATH_OUTPUT,figname));
end




