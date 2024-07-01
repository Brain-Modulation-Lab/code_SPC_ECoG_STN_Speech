function fh = plot_ClusterTime(cfg,Clusters);

T = cfg.T;
EvtTimes = cfg.EvtTimes;
bands = cfg.bands;
EvtTypes = cfg.cfg.plot.EvtTypes;
Unit_type = cfg.unit_type;

if ~isfield(cfg,'order')
    order = 'none';
else
    order = cfg.order;
end

if strcmpi(Unit_type,'all') || isempty(Unit_type)
    [ClusterCentroid,ClusterCentroidBand, ClusterOnOff] = get_SPCClustersProperties(cfg, Clusters);


    % order plot by frequency bands
    [~, idx_forder] = sort(ClusterCentroid(:,2),'descend');

    fh = figure('Position', [200 200 500 650],'Renderer','painters');
    tiledlayout(2,3)
    nexttile(1,[2 2])
    for clus_i = 1 : numel(ClusterCentroidBand)
        plot(ClusterOnOff(idx_forder(clus_i),:),clus_i*[1 1],'color',hex2rgb(bands.color(ClusterCentroidBand(idx_forder(clus_i)))))
        hold on
    end
    for evt = [2 3 4 5]'
        plot(EvtTimes(evt)' *[1,1], [0 numel(ClusterCentroidBand) + 1], '--', 'Color', 'k', 'LineWidth', 2)
    end
    set(gca, 'XTick', EvtTimes([2 3 4 5]))
    xLabeling = EvtTypes(2:end-1);
    set(gca, 'XTickLabel', xLabeling)
    set(gca,'YDir','Normal')
    xtickangle(30)
    xlim([T(1) T(end)])
    ylabel('t-SPC [id]')
    box off

    nexttile(3,[2 1])
    scatter(ClusterCentroid(idx_forder,2),1:numel(ClusterCentroidBand),5,'k','filled')
    box off
    xlabel('f_{t-SPC}')

else
    [ClusterCentroid,ClusterCentroidBand, ClusterOnOff,ClusterS_FRMod] = get_SPCClustersProperties(cfg, Clusters);

    for uoi = 1  : numel(Unit_type)

        switch Unit_type(uoi)
            case 'D'
                unit_type = -1;
                pair_list = find(ClusterS_FRMod == unit_type);

            case 'N'
                unit_type = 0;
                pair_list = find(ClusterS_FRMod == unit_type);

            case 'I'
                unit_type = 1;
                pair_list = find(ClusterS_FRMod == unit_type);

            case 'M'
                unit_type = 2;
                pair_list = find(ClusterS_FRMod == unit_type);
           

                % encoding
            case 'C'
                Unit_Info = get_EncodingTypes(Clusters,'consonant');
                pair_list = (Unit_Info == 1);
            case 'V'
                Unit_Info = get_EncodingTypes(Clusters,'vowel');
                pair_list = (Unit_Info == 1);
            case 'P'
                Unit_Info = get_EncodingTypes(Clusters,'position');
                pair_list = (Unit_Info == 1);

        end
        

        if size(pair_list,2) == 1
            pair_list = pair_list';
        end
        

        if strcmpi(order,'freq-wise')
            tmp = pair_list;
           [Units, idxUnits] = getUniquePairs(Clusters);
           idxUnitsOfInt = idxUnits(tmp);
           Frequency = ClusterCentroidBand(tmp);
           UniqueidxUnitsOfInt  =     unique(idxUnitsOfInt);
           modefrequency = arrayfun(@(x) mode(Frequency(idxUnitsOfInt == x)),UniqueidxUnitsOfInt);
           [~, orderfrequency] = sort(modefrequency,'descend');
           UniqueidxUnitsOfInt_sorted = UniqueidxUnitsOfInt(orderfrequency);

           pair_list = [];
           unit_list = [];
           for ui = 1 : numel(UniqueidxUnitsOfInt_sorted)
               pair_list = [pair_list find(idxUnits == UniqueidxUnitsOfInt_sorted(ui))'];
               unit_list = [unit_list [1 zeros(1, sum(idxUnits == UniqueidxUnitsOfInt_sorted(ui)) -1)]];
           end
           unit_list = find(unit_list == 1);
           ClusterCentroid_tmp = ClusterCentroid(pair_list,:);
           ClusterCentroidBand_tmp = ClusterCentroidBand(pair_list);
           ClusterOnOff_tmp = ClusterOnOff(pair_list,:);
        else

            % order plot by frequency bands
            [~, idx_forder] = sort(ClusterCentroid(pair_list,2),'descend');
            ClusterCentroid_tmp = ClusterCentroid(pair_list(idx_forder),:);
            ClusterCentroidBand_tmp = ClusterCentroidBand(pair_list(idx_forder));
            ClusterOnOff_tmp = ClusterOnOff(pair_list(idx_forder),:);
        end

        fh{uoi} = figure('Position', [200 200 500 650],'Renderer','painters');
        tiledlayout(2,3)
        nexttile(1,[2 2])
        for clus_i = 1 : numel(ClusterCentroidBand_tmp)
            plot(ClusterOnOff_tmp(clus_i,:),clus_i*[1 1],'color',hex2rgb(bands.color(ClusterCentroidBand_tmp(clus_i))))
            hold on
        end
        for evt = [2 3 4 5]'
            plot(EvtTimes(evt)' *[1,1], [0 numel(ClusterCentroidBand_tmp) + 1], '--', 'Color', 'k', 'LineWidth', 2)
        end
        set(gca, 'XTick', EvtTimes([2 3 4 5]))
        xLabeling = EvtTypes(2:end-1);
        set(gca, 'XTickLabel', xLabeling)
        set(gca,'YDir','Normal')
        xtickangle(30)
        xlim([T(1) T(end)])
        ylabel('t-SPC [id]')
        box off

        nexttile(3,[2 1])
        scatter(ClusterCentroid_tmp(:,2),1:numel(ClusterCentroidBand_tmp),5,'k','filled')
        if strcmpi(order,'freq-wise')
            hold on
            plot(repmat([0; 140],1,numel(UniqueidxUnitsOfInt_sorted)),[unit_list; unit_list],'--k')
        end
        box off
        xlabel('f_{t-SPC}')
    end
end