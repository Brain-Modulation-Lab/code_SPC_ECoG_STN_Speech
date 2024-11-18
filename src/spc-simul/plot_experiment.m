function plot_experiment(SPCs, spc_factor, plot_title)
    % Plot results from the simulation sweeps
    %
    % INPUTS:
    %   SPCs        - cell array of SPC structs from different sweep conditions
    %   spc_factor  - array of SPC factors used in the simulation
    %   ES          - effect sizes for the SPC factors
    

    
    % Define legend entries
    legendEntries = arrayfun(@(x) sprintf('SPC: %.2f', x), spc_factor, 'UniformOutput', false);
    legendEntries{end+1} = "true SPC";
    legendEntries{end+1} = "Event1";
    legendEntries{end+1} = "Event2";

    % Create the figure
    figure("Position",[200 200 1000 800])
    tiledlayout(5, numel(SPCs));  % Define the grid layout

    % Iterate over the SPC conditions (SPC1, SPC2, etc.)
    for cond = 1:numel(SPCs)
        SPC = SPCs{cond};
        cfg = SPC.cfg;
        %title_text = sprintf('Condition %d: %.2fx firing rate mod', cond, spc_factor(cond));

        % Plot firing rate
        nexttile(cond);
        for ii = 1:numel(SPC)
            plot(SPC(ii).time_bin, SPC(ii).FR, '-'); hold on;
        end
        xline(median(cfg.event_time), '--k');
        xline(median(cfg.event_time2), '--m');
        xlabel("time [s]"); ylabel("Firing rate [Hz]");
        patch([cfg.coupling_window fliplr(cfg.coupling_window)], [0 0 40 40], 'r', 'FaceAlpha', 0.1);
        xlim([0.7 2.4]); title(plot_title{cond});

        % Plot PPC
        nexttile(cond + numel(SPCs), [2 1]);
        for ii = 1:numel(SPC)
            plot(SPC(ii).time_bin, SPC(ii).PPC_time, '-'); hold on;
        end
        xlabel("time [s]"); ylabel("PPC");
        patch([cfg.coupling_window fliplr(cfg.coupling_window)], [0 0 1 1], 'r', 'FaceAlpha', 0.1);
        xline(median(cfg.event_time), '--k');
        xline(median(cfg.event_time2), '--m');
        xlim([0.7 2.4]);
        ylim([0 0.13])

        % Plot PPC with WIN_METHOD
        nexttile(cond + numel(SPCs)*3, [2 1]);
        for ii = 1:numel(SPC)
            timeVec = median(SPC(ii).WIN_METHOD.winCenters - SPC(ii).WIN_METHOD.winCenters(:, 1)) + median(cfg.event_time2 - 0.8);
            plot(timeVec, SPC(ii).WIN_METHOD.PPC, '-'); hold on;

            if cfg.nperms > 0
                clust = bwconncomp(SPC(ii).WIN_METHOD.Sign);
                clust = clust.PixelIdxList;
                [~,idx_clust] = max(cellfun(@(x) numel(x),clust));

                plot(timeVec(clust{idx_clust}), repmat(0.15+ 0.02*(ii-1),1,numel(clust{idx_clust})), 'LineWidth',2,'color','k')      
            end

            
        end
        xlabel("time [s]"); ylabel("PPC win-method");
        patch([cfg.coupling_window fliplr(cfg.coupling_window)], [0 0 1 1], 'r', 'FaceAlpha', 0.1);
        xline(median(cfg.event_time), '--k');
        xline(median(cfg.event_time2), '--m');
        xlim([0.7 2.4]);
        ylim([0 0.13])


        
    end

    % Add legend to the final plot
    legend(legendEntries, "Location", "best");
end
