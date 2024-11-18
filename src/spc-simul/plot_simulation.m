function plot_simulation(time_bin, FR, PPC_across_trials,SPC, cfg)
    % Plot firing rate and PPC results
    figure;
    tiledlayout(3, 1);

    % Plot firing rate
    nexttile;
    plot(time_bin, FR);
    xlabel('Time window (s)');
    ylabel('Firing rate [Hz]');
    title(sprintf('Firing Rate Over Time (Es = %1.3f; SPC = %1.3f)',cfg.ES, cfg.spc_factor));
    patch([cfg.coupling_window fliplr(cfg.coupling_window)], [0 0 cfg.fr_baseline*1.8 cfg.fr_baseline*1.8], 'r', 'FaceAlpha', 0.1);
    xlim([0.7 2.4]);
    box off
    % Plot PPC
    nexttile;
    bar(time_bin, PPC_across_trials);
    xlabel('Time window (s)');
    ylabel('PPC');
    patch([cfg.coupling_window fliplr(cfg.coupling_window)], [0 0 1.2*max(PPC_across_trials)*ones(1,2)], 'r', 'FaceAlpha', 0.1);
    xlim([0.7 2.4]);
    ylim([0 0.01])
    box off
    % Plot PPC win-method
    nexttile;
    timeVec = median(SPC.winCenters - SPC.winCenters(:, 1)) + median(cfg.event_time2 - 0.8);
    bar(timeVec, SPC.PPC)
    xlabel('Time window (s)');
    ylabel('PPC win-method');
    patch([cfg.coupling_window fliplr(cfg.coupling_window)], [0 0 1.2*max(SPC.PPC)*ones(1,2)], 'r', 'FaceAlpha', 0.1);
    xlim([0.7 2.4]);
    box off
     ylim([0 0.01])


end