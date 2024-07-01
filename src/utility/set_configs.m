function cfg = set_configs(mode);

switch mode
    case 'default'

        cfg = struct();

        cfg.flip.Fhp = 300;
        cfg.flip.Flp = 100;
        cfg.flip.ordFhp = 4;
        cfg.flip.ordFlp = 4;
        cfg.flip.nPhasebins = 126;
        cfg.flip.smooth_win = 20;
        cfg.flip.phase_win = 0.4*pi;
        cfg.flip.foi = [60 80];
        cfg.flip.ordFfoi = 4;
        cfg.flip.fig_flag  = false; % set to true to plot fig
        cfg.flip.flip_check = true; % set to true to check flip sign

        cfg.plv.WIN_GET_SPIKENR  = 0.15; % [sec], if you do not have lots of trials and thus not so many spikes you will want to make this bigger
        cfg.plv.tWidth_avgSpike  = 2;  % [sec] window size and..
        cfg.plv.tOffset_avgSpike = 1;  % [sec] ...offset used to calculate the average spike number

        cfg.plv.NUM_PERMS = 80;       % number of permutations, set at least to 200 for analysis (500 ideal but it is very slow!)
        cfg.plv.PERM_WIN_LEN_SEC = 2; % this is used to permute the phase-signal, the exact number does not matter,
        % it only needs to be longer than the largest bin that you'll use to compute the PLV with


        cfg.plv.LENGTH_WINDOW = 0.5;  % this is to set NBINS_WINDOW of bins in a window LENGTH_WINDOW
        cfg.plv.NBINS_WINDOW = 10;

        cfg.plv.THR_NAN = 0; % skip S epochs with >= THR_NAN of nan sample in ECoG data
        cfg.plv.MIN_TRIALS = 10;
        cfg.plv.MAX_NAN = 40; % set max percentage of nan samples in ECoG (whole session);
        cfg.plv.freqsOfInt = [4:4:26 30:10:140; 8:4:30 40:10:150]';
        cfg.plv.gamma_env = false;

        cfg.stat.DEACTIVATE_MPC = false; %  deactivate multiple comparison correction if true
        cfg.stat.FOCUS_SIGTEST = false; % default false




        % plot configuration
        cfg.plot.EvtTypes = {'-0.75', 'Cue-Onset', 'Cue-Offset', 'Speech-Onset', 'Speech-Offset','+0.75'}; % insert your relevant events here


        white     = [1, 1, 1] * 1;
        blue      = [0, 0.5, 1];
        darkBlue  = [0, 98, 179] / 255;
        grey      = [1, 1, 1] * 0.3;
        red       = [220, 72, 6] / 255;

        cfg.plot.whiteBlue_CM = [linspace(white(1), blue(1), 80); ...
            linspace(white(2), blue(2), 80); ...
            linspace(white(3), blue(3), 80)]';

        cfg.plot.whiteBlack_CM = [linspace(white(1), 0, 80); ...
            linspace(white(2), 0, 80); ...
            linspace(white(3), 0, 80)]';

        cfg.plot.whiteGrey_CM = [linspace(white(1), grey(1), 80); ...
            linspace(white(2), grey(2), 80); ...
            linspace(white(3), grey(3), 80)]';

        cfg.plot.whiteDarkBlue_CM = [linspace(white(1), darkBlue(1), 80); ...
            linspace(white(2), darkBlue(2), 80); ...
            linspace(white(3), darkBlue(3), 80)]';

        cfg.plot.whiteRed_CM = [linspace(white(1), red(1), 80); ...
            linspace(white(2), red(2), 80); ...
            linspace(white(3), red(3), 80)]';

        cfg.plot.greyBlue_div_CM = [flipud(cfg.plot.whiteGrey_CM); cfg.plot.whiteDarkBlue_CM];

        cfg.plot.greyRed_div_CM = [flipud(cfg.plot.whiteGrey_CM); cfg.plot.whiteRed_CM];

        cfg.plot.cols_goTrials = [150, 150, 150; ...
            27, 153, 139; ...
            64, 55, 107; ...
            134, 123, 187] / 255;

        cfg.plot.var = 'PPC';
        cfg.typedata = 'default';
        cfg.locked = false;
        cfg.powertrim = false;

 case "default-5_30Hz"
        cfg.flip.Fhp = 300;
        cfg.flip.Flp = 100;
        cfg.flip.ordFhp = 4;
        cfg.flip.ordFlp = 4;
        cfg.flip.nPhasebins = 126;
        cfg.flip.smooth_win = 20;
        cfg.flip.phase_win = 0.4*pi;
        cfg.flip.foi = [60 80];
        cfg.flip.ordFfoi = 4;
        cfg.flip.fig_flag  = false; % set to true to plot fig
        cfg.flip.flip_check = true; % set to true to check flip sign

        cfg.plv.WIN_GET_SPIKENR  = 0.15; % [sec], if you do not have lots of trials and thus not so many spikes you will want to make this bigger
        cfg.plv.tWidth_avgSpike  = 2;  % [sec] window size and..
        cfg.plv.tOffset_avgSpike = 1;  % [sec] ...offset used to calculate the average spike number

        cfg.plv.NUM_PERMS = 80;% 80;       % number of permutations, set at least to 200
        cfg.plv.PERM_WIN_LEN_SEC = 2; % this is used to permute the phase-signal, the exact number does not matter,
        % it only needs to be longer than the largest bin that you'll use to compute the PLV with


        cfg.plv.LENGTH_WINDOW = 0.5;  % this is to set NBINS_WINDOW of bins in a window LENGTH_WINDOW
        cfg.plv.NBINS_WINDOW = 10;

        cfg.plv.THR_NAN = 0; % skip S epochs with >= THR_NAN of nan sample in ECoG data
        cfg.plv.MIN_TRIALS = 10;
        cfg.plv.MAX_NAN = 40; % set max percentage of nan samples in ECoG (whole session);
        cfg.plv.freqsOfInt = [4:2.5:34; 8:2.5:38]';
        cfg.plv.gamma_env = false;

        cfg.stat.DEACTIVATE_MPC = false; %  deactivate multiple comparison correction if true
        cfg.stat.FOCUS_SIGTEST = false; % default false




        % plot configuration

        cfg.plot.EvtTypes = {'-0.75', 'Cue-Onset', 'Cue-Offset', 'Speech-Onset', 'Speech-Offset','+0.75'};



        white     = [1, 1, 1] * 1;
        blue      = [0, 0.5, 1];
        darkBlue  = [0, 98, 179] / 255;
        grey      = [1, 1, 1] * 0.3;
        red       = [220, 72, 6] / 255;

        cfg.plot.whiteBlue_CM = [linspace(white(1), blue(1), 80); ...
            linspace(white(2), blue(2), 80); ...
            linspace(white(3), blue(3), 80)]';

        cfg.plot.whiteBlack_CM = [linspace(white(1), 0, 80); ...
            linspace(white(2), 0, 80); ...
            linspace(white(3), 0, 80)]';

        cfg.plot.whiteGrey_CM = [linspace(white(1), grey(1), 80); ...
            linspace(white(2), grey(2), 80); ...
            linspace(white(3), grey(3), 80)]';

        cfg.plot.whiteDarkBlue_CM = [linspace(white(1), darkBlue(1), 80); ...
            linspace(white(2), darkBlue(2), 80); ...
            linspace(white(3), darkBlue(3), 80)]';

        cfg.plot.whiteRed_CM = [linspace(white(1), red(1), 80); ...
            linspace(white(2), red(2), 80); ...
            linspace(white(3), red(3), 80)]';

        cfg.plot.greyBlue_div_CM = [flipud(cfg.plot.whiteGrey_CM); cfg.plot.whiteDarkBlue_CM];

        cfg.plot.greyRed_div_CM = [flipud(cfg.plot.whiteGrey_CM); cfg.plot.whiteRed_CM];

        cfg.plot.cols_goTrials = [150, 150, 150; ...
            27, 153, 139; ...
            64, 55, 107; ...
            134, 123, 187] / 255;

        cfg.plot.var = 'PPC';
        cfg.typedata = 'default';
        cfg.locked = false;
        cfg.powertrim = false;



end