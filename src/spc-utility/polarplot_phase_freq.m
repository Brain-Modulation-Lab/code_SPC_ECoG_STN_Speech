function fh = polarplot_phase_freq(phase, radius, freqs, sz , colors)
freqs_ = unique(freqs);
phase_mean = arrayfun(@(x) circ_mean(phase(freqs  == x)), freqs_);
radius_mean = arrayfun(@(x) mean(radius(freqs  == x),'omitnan'), freqs_);


fh = figure('Position',[300 300 400 400]);
polarscatter(phase_mean,radius_mean,sz,freqs_,'filled')
hold on
polarplot(phase_mean,radius_mean,'k')
cb = colorbar;

if isempty(colors) ||  size(colors,1) == 1
    blueRamp = uint8(linspace(0, 255, numel(phase_mean)));
    cmap = zeros(numel(phase_mean),3,'uint8');
    cmap(:,2) = blueRamp;
elseif size(colors,3)

    cmap = multigradient(colors, [3 5 10]);

end
colormap(cmap) 
ylabel(cb, 'Frequency [Hz]')
thetaticks(0 : 90 : 360);
thetaticklabels({"0","\pi/2","\pi", "-\pi/2"})

