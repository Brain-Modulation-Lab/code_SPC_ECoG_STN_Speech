function [r,delay, f, pval] = compute_delay(freqs,phase_mean)


phase_mean_unwrap = unwrap(phase_mean);

p = polyfit(freqs,phase_mean_unwrap*180/pi,1);
delay = p(1)/360;
phase0 = p(2);

if nargout > 3
    [r,pval] = corr(freqs, phase_mean_unwrap);
else
    r = corr(freqs, phase_mean_unwrap);
end

if nargout > 2
    f = figure('Position',[200 200 400 400],'visible','off');
    %cmap = multigradient(hex2rgb({'#d7191c'; '#fdae61'; '#ffef5c'}), [3 5 10]);
    scatter(freqs,phase_mean_unwrap*180/pi,25,'k','filled')
    hold on
    plot(freqs, phase0 + delay*freqs*360)
    title(sprintf('delay = %1.2f ms, R^{2} = %1.2f', delay*1000,r^2))
    xlabel(' SPC frequency [Hz] ')
    ylabel(' Unwrapped SPC phase [°]')
    box off



end