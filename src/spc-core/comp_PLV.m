function  [PLV, meanPhase] = comp_PLV(phase)
% COMP_PLV computes the phase locking value
% INPUTS
% phase - vector of phases (or phase differences)
% OUTPUTS 
% PLV - Phase-locking value, which is a biased estimate depending on
%       N, if N is small, the PLV is inflated

PLV    = abs(nanmean(exp(1i*phase)));
 

if sum(isnan(phase(:))) > 0
    phase = (phase(~isnan(phase)));
end

if size(phase,2) == 1
   phase = phase'; 
end

meanPhase = circ_mean(phase, [], 2); % requires the circ_stats toolbox
% https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics

