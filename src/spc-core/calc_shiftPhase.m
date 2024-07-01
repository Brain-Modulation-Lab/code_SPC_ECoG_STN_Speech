function phaseShift = calc_shiftPhase(phase, shift)
%CALC_SHIFTPHASE Summary of this function goes here
%   Detailed explanation goes here
phaseShift = angle(exp(1i*(phase + shift)));
end

