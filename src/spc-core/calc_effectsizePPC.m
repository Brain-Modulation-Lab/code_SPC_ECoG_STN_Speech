function ES = calc_effectsizePPC(PPC)
% Calculate effect size PPC following the paper Voloh 2020
%  ES = (1 + 2*sqrt(PPC))/(1 - 2*sqrt(PPC))

PPC(PPC< 0) = 0;
ES = (1 + 2*sqrt(PPC))./(1 - 2*sqrt(PPC));

end

