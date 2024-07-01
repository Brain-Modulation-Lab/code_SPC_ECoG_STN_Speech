function PPC = correct_PLV(PLV,n)
% ##########################################
% correct plv with ppc0 using the formula Vinck 2010:
%  PPC = n/(n-1)*(PLV^2 - 1/n)
% ##########################################
PPC = n/(n-1)*(PLV.^2 - 1/n);
