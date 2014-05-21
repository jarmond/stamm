function n=stammMinGenes(k,kappa,T)
% STAMMMINGENES Calculate minimum number of genes for k-state linear model
%
%   STAMMMINGENES(K,KAPPA,T) Minimum number of genes to determine K-state linear
%   model, with T timepoints. KAPPA is factor for over-determining, set to 1
%   for true minimum.

n=2*kappa*(k-1)/(T-kappa*k);
