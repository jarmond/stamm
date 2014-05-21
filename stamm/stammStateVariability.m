function [fractions]=stammStateVariability(beta,foldfactor)
% STAMMSTATEVARIABILITY Calculate fraction of changing between states
%
%  [FRACTIONS]=STAMMSTATEVARIABILITY(BETA,FOLDFACTOR) Calculates fold change
%  in BETA between successive states and reports the fraction of genes
%  whose fold change exceeds FOLDFACTOR.

if nargin<2
    foldfactor=2;
end

[n k]=size(beta);

fractions=zeros(1,k-1);
for i=1:k-1
    % Fold change between state i and i+1.
    change=beta(:,i+1)./beta(:,i);
    % Count number exceeding FOLDFACTOR, in either direction.
    fractions(i)=sum(change>foldfactor | change<(1/foldfactor))/n;
end
