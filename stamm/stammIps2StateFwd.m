function [S,P,pinf,lambda]=stammIps2StateFwd(w,beta,t)
% STAMMIPS2STATEFWD iPS reprogramming 2-state model, forward transitions only
% Input: w = transition matrix (1x1) [w_21]
%        beta = characteristic expr level for gene each state (nx2)
% Output: S = timecourse mean expression levels for genes (nxT)
%         P = fraction in state (2xT)
% where n=# genes, T=# timepoints

T=length(t);

p0=1;
P=zeros(2,T);
P(2,:)=p0*exp(-w*t);
P(1,:)=1-P(2,:);

S=beta*P;

if nargout>2
    pinf=[1;0];
    lambda=1/w;
end
