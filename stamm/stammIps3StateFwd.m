function [S,P,pinf,lambda]=stammIps3StateFwd(w,beta,t)
% STAMMIPS3STATEFWD iPS reprogramming linear 3-state model, forward only
% Input: w = transition rate matrix (3x3)
%        beta = characteristic expr level for gene each state (nx3)
% Output: S = timecourse mean expression levels for genes (nxm)
%         P = fraction in state (3xm)
% where n=# genes, m=# timepoints
m=length(t);

% state probabilities

P=zeros(3,m);
P(3,:)=exp(-w(2,3)*t);
if abs(w(1,2)-w(2,3))<eps
    P(2,:)=0.;
else
    P(2,:)=(w(2,3)/(w(1,2)-w(2,3)))*(exp(-w(2,3)*t)-exp(-w(1,2)*t));
end
P(1,:)=1-sum(P);

% mean gene expression
S=beta*P;

if nargout>2
    pinf=[1; 0; 0];
    lambda=[1/w(2,3); 1/w(1,2)];
end
