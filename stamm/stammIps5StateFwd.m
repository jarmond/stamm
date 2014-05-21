function [S,P,pinf,l]=stammIps5StateFwd(w,beta,t)
% STAMMIPS5STATEFWD iPS reprogramming 5-state linear model, fwd transistions only
%
% Input: w = transition rate matrix (kxk)
%        beta = characteristic expr level for gene each state (nxk)
%          t = time points (1xm);
% Output: S = timecourse mean expression levels for genes (nxm)
%         P = fraction in state (kxm)
% where n=# genes, m=# timepoints

pinf=[];
l=[];
m=length(t);

W=[0 w(1,2)  0       0       0;
   0 -w(1,2) w(2,3)  0       0;
   0 0       -w(2,3) w(3,4)  0;
   0 0       0       -w(3,4) w(4,5)
   0 0       0       0       -w(4,5)];

% initial conditions
P0=[0; 0; 0; 0; 1];

P=zeros(5,m);
for i=1:m
    P(:,i)=expm(W*t(i))*P0;
end

% mean gene expression
S=beta*P;
