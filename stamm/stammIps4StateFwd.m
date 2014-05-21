function [S,P,pinf,l]=stammIps4StateFwd(w,beta,t)
% STAMMIPS4STATEFWD iPS reprogramming 4-state linear model, fwd transitions only
% 
% Input: w = transition rate matrix (4x4)
%        beta = characteristic expr level for gene each state (nx4)
%          t = time points (1xm);
% Output: S = timecourse mean expression levels for genes (nxm)
%         P = fraction in state (4xm)
% where n=# genes, m=# timepoints
pinf=[];
l=[];
m=length(t);

W=[0 w(1,2)  0       0;
   0 -w(1,2) w(2,3)  0;
   0 0       -w(2,3) w(3,4);
   0 0       0      -w(3,4)];

% initial conditions
P0=[0; 0; 0; 1];

P=zeros(4,m);
for i=1:m
    P(:,i)=expm(W*t(i))*P0;
end

% mean gene expression
S=beta*P;
