function [S,P,pinf,l]=stammIps4StateFwdAna(w,beta,t)
% STAMMIPS4STATEFWDANA iPS reprogramming 4-state linear model, forward transitions only
% Analytical solution.
% Input: w = transition rate matrix (4x4)
%        beta = characteristic expr level for gene each state (nx4)
%          t = time points (1xm);
% Output: S = timecourse mean expression levels for genes (nxm)
%         P = fraction in state (4xm)
% where n=# genes, m=# timepoints
pinf=[];
l=[];
m=length(t);

P=zeros(4,m);
P(4,:)=exp(-w(3,4)*t);
lambda=[w(2,3)-w(3,4); w(1,2)-w(3,4); w(1,2)-w(2,3)];

P(3,:)=w(3,4)*(exp(-w(3,4)*t)-exp(-w(2,3)*t))/lambda(1);
P(2,:)=(w(2,3)*w(3,4)/lambda(1))*(((exp(-w(3,4)*t)-exp(-w(1,2)*t))/lambda(2))-((exp(-w(2,3)*t)-exp(-w(1,2)*t))/lambda(3)));
P(1,:)=1-P(2,:)-P(3,:)-P(4,:);

% mean gene expression
S=beta*P;
