function b=bic(T,m,s)
% BIC Bayesian Information Criterion
%
%   B = BIC(T,M,S) Calculate Bayesian Information Criterion assuming normally
%   distributed errors. T is number of observations, M is number of
%   parameters and S is the sum of square residuals.

b=-2*s-m*log(T);
