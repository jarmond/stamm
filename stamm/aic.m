% Akaike's Information Criterion
% assuming normally distributed erros
% T : observations
% m : parameters
% s : sum of square residual
function b=aic(T,m,s)

b=-2*s-2*m;
