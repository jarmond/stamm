% Akaike's Information Criterion Corrected
% assuming normally distributed erros
% T : observations
% m : parameters
% s : sum of square residual
function b=aicc(T,m,s)

b=-2*s-T*(1+m/T)/(1-(m+2)/T);
