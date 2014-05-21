function stammSummary(data,filename,head,latex)
% STAMMSUMMARY Output summary data from fitting
%
%    STAMMSUMMARY(DATA,FILENAME) Output summary data from fitting in FILENAME
%    corresponding to DATA.
%
%    STAMMSUMMARY(DATA,FILENAME,HEAD,LATEX) Optional parameters to control
%    output: HEAD set to 1 prints header, LATEX set to 1 prints in LaTeX format.

if nargin<3
    head=true;
end
if nargin<4
    latex=false;
end

load(filename);
minS=1e-12;
% when doing fixed param runs, the err in the output is only for the
% k new parameters each iterations, so recalculate the lsq residuals
% now
S=eval([result.model '(result.W,result.beta,data.t)']);
y=(log2(max(S,minS))-data.g(result.ind,:)).^2;
err=sum(y(:));

modelSummary(result.model,length(data.t),result.beta,result.W,err,head,latex);
