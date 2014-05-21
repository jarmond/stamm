function stammPlotFit(data,filename,outdir,maxn)
% STAMMPLOTFIT Plots fitted model against experimental data
%
%    STAMMPLOTFIT(DATA,FILENAME) Uses experimental DATA and fit data in
%    FILENAME to plot.
%
%    STAMMPLOTFIT(DATA,FILENAME,OUTDIR) Saves PDF to OUTDIR.
%
%    STAMMPLOTFIT(DATA,FILENAME,OUTDIR,MAXN) Limits plotting to first MAXN genes.

if nargin<3
    outdir=[];
end

load(filename);

if isempty(outdir)
    prefix=[];
else
    prefix=[outdir '/' result.model];
end

if nargin>=4
    stammCompareModel(data,result.ind(1:maxn),result.beta(1:maxn,:),result.W,...
                     result.model,prefix);
else
    stammCompareModel(data,result.ind,result.beta,result.W,result.model,prefix);
end
