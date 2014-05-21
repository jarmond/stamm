function [W,beta]=stammPlotConf(data,m,filepref,path,m1)
% STAMMPLOTCONF Plots fitted models with confidence intervals against experimental data
%
%    STAMMPLOTCONF(DATA,M,FILEPREF,PATH) Uses experimental DATA and sequence of M fit
%    datafiles with prefix FILEPREF in folder PATH to plot. The loaded file
%    is the concatenation of FILEPREF and a number from 1 to M with two
%    digits. PDF of plot is saved in PATH.
%
%    STAMMPLOTCONF(DATA,M,FILEPREF,PATH,M1) Limits sequence of files read to
%    start from M1 instead of 1.
%
%    [W,BETA] = STAMMPLOTCONF(DATA,M,FILEPREF,PATH) Returns all the W and BETA
%    matrixs read as cell arrays.

if nargin<5
    m1=1;
end

n=m-m1+1;
W=cell(n,1);
beta=cell(n,1);
ind=cell(n,1);
% Load each file.
for i=1:n
    tfile=[path '/' filepref num2str(i+m1-1,'%02d')];
    fprintf('Reading: %s\n',tfile);
    load(tfile);
    ind{i}=result.ind;
    W{i}=result.W;
    beta{i}=result.beta;
end

prefix=[path '/' filepref];
stammCompareModelConf(data,ind,beta,W,result.model,prefix);
