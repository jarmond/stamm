function stammDumpFit(data,filename,outdir)
% STAMMDUMPFIT Dump model fit to textfiles for external plotting

load(filename);

t_model=linspace(0,max(data.t),100);
[S,P]=eval([result.model '(result.W,result.beta,t_model)']);
S=log2(S);

prefix=[outdir '/' result.model];

% Dump data
dlmwrite([prefix '_data.dat'],[data.t' data.g(result.ind,:)'],' ');

% Dump expr fit
dlmwrite([prefix '_exprfit.dat'],[t_model' S'],' ');

% Dump state fit
dlmwrite([prefix '_statefit.dat'],[t_model' P'],' ');
