function rc=stammPlotCondNumber(models,datadir,outdir)
% STAMMPLOTCONDNUMBER Plots reciprocal condition number of beta matrix of a set of models.
%
%   RC = STAMMPLOTCONDNUMBER(MODELS,DATADIR) For each model in MODELS loads the
%   corresponding fit file in DATADIR and plots the reciprocal condition number
%   of beta matrix and returns in RC.
%
%   STAMMPLOTERRS(MODELS,DATADIR,OUTDIR) Saves the plots as PDF in OUTDIR.

if nargin<2
    outdir=[];
end

n=length(models);
k=2:(n+1);

rc=zeros(n,1);

for i=1:n
    load([datadir '/' models{i} '.mat']);
    rc(i)=1/(cond(result.beta));
end

figure(1);
clf;
semilogy(k,rc,'-x');
xlabel('States');
ylabel('Reciprocal condition number');

if ~isempty(outdir)
    outfile=[outdir '/rcond'];
    saveas(gcf,[outfile '.pdf']);
    dlmwrite([outfile '.dat'],[k' rc],'delimiter',' ');
end
