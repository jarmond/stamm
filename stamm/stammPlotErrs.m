function errs=stammPlotErrs(models,datadir,outdir,tuning,timepoints)
% STAMMPLOTERRS Plots the sum squared error of a set of models.
%
%   ERRS = STAMMPLOTERRS(MODELS,DATADIR) For each model in MODELS loads the
%   corresponding fit file in DATADIR and plots the sum squared error and
%   returns in ERRS.
%
%   STAMMPLOTERRS(MODELS,DATADIR,OUTDIR) Saves the plots as PDF in OUTDIR.
%
%   STAMMPLOTERRS(MODELS,DATADIR,OUTDIR,TUNING) Accounts for TUNING prior
%   parameter.
%
%   STAMMPLOTERRS(MODELS,DATADIR,OUTDIR,TUNING,TIMEPOINTS) Plot as error per TIMEPOINT.

if nargin<3
    outdir=[];
end
if nargin<4
    tuning=0;
end
if nargin<5
    timepoints=0;
end

n=length(models);
k=2:(n+1);

errs=zeros(n,1);

for i=1:n
    load([datadir '/' models{i} '.mat']);
    % Subtract off penalty term
    penalty=tuning*sum(result.beta(:));
    if timepoints==0
        errs(i)=result.totalErr-penalty;
    else
        errs(i)=(result.totalErr-penalty)/(size(result.beta,1)*timepoints);
    end
end

figure(1);
clf;
plot(k,errs,'-x');
xlabel('States');
ylabel('Residual sum of squares');

if ~isempty(outdir)
    outfile=[outdir '/errs'];
    saveas(gcf,[outfile '.pdf']);
    dlmwrite([outfile '.dat'],[k' errs],'delimiter',' ');
end
