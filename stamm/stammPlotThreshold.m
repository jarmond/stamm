function stammPlotThreshold(data,sel,filename,threshold,outdir,state)
% STAMMPLOTTHRESHOLD Plot fraction of cells exceeding a given expression threshold
%
%    STAMMPLOTTHRESHOLD(DATA,SEL,FILENAME,THRESHOLD) Plots predicted fraction
%    of cells exceeding expression values in vector THRESHOLD, for time range
%    in DATA. SEL selects the gene to compare and beta values are read from FILENAME.

if nargin<5
    outdir=[];
end

load(filename);

k=size(result.beta,2);
t_model=linspace(0,max(data.t),100);
[S,P]=eval([result.model '(result.W,result.beta,t_model)']);

Pt=zeros(length(threshold),length(t_model));
if nargin<6
    krange=1:k;
else
    krange=state;
end

% Expects beta columns with final state first.
bbeta=fliplr(result.beta);
for j=1:length(threshold)
    for i=krange
        Pt(j,:)=Pt(j,:)+P(i,:)*(1-expcdf(threshold(j),bbeta(sel,i)));
    end
end

figure(1);
clf;
plot(repmat(t_model,length(threshold),1)',Pt');
ylim([0 1]);
xlabel('Time');
ylabel('Fraction of cells over threshold');

if ~isempty(outdir)
    saveas(gcf,[outdir '/threshold.pdf']);
    dlmwrite([outdir '/threshold.dat'],[t_model' Pt'],'delimiter',' ');
end

