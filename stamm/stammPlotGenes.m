function stammPlotGenes(ind,t,g,g_names,outname)
% STAMMPLOTGENES Plot gene expression data
%
%   STAMMPLOTGENES(IND,T,G,G_NAMES) Plots gene expression data G for genes
%   indexed by IND, over timepoints T. Labels with names from G_NAMES, again
%   indexed by IND.
%
%   STAMMPLOTGENES(IND,T,G,G_NAMES,OUTNAME) Saves PDF to OUTNAME.

if nargin<5
    outname=[];
end

n=length(ind);
fig_n=ceil(sqrt(n));
fig_m=ceil(n/fig_n);

figure(1);
clf;
for i=1:n
    subplot(fig_m,fig_n,i);
    plot(t,g(ind(i),:),'.-');
    if mod(i-1,fig_n)==0
        ylabel('Expression');
    end
    if i>(fig_m-1)*fig_n
        xlabel('Time (days)');
    end
    xlim([min(t) max(t)]);
    title(g_names(ind(i)));
end
if ~isempty(outname)
  saveas(gcf,[outname '_expr.pdf']);
end
