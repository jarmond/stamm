function stammCompareModel(data,ind,beta,W,model,outname)
% STAMMCOMPAREMODEL Make plots of model compared to data. Called by stammPlotFit.

if nargin<6
  outname=[];
end

% Restore beta order to final state first.
beta=fliplr(beta);
% Evaluate model.
t_model=linspace(0,max(data.t),100);
[S,P]=eval([model '(W,beta,t_model)']);
S=log2(S);

% Make figure with subplot for each gene and plot expression timecourse.
n=length(ind);
k=size(beta,2);
fig_n=ceil(sqrt(n));
fig_m=ceil(n/fig_n);
figure(1);
clf;
for i=1:n
    subplot(fig_m,fig_n,i);
    plot(data.t,data.g(ind(i),:),'.-',t_model,S(i,:));
    if mod(i-1,fig_n)==0
        ylabel('Expression');
    end
    if i>(fig_m-1)*fig_n
        xlabel('Time (days)');
    end
    xlim([min(data.t) max(data.t)]);
    title(data.g_names(ind(i)));
end
if ~isempty(outname)
  saveas(gcf,[outname '_expr.pdf']);
end

% Make figure and plot model state probabilities.
figure(2);
clf;
plot(t_model,P);
graph_names=cellfun(@(x) num2str(x,'p_%d'),num2cell(1:k),'UniformOutput',false);
legend(graph_names);
xlabel('Time / days');
ylabel('Probability');
ylim([0 1]);
if ~isempty(outname)
  saveas(gcf,[outname '_states.pdf']);
end
