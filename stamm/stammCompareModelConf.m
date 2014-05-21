function stammCompareModelConf(data,ind,beta,W,model,outname)
% STAMMCOMPAREMODELCONF Compare model and data, with confidence intervals
%
%   Called by STAMMPLOTCONF.

if nargin<6
  outname=[];
end

[k,nw]=stammNumSP(stammModelToId(model));
t_model=linspace(min(data.t),max(data.t),100);
t_model_g=linspace(min(data.t),max(data.t),size(data.g,2));
repl=size(ind,1);

n=size(beta{1},1); % number of genes
ind0=ind{1};

Sdist=zeros([n length(t_model_g) repl]);
Pdist=zeros([k length(t_model) repl]);
for i=1:repl
    beta{i}=fliplr(beta{i});
    [~,Pdist(:,:,i)]=eval([model '(W{i},beta{i},t_model)']);
    S=eval([model '(W{i},beta{i},t_model_g)']);
    Sdist(:,:,i)=log2(S);
end
Smean=mean(Sdist,3);
Sstd=std(Sdist,0,3);
Pmean=mean(Pdist,3);
Pstd=std(Pdist,0,3);

fig_n=ceil(sqrt(n));
fig_m=ceil(n/fig_n);
figure(1);
clf;
for i=1:n
    subplot(fig_m,fig_n,i);
    hold on
    plot(data.t,data.g(ind0(i),:),'.-');
    errorbar(data.t,Smean(i,:),Sstd(i,:),'r');
    
    if mod(i-1,fig_n)==0
        ylabel('Expression');
    end
    if i>(fig_m-1)*fig_n
        xlabel('Time (days)');
    end
    xlim([min(data.t) max(data.t)]);
    title(data.g_names(ind0(i)));
end
if ~isempty(outname)
  saveas(gcf,[outname '_expr.pdf']);
end

figure(2);
clf;
plot(t_model,Pmean);
hold on;
cind=lines;
for i=1:size(Pmean,1)
    plot(t_model,Pmean(i,:)+Pstd(i,:),':','color',cind(i,:))
    plot(t_model,Pmean(i,:)-Pstd(i,:),':','color',cind(i,:))
end
hold off

graph_names=cellfun(@(x) num2str(x,'p_%d'),num2cell(1:k),'UniformOutput',false);
legend(graph_names);
xlabel('Time (days)');
ylabel('Probability');
ylim([0 1]);
if ~isempty(outname)
  saveas(gcf,[outname '_states.pdf']);
  dlmwrite([outname '_states.dat'],[t_model' Pmean' Pstd'],'delimiter',' ');
end
