function stammPlotPerturbations(data,m,ind,filepref,path)
% STAMMPLOTPERTURBATIONS Plots perturbations to data
%
%    STAMMPLOTFIT(DATA,M,IND,FILEPREF,PATH) Plots perturbations to experimental
%    DATA in sequence of M datafiles with prefix FILEPREF in folder PATH. The
%    loaded file is the concatenation of FILEPREF and a number from 1 to M with
%    two digits. PDF of plot is saved in PATH. Only genes in IND indexed into
%    DATA.G are ploted.

% Read each file.
perturbG=cell(m,1);
for i=1:m
    tfile=[path '/' filepref num2str(i,'%02d')];
    fprintf('Reading: %s\n',tfile);
    load(tfile);
    perturbG{i}=result.perturbG;
end

% Convert to 3d array.
perturbG=cell2mat(perturbG');
perturbG=reshape(perturbG,[size(perturbG,1) length(data.t) m]);

n=length(ind);
fig_n=ceil(sqrt(n));
fig_m=ceil(n/fig_n);
figure(1);
clf;
% Plot each gene.
for i=1:n
    subplot(fig_m,fig_n,i);
    plot(data.t,squeeze(perturbG(ind(i),:,:))');
    if mod(i-1,fig_n)==0
        ylabel('Expression');
    end
    if i>(fig_m-1)*fig_n
        xlabel('Time (days)');
    end
    xlim([min(data.t) max(data.t)]);
    title(data.g_names(ind(i)));
end

saveas(gcf,[path '/' filepref '_perturbations.pdf']);
