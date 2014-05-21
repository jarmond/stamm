function [ind,sumd]=stammClusterGenes(data,n,outfile)
% CLUSTERGENES k-means cluster genes and pick representatives from clusters
%
%   [IND,SUMD] = CLUSTERGENES(DATA,N,OUTFILE) Runs k-means on genes specified by
%   DATA.GENE_IND and returns index of a representative gene from each
%   cluster in IND. The sum of distances from each cluster member from their
%   centroid is returned in SUMD.

if nargin<3
    outfile=[];
end

t=data.t;
g=data.g;
sel=data.gene_ind;
g_names=data.g_names;

% Do the k-means.
[idx,c,sumd,D]=kmeans(g(sel,:),n,'replicates',1000);
sumd=sum(sumd);

% Plot cluster centroids.
figure(1);
clf;
plot(t,c);
graph_names=cellfun(@(x) num2str(x,'Centroid %d'),num2cell(1:n),'UniformOutput',false);
legend(graph_names);
xlabel('Time, t / days');
ylabel('log_2 expression');

% Print genes names in each cluster.
ind=zeros(n,1);
cind=cell(n,1);
cind_names=cell(n,1);
for i=1:n
    fprintf('Cluster %d:\n',i);
    cind{i}=sel(idx==i);
    nl=0;
    cind_names{i}=cell(length(cind{i}),1);
    for j=1:length(cind{i})
        cind_names{i}{j}=g_names{cind{i}(j)};
        fprintf('%s ',cind_names{i}{j});
        nl=nl+1;
        if nl>=12
            fprintf('\n');
            nl=0;
        end
    end
    fprintf('\n');
    % Pick gene in this cluster closest to centroid.
    [v,k]=min(D(idx==i,i));
    ind(i)=cind{i}(k);
end

% Print the cluster representatives.
rep_genes=cell(n,1);
fprintf('\nRepresentative cluster genes: ');
for i=1:n
    rep_genes{i}=g_names{ind(i)};
    fprintf('%s ',rep_genes{i});
end

% Save figure PDF and clusters in mat.
if ~isempty(outfile)
    save(outfile,'rep_genes','ind','cind','cind_names','sumd','idx','c','D');
    saveas(gcf,strrep(outfile,'mat','pdf'));
end
