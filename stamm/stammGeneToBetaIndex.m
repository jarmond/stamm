function index=stammGeneToBetaIndex(names,g_names,ind)
% STAMMGENETOBETAINDEX Identify (row) index in beta matrix for a set of gene names
% Returns zero for not found.

if ~iscell(names)
    names={names};
end

n=length(names);
index=zeros(n,1);
for i=1:n
    datasetindex=stammGeneToIndex(names{i},g_names);
    index(i)=find(ind==datasetindex);
end

