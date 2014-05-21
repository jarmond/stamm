function index=stammGeneToIndex(name,g_names)
% STAMMGENETOINDEX Identify dataset index for a gene name
% Return zero for not found.

for i=1:length(g_names)
    if strcmpi(name,g_names(i))==1
        index=i;
        return;
    end
end
index=0;
