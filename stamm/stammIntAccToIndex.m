function idx=stammIntAccToIndex(query,acc)
% STAMMINTACCTOINDEX Convert integer accession number to dataset index

idx=zeros(length(query),1);
for i=1:length(query)
    f=find(acc == query(i), 1); % find first match
    if ~isempty(f)
        idx(i)=f;
    end
end
