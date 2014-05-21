function index=stammAccToIndex(acc,accArray)
% STAMMACCTOINDEX  Identify index for a gene accession number.
%
%    INDEX = STAMMACCTOINDEX(ACC,ACCARRAY) Returns the INDEX of a gene accession
%    number ACC in cell array of accession numbers ACCARRAY.

if iscell(acc)
    index=zeros(size(acc));
    for i=1:length(accArray)
        % Remove trailing .n numbers
        sp=regexp(accArray{i},' /// ','split'); % some contain multiple
        for k=1:length(sp)
            for j=1:length(acc) % assumes acc contains unique entries
                sp2=regexp(acc{j},' /// ','split');
                for l=1:length(sp2) % acc might contain ///
                    if ~isempty(strfind(sp{k},sp2{l}))
                        index(j)=i;
                        break;
                    end
                end
            end
        end
    end
else
    for i=1:length(accArray)
        % Remove trailing .n numbers
        sp=regexp(accArray{i},' /// ','split'); % some contain multiple
        for k=1:length(sp)
            sp2=regexp(acc,' /// ','split');
            for l=1:length(sp2) % acc might contain ///
                if ~isempty(strfind(sp{k},sp2{l}))
                    index=i;
                    return;
                end
            end
        end
    end
    index=0;
end
