function idx=stammGoFilterBingoData(pval,distThreshold)
% STAMMGOFILTERBINGODATA Filter Bingo data on distance threshold
%
%   IDX = STAMMGOFILTERBINGODATA(PVAL,DISTTHRESHOLD) Filters Bingo p-values on
%   differences between states. Difference between PVAL for each category is
%   summed over each possible state pairing and compared to DISTTHRESHOLD. Index
%   of retained categories returned in IDX.

[n k]=size(pval);

% Calculate distance between states
pval(isnan(pval))=0;
D=zeros(n,1);
for ii=1:n
    for i=1:k
        for j=i+1:k
            D(ii)=D(ii)+abs(pval(ii,i)-pval(ii,j));
        end
    end
end

idx=(D>distThreshold);
