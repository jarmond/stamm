function [idx,D,y]=stammGoSortBingoData(pval,mode)
% STAMMGOSORTBINGODATA Sort data from Bingo analysis
%
%   [IDX,D,Y] = STAMMGOSORTBINGODATA(PVAL,MODE) With MODE 1 sorts by k-means and
%   returns cluster distances in D. With MODE 2 sorts by distance between
%   states, and returns distance in D.  Returns IDX for sorting and sorted
%   p-values in Y.

[n k]=size(pval);

switch mode
case 1
  % sort by clustering
  pval(isnan(pval))=0;
  nc=16;
  [idx,C,D]=kmeans(pval,nc,'replicates',100,'distance','sqEuclidean');
  [y,idx]=sort(idx,1,'ascend');
  pval(pval==0)=nan;
case 2
  % sort by distance between states
  pval(isnan(pval))=0;
  D=zeros(n,1);
  for ii=1:n
      for i=1:k
          for j=i+1:k
              D(ii)=D(ii)+abs(pval(ii,i)-pval(ii,j));
          end
      end
  end
  pval(pval==0)=nan;
  [y,idx]=sort(D,1,'descend');
end  

