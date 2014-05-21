function [names,score,ind]=stammRankGenes(data,ind,beta,k,n,direction,outfile)
% STAMMRANKGENES Find top n genes, ranked by relative expression in state k
%
%   [NAMES,SCORE,IND] = STAMMRANKGENES(DATA,IND,BETA,K,N,DIRECTION,OUTFILE)
%   Ranks top N genes by relative expression in state K, given by subset of BETA
%   indexed by IND, if DIRECTION is 'descend', or from bottom if DIRECTION is
%   'ascend'. Saves CSV to OUTFILE, with gene names from DATA.

if nargin<7
    outfile=[];
end

% Compute relative expression scores.
score=beta(:,k)./mean(beta,2);
[score,idx]=sort(score,1,direction);
ind=ind(idx);

if ~isempty(outfile)
    f=fopen(outfile,'w+t');

    fprintf(f,'Gene,Rel expr\n');
    for i=1:n
        fprintf(f,'%s,%.4f\n',data.g_names{ind(i)},score(i));
    end
    fclose(f);
end

score=score(1:n);
ind=ind(1:n);
names=data.g_names(ind);
