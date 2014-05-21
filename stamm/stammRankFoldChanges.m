function [score,fold,rank,ind]=stammRankFoldChanges(beta,ind,rawG,k,n,dir)
% STAMMRANKFOLDCHANGES Report fold change rank for top state genes
%
%   [SCORE,FOLD,RANK]=STAMMRANKFOLDCHANGES(BETA,RAWG,K,N,DIR) Ranks top N each
%   gene in BETA by relative expression for state K, and also returns the FOLD
%   change of those genes based on raw gene data RAWG. DIR if specified can
%   be 'descend' or 'ascend' for ranking from top or bottom; default from
%   top. Returns relative expression SCORE, FOLD change, and fold change RANK.
%
%   [SCORE,FOLD,RANK,IND]=STAMMRANKFOLDCHANGES(BETA,RAWG,K,N,DIR,IND) Also
%   returns IND sorted corresponding to SCORE.

if nargin<6
    dir='descend';
end

% Sort by state specificity score first (report fold rank for genes in state
% score order).
score=beta(:,k)./mean(beta,2);
[score,idx]=sort(score,dir);
if nargin>=6
    ind=ind(idx);
end

Tend=size(rawG,2);
Tstart=1;

fold=rawG(ind,Tend)-rawG(ind,Tstart); % log2 fold change.

rank=stammRankGenes(abs(fold));

% Trim results.
score=score(1:n);
fold=fold(1:n);
rank=rank(1:n);
if nargin>=6
    ind=ind(1:n);
end
