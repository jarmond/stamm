function [stateRank,foldRankUp,foldRankDown]=stammPlotRankBars(rawdata,ind,beta,sel,outname,invert)
% PLOTRANKBARS Plots the indiviual state rank and fold up/down rank for gene
%
%   PLOTRANKBARS(RAWDATA,IND,BETA,SEL) Plots the relative expression rank in
%   each state in BETA, and the differential expression fold change up and down
%   in RAWDATA as bars for selected genes SEL indexed into IND.
%
%   PLOTRANKBARS(RAWDATA,IND,BETA,SEL,OUTNAME) Saves plot as PDF in OUTNAME.
%
%   PLOTRANKBARS(RAWDATA,IND,BETA,SEL,OUTNAME,INVERT) Setting INVERT to 1
%   plots state bars with ranking reversed.
%
%   [STATERANK,FOLDRANKUP,FOLDRANKDOWN] = PLOTRANKBARS(RAWDATA,IND,BETA,SEL)
%   Returns the state ranking, DE fold change up and down rankings in
%   STATERANK, FOLDRANKUP and FOLDRANKDOWN, respectively.

if nargin<5
    outname=[];
end
if nargin<6
    invert=0;
end

% Calculate fold change rank.
Tend=size(rawdata.g,2);
Tstart=1;
fold=rawdata.g(ind,Tend)-rawdata.g(ind,Tstart); % log2 fold change.
foldRankUp=stammComputeRank(fold);

fold=rawdata.g(ind,Tend)-rawdata.g(ind,Tstart); % log2 fold change.
foldRankDown=stammComputeRank(fold,'ascend');

% Compute state specificity scores.
k=size(beta,2);
stateRank=zeros(size(beta));
for i=1:k
    score=beta(:,i)./mean(beta,2);
    [score,idx]=sort(score,1,'descend');
    % where in sorted ranking was top scoring rows
    [idx,stateRank(:,i)]=sort(idx,'ascend');
end

n=length(sel);
fig_n=ceil(sqrt(n));
fig_m=ceil(n/fig_n);

maxRank=max([foldRankUp; foldRankDown])
invert
figure(1);
clf;
for i=1:n
    subplot(fig_m,fig_n,i);
    if invert
        y=[maxRank-stateRank(sel(i),:) foldRankUp(sel(i)) foldRankDown(sel(i))];
    else
        y=[stateRank(sel(i),:) foldRankUp(sel(i)) foldRankDown(sel(i))];
    end
    bar(k:-1:-1,log10(y));
    title(rawdata.g_names(ind(sel(i))));
end

% Limit to selected genes.
stateRank=stateRank(sel,:);
foldRankUp=foldRankUp(sel);
foldRankDown=foldRankDown(sel);

if ~isempty(outname)
    saveas(gcf,outname);
    dlmwrite(strrep(outname,'pdf','csv'),[ind(sel) foldRankUp foldRankDown stateRank]);
end
