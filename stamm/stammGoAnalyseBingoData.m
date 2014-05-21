function [terms,names]=stammGoAnalyseBingoData(datadir,outprefix,varargin)
% STAMMGOANALYSEBINGODATA Analyse GO results from Bingo in context of STAMM results
%
%   [TERMS,NAMES] = STAMMGOANALYSEBINGODATA(OUTPREFIX,...) 
%      Options specified as string followed by parameter:
%        'filter', level : filter GO terms to level 2 or level 3, or Kris' genes (1)
%        'sort', mode : sort by kmeans (1) or by state discrimination (2)
%        'pvalThreshold', [lowerThreshold, upperThreshold] :
%        'distThreshold', threshold :
%        'termThreshold', [lowerThreshold, upperThreshold] :
%        'goObject', goObj : pass in Matlab GO object
%        'topInState', n : plot only the top n in each state
%        'transform', mode : log(log(1/p)) (1:default) or log(1/p) (2)
%
%      Saves results with OUTPREFIX as a prefix.

% Set default values.
filter=0;
sorting=0;
transform=2;
pvalThreshold=[];
distThreshold=[];
termThreshold=[];
goObj=[];
topInState=[];
readDownData=0;

% Process options.
for i=1:2:length(varargin)
    name=varargin{i};
    value=varargin{i+1};
    switch name
      case 'filter'
        filter=value;
      case 'sort'
        sorting=value;
      case 'pvalThreshold'
        pvalThreshold=value;
      case 'distThreshold'
        distThreshold=value;
      case 'termThreshold'
        termThreshold=value;
      case 'goObject'
        go=value;
      case 'topInState'
        topInState=value;
      case 'transform'
        transform=value;
      case 'readDownData';
        readDownData=value;
      otherwise
        error(['Unknown option ' name]);
    end
end

% Read Bingo data.
bingo=stammPrepareBingoData(readDownData,datadir);

terms=union(bingo{1}.id,bingo{2}.id);
terms=union(terms,bingo{3}.id);
terms=union(terms,bingo{4}.id);

switch filter
  case 1
    kterms=stammGoKrisGenes();
    % Map kterms to terms, to preserve order avoid intersect.
    map=zeros(length(kterms),1);
    for i=1:length(kterms)
        map(i)=find(terms==kterms(i));
    end
    if nnz(map)~=length(kterms)
        error('Some Kris terms missing?');
    end

    % Select only the kterms.
    terms=terms(map);
  case 2
    go2ndlevel=[22610; 65007; 9758; 15976; 1906; 8283; 71840; 9987; 16265;...
                32502; 51234; 40007; 2376; 51179; 40011; 8152; 51704; 32501;...
                48519; 19740; 6794; 43473; 48518; 50789; 3; 22414; 50896; ...
                48511; 23052; 6791; 16032];
    terms=intersect(terms,go2ndlevel);
  case 3
    load('go/go3rdlevel.mat');
    terms=intersect(terms,go3rdlevel);
  case {4, 5}
    % Mario's genes from term thresholding.
    terms=stammGoMarioGenes(filter==5);
  otherwise
end

k=4;
n=length(terms);
pval=zeros(n,k);
names=cell(n,1);
for i=1:n
    for j=1:k
        ind=find(bingo{j}.id==terms(i));
        if ~isempty(ind)
            pval(i,j)=bingo{j}.pval(ind);
            names(i)=bingo{j}.names(ind);
        else
            pval(i,j)=nan;
        end
    end
    if sum(isnan(pval(i,:)))==k
        names{i}='missing';
    end
end

% Threshold p-values.
if ~isempty(pvalThreshold)
    pval(pval>pvalThreshold(1))=nan;
    pval(pval<pvalThreshold(2))=pvalThreshold(2);
end

if ~isempty(topInState)
    for i=1:k
        % Sort by p-value for state i.
        [y,idx]=sort(pval(:,i),1,'ascend');
        stateOutprefix=[outprefix '_s' num2str(i)];
        stammGoPlotBingoData(-log10(pval(idx(1:topInState),:)),names(idx(1:topInState)),...
                                                        stateOutprefix);
    end
    
    % Finish here.
    return
end

if ~isempty(distThreshold)
    % Filter out p-values that are too similar across states.
    idx=stammGoFilterBingoData(pval,distThreshold);
    pval=pval(idx,:);
    names=names(idx);
    terms=terms(idx);
end

if ~isempty(termThreshold)
    % Filter out GO terms which have too few or too many descendants.
    d=arrayfun(@(x) length(getdescendants(go,x,'exclude','true')),terms);
    idx=find(d>termThreshold(1) & d<termThreshold(2));
    pval=pval(idx,:);
    names=names(idx);
    terms=terms(idx);
end

if sorting>0
    idx=stammGoSortBingoData(pval,sorting);
    pval=pval(idx,:);
    names=names(idx);
    terms=terms(idx);
end

% Transform for better color scale.
switch transform
  case 1
    pval=log10(-log10(pval));
  case 2
    pval=-log10(pval);
end

if ~isempty(outprefix)
    stammGoPlotBingoData(pval,names,outprefix);
end

fprintf('%d categories\n',size(pval,1));
