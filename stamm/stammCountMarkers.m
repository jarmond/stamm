function [counts,states,cellquad]=stammCountMarkers(pairs,g)
% STAMMCOUNTMARKERS Count cells in each quadrant (Fluidigm or Tang data)
%
%   [COUNTS,STATES,CELLQUAD] = STAMMCOUNTMARKERS(PAIRS,G) Counts number of cells
%   in each quadrant of pair expression G viewed as binary
%   (i.e. on/off). Returns COUNTS, number of STATES for each PAIR, and matrix
%   CELLQUAD indicating which cells are in which quadrant. PAIRS should be
%   indexed into G.

% Labelling: offoff 1, offon 2, onoff 3, onon 4.

epsilon=1; % Small cutoff to define zeros
k=4; % assume 4 states
m=size(pairs,1); % number marker pairs
n=size(g,1); % number of cells
counts=zeros(m,k); % # cells in each quadrant for pair
states=zeros(m,1); % # states represented in pair
cellquad=zeros(n,m);
for i=1:m
    % Logical indexing on quadrants.
    offoff=squeeze(g(:,pairs(i,1))<epsilon & g(:,pairs(i,2))<epsilon);
    offon=squeeze(g(:,pairs(i,1))<epsilon & g(:,pairs(i,2))>=epsilon);
    onoff=squeeze(g(:,pairs(i,1))>=epsilon & g(:,pairs(i,2))<epsilon);
    onon=squeeze(g(:,pairs(i,1))>=epsilon & g(:,pairs(i,2))>=epsilon);
    notnan=squeeze(~isnan(g(:,pairs(i,1))) & ~isnan(g(:,pairs(i,2))));
    counts(i,:)=[sum(offoff(:)) sum(offon(:)) sum(onoff(:)) sum(onon(:))];
    if sum(counts(i,:))~=sum(notnan(:))
        error('missing cells');
    end
    % Count number of states.
    states(i)=nnz(counts(i,:));
    
    % Record which quadrant each cell is in.
    cellquad(offoff,i)=1;
    cellquad(offon,i)=2;
    cellquad(onoff,i)=3;
    cellquad(onon,i)=4;
end
