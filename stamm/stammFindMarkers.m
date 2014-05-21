function markers=stammFindMarkers(beta,ind)
% STAMMFINDMARKERS Test each gene pairing for suitability to identify states
%
%   MARKERS = STAMMFINDMARKERS(BETA,IND) Calculates a discriminatory score
%   between each pair of genes in subset of BETA. Returns MARKERS structure
%   with .PAIRS: two column matrix of pair indexes into BETA, .PAIRSIND: same
%   as .PAIRS but indexed into original dataset DATA.G, .DISC: discriminatory
%   scores.

% pairs is indexed into beta
% pairsind is indexed into dataset

[n k]=size(beta);
np=nchoosek(n,2);
fprintf('%d pairs\n',np);

disc=zeros(np,1,'single');
if n>2^16-1 
    error('Data size too small');
end
pairs=zeros(np,2,'uint16');

% Loop of each pair.
c=0;
for i=1:n
    for ii=i+1:n
        if mod(c,1000000)==0 && c>0
            fprintf('%d pairs checked\n',c);
        end

        % Calculate discriminatory score.
        D=0;
        for j=1:k
            for l=1:k
                D=D+abs(beta(i,j)-beta(ii,l));
            end
        end
        c=c+1;
        disc(c)=D;
        pairs(c,:)=[i ii];
    end
end

[markers.disc,idx]=sort(disc,'descend');
markers.pairs=pairs(idx,:);
markers.pairsind=zeros(size(markers.pairs),'uint16');
markers.pairsind=ind(markers.pairs);


