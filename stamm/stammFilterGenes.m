function ind=stammFilterGenes(g,percentile)
% STAMMFILTERGENES Filter gene data by average level, then by variance
%
%   IND = STAMMFILTERGENES(G) Filters out genes from G with mean lower than
%   25th percentile, and then those with variance lower than 25th
%   percentile. Returns indexes of surviving genes in IND.
%
%   STAMMFILTERGENES(G,PERCENTILE) Filters at PERCENTILE level.

if nargin<2
    percentile=25;
end

% Filter by average levels.
mu=mean(g,2);
cutoff=prctile(mu,percentile);
fprintf('Mean cutoff = %f\n',cutoff);
ind=find(mu>=cutoff);

% Filter by variance.
sigma=var(g(ind,:),0,2);
cutoff=prctile(sigma,percentile); 
fprintf('Variance cutoff = %f\n',cutoff);
ind=ind(sigma>=cutoff);

fprintf('%d genes remaining\n',length(ind));
