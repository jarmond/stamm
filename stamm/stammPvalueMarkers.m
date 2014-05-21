function pvalue=stammPvalueMarkers(pairs,g,m,iters)
% STAMMPVALUEMARKERS Generate random rankings to calculate p-value for ranking
%
%   PVALUE = STAMMPVALUEMARKERS(PAIRS,G,M,ITERS) To calculate a p-value for
%   ranking for number of states under marker pair discriminatory score, ITERS
%   random permutations of the ranking are generated and the mean number of
%   states of the first M pairs is compared.

n=size(pairs,1);
% Count number of states in the input pairs.
[counts,states]=stammCountMarkers(pairs,g);

% Mean number of states of first m pairs.
ranked=mean(states(1:m));

geq=0;
for i=1:iters
    % Generate random permutation.
    randstates=states(randperm(length(states)));
    
    % Check against input pairs.
    random=(mean(randstates(1:m)));
    if random>=ranked
        geq=geq+1;
    end
end

% Calculate p-value.
pvalue=geq/iters;
