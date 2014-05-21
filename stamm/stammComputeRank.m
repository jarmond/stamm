function ranks=stammComputeRank(values,direction)
% STAMMCOMPUTERANKS Computes the rank order of values.

if nargin<2
    direction='descend';
end

% Sort.
[~,idx]=sort(values,direction);
% Where in sorted ranking was top scoring rows?
[~,ranks]=sort(idx,'ascend');
