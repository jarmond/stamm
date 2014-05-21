function x1=stammStartPointRandom(rstd,x0,modelId)
% STAMMSTARTPOINTRANDOM Generate random start point for optimization
%
%   X1 = STAMMSTARTPOINTRANDOM(RSTD,X0,MODELID) Random start point is
%   random normal centred on X0 with standard deviation RSTD. MODELID
%   specifies model type, for checking pathological cases.

epsilon=1e-3;
not_suitable=true;
while not_suitable
    % Generate proposed starting point.
    x1=max(epsilon*ones(size(x0)),x0+rstd*randn(size(x0)));
    not_suitable=false;
    
    % Check pathological cases (zero divisors).
    if modelId==4 && length(x1)>numSP(modelId) && abs(x1(1)+x1(2)-x1(3)-x1(4))<epsilon
        not_suitable=true;
    end
    if modelId==17 
        lambda=[x1(2)-x1(3); x1(1)-x1(3); x1(1)-x1(2)];
        if any(abs(lambda)<epsilon)
            not_suitable=true;
            % Sometimes discards valid starting points but doesn't matter.
        end
    end
end

        
