function [beta,W,best,summary]=stammFitStateModel(ind,g,t,modelId,reps,fixParams,W,...
                                                 beta,tuning)
% STAMMFITSTATEMODEL Fit model to a subset of gene expression timecourse
%
%   [BETA,W,BEST,SUMMARY] = STAMMFITSTATEMODEL(IND,G,T,MODELID,REPS) Fit model
%   to a subset of gene expression timecourse G indexed into by IND
%   (i.e. G(:,IND)). T is a vector of time points. MODELID is the model
%   identifying number as returned by modelToId. REPS specifies the number of
%   random starting points to use in fitting. Returns the fitted BETA vector and
%   W matrix, as well as squared error of the BEST fit and SUMMARY stats of
%   all the fits.
%
%   STAMMFITSTATEMODEL(IND,G,T,MODELID,REPS,FIXPARAMS,W,BETA) Setting FIXPARAMS
%   to 1 keeps W fixed and fits only the BETA's for the last entry in IND, which
%   are concatenated onto passed in BETA's.
%
%   STAMMFITSTATEMODEL(IND,G,T,MODELID,REPS,FIXPARAMS,W,BETA,TUNING) To run
%   with a prior constraint on BETA, set the TUNING parameter. See article
%   for details.

% Get model string and number of parameters required.
model=stammIdToModel(modelId);
[k,nw]=stammNumSP(modelId);
n=length(ind); % # of genes to fit.

if nargin<6
    fixParams=false;
end
if nargin<8 || isempty(W)
    if fixParams
        error('Must specify W and beta to use fixed mode');
    else
        x0=ones(n*k+nw,1);
    end
end
if nargin<9
    tuning=0
end

if fixParams
    x0=ones(1,k);
end

genes=g(ind,:);
opts=optimset('Display','off','TolFun',1e-6,'TolX',1e-6,'MaxIter',5000, ...
              'MaxFunEvals',50000);

% Fitting bounds.
if fixParams
    lb(1:k)=0.;
    ub(1:k)=1000.;
else
    lb(nw+1:n*k+nw)=0.; % min gene expr
    lb(1:nw)=0.; % min transition rate
    ub(nw+1:n*k+nw)=1000.; % max gene expr
    ub(1:nw)=5.; % max rate
end

rstd=5.; % Std. dev. for normal random starting points centred at x0.
best=inf;
best_x=[];
errs=zeros(reps,1);
for i=1:reps
    % Generate random starting point.
    x1=stammStartPointRandom(rstd,x0,modelId);
    % Fit with nonlinear least squares.
    if fixParams
        [x,errs(i),res,exitflag]=lsqnonlin(@wrap_statemodel_fix,x1,lb,ub, ...
                                           opts);
    else
        [x,errs(i),res,exitflag]=lsqnonlin(@wrap_statemodel,x1,lb,ub,opts);
    end
    if exitflag<1
        fprintf('Warning: exitflag was %d in rep %d\n',exitflag,i);
    elseif errs(i)<best
        % Record best fit.
        best=errs(i);
        best_x=x;
    end
end

if fixParams
    beta=[beta; best_x];
else
    % Unpack x.
    [W,beta]=stammUnpack(modelId,best_x,n);
end

if nargout>=4
    summary=[min(errs) max(errs) mean(errs) std(errs)];
end

%% NESTED FUNCTIONS
    function y=wrap_statemodel(x)
        % Run model using W and beta packed in x.
        [W,beta]=stammUnpack(modelId,x,n);
        S=eval([model '(W,beta,t)']);
        rss=log2(S)-genes;
        penalty=sqrt(tuning*sum(abs(beta),2)); % Prior mean for beta=1.
        y=[rss(:); penalty];
    end

    function y=wrap_statemodel_fix(x)
        % Run model using beta for new gene to fit.
        beta_new=x;
        S=eval([model '(W,beta_new,t)']);
        rss=log2(S)-genes(end,:);
        penalty=sqrt(tuning*sum(abs(beta_new)));
        y=[rss penalty];
    end
end
