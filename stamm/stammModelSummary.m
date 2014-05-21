function stammModelSummary(model,sizet,beta,W,err,head,latex)
% STAMMMODELSUMMARY Output summary results of a model fit
%
%    Called by stammSummary.

if nargin<5
    head=true;
end
if nargin<6
    latex=false;
end

mid=stammModelToId(model);
[k,nw]=stammNumSP(mid);
x=stammPackParams(beta,W,mid);
[S,P,pinf,lambda]=eval([model '(W,0,0)']);

if head
    % Print header.
    fmtstr='%3s';
    vars={'n'};
    % Transition parameters.
    wn=wnames(mid);
    for i=1:nw
        fmtstr=[fmtstr ' %7s'];
        vars=[vars wn{i}];
    end
    % Asymptotic probabilities.
    for i=1:length(pinf)
        fmtstr=[fmtstr ' %7s'];
        vars=[vars num2str(i,'p%di')];
    end
    % Timescales.
    for i=1:length(lambda)
        fmtstr=[fmtstr ' %7s'];
        vars=[vars num2str(i,'t%d')];
    end
    % Misc. statistics.
    fmtstr=[fmtstr ' %3s %3s %7s %8s %8s'];
    vars=[vars {'m','mu','s^2','BIC','AIC'}];
    if latex
        fmtstr=latexify(fmtstr);
    end
    fprintf([fmtstr '\n'],vars{:});
end

% Print results.
fmtstr='%3d';
n=size(beta,1);
vars={n};
for i=1:nw
    fmtstr=[fmtstr ' %7.2g'];
    vars=[vars x(i)];
end

for i=1:length(pinf)
    fmtstr=[fmtstr ' %7.2g'];
    vars=[vars pinf(i)];
end
for i=1:length(lambda)
    fmtstr=[fmtstr ' %7.2g'];
    vars=[vars lambda(i)];
end

m=nw+k*n;
mu=sizet*n;
a=aic(mu,m,err);
b=bic(mu,m,err);
fmtstr=[fmtstr ' %3d %3d %7.3g %8.1f %8.1f'];
vars=[vars {m,mu,err,b,a}];

if latex
    fmtstr=latexify(fmtstr);
end
fprintf([fmtstr '\n'],vars{:});

%% NESTED FUNCTIONS
    function wn=wnames(mid)
    switch mid
      case 1
        wn={'w_1','w_2'};
      case 2
        wn={'w_12','w_21','w_23','w_32'};
      case 3
        wn={'w_13','w_31','w_23','w_32'};
      case 4
        wn={'w_12','w_21','w_24','w_34'};
      case 5
        wn={'w_12','w_21','w_24','w_34','w_42','w_43'};
      case 6
        wn={'w_12','w_21','w_23','w_32','w_34','w_43'};
      case 7
        wn={'w_12','w_21','w_23','w_24','w_32','w_42'};
      case 8
        wn={'w_12','w_21','w_23','w_32','w_34','w_43','w_45','w_54'};
      case 9
        wn={'w_12','w_13','w_21','w_23','w_31','w_32','w_34','w_43'};
      case 10
        wn={'w_12'};
      case 11
        wn={'w_12','w_23','w_34'};
      case 12
        wn={'w_12','w_23','w_34','w_45'};
      case 13
        wn={'w_12','w_23'};
      case 14
        wn={'w_13','w_23'};
      case 15
        wn={'w_12','w_24','w_32'};
    end
    end

    function s=latexify(s)
        s=strrep(s,' ',' & ');
        s=[s '\\\\'];
    end
end
