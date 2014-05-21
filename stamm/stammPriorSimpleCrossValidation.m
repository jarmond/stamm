function optimality=stammPriorSimpleCrossValidation(data,model,reps,clusterFile,outdir,lambda)
% STAMMPRIORSIMPLECROSSVALIDATION Cross-validation to choose MAP penalization tuning parameter lambda
%
%   OPTIMALITY =
%   STAMMPRIORSIMPLECROSSVALIDATION(DATA,MODEL,REPS,CLUSTERFILE,OUTDIR,LAMBDA)
%   Fits cluster representative genes from CLUSTERFILE to DATA for given MODEL
%   with timepoint knockouts, repeating fits REPS times for each MAP parameter
%   LAMBDA, then validates against excluded timepoint. Saves results to
%   OUTDIR. Plots residual sum of squared error vs LAMBDA

% Initial cluster representative genes for first fitting.
c=load(clusterFile);
result.ind=c.ind;
cn=length(result.ind);

modelId=stammModelToId(model);
m=length(data.t);
predictors=1:m;
optimality=zeros(length(lambda),2); % cols: mean, std
trainingErr=zeros(length(lambda),2); % cols: mean, std

% Loop over each tuning parameter.
for a=1:length(lambda)
    l=lambda(a);
    fprintf('Cross-validating lambda = %f\n',l);
    
    for j=predictors
        % Full fit clustered genes.
        fprintf('Timepoint knockout %d of %d:\n',j,m);
        fprintf('Fitting %d genes from clustering:\n',cn);
        fprintf('%s ',data.g_names{result.ind});
        fprintf('\n');
        % Timepoint to knock-out.
        ko_ind=setdiff(predictors,j);
        result.tko_g=data.g(:,ko_ind);
        result.tko_t=data.t(ko_ind);
        
        outfile=[outdir '/' model '-t' num2str(j,'%02d') '-a' num2str(a,'%02d') ...
                 '.mat'];

        [result.beta,result.W,result.err]=stammFitStateModel(result.ind,result.tko_g,result.tko_t,...
                                                          modelId,reps,false,[],[],l);
        fprintf('RSS = %f\n',result.err);

        % Make beta column order match Markov-chain order.
        result.beta = fliplr(result.beta);
        result.model = model;
        result.lambda = l;
        result.j = j;
        save(outfile,'result');
    end
    
    % Calculate average RSS over predictors
    rss=[];
    terr=[];
    for j=predictors
        outfile=[outdir '/' model '-t' num2str(j,'%02d') '-a' num2str(a,'%02d') ...
                 '.mat'];
        load(outfile);
        S=eval([model '(result.W,result.beta,data.t(j))']);
        rss(end+1)=sum((log2(S)-data.g(result.ind,j)).^2);
        terr(end+1)=result.err;
    end
    optimality(a,:)=[mean(rss) std(rss)];
    trainingErr(a,:)=[mean(terr) std(terr)];
end

figure(1);
clf;
errorbar(lambda,optimality(:,1),optimality(:,2));
xlabel('Lambda');
ylabel('RSS of excluded timepoint');
saveas(gcf,[outdir '/' model '-crossval.pdf'])
save([outdir '/' model '-crossval.mat'],'lambda','optimality','trainingErr');

