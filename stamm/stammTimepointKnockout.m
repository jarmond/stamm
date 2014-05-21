function stammTimepointKnockout(data,model,reps,clusterFile,outdir,tuning)
% STAMMTIMEPOINTKNOCKOUT Fit genes to perturbed data
%
%    STAMMTIMEPOINTKNOCKOUT(DATA,MODEL,REPS,CLUSTERFILE,PREPS,OUTDIR) Fits genes
%    indexed by DATA.GENE_IND using MODEL and returns RESULT structure. Number
%    of random starting points for fitting specified by integer
%    REPS. Representative genes to use for initial fitting are taken from
%    CLUSTERFILE. Fitting is repeated LENGTH(DATA.T) times, each time with the
%    ith timepoint removed. RESULTS structure is saved with MODEL name prefix in
%    OUTDIR.
%
%    STAMMTIMEPOINTKNOCKOUT(DATA,MODEL,REPS,CLUSTERFILE,OUTDIR,TUNING) TUNING
%    parameter applies a prior to the beta fit, penalizing excessively large
%    beta's.

if nargin<6
    tuning=0;
end

% Initial cluster representative genes for first fitting.
c=load(clusterFile);
result.ind0=c.ind;
cn=length(result.ind0);

% Genes for sequential fitting.
gene_other=setdiff(data.gene_ind,result.ind0);
modelId=stammModelToId(model);
m=length(data.t);

% Loop of number of timepoints.
for j=1:m
    % Full fit clustered genes
    fprintf('Timepoint knockout %d of %d:\n',j,m);
    fprintf('Fitting %d genes from clustering (%d remaining):\n',cn, ...
            length(gene_other));
    fprintf('%s ',data.g_names{result.ind0});
    fprintf('\n');
    ko_ind=[1:(j-1) (j+1):m];
    result.tko_g=data.g(:,ko_ind);
    result.tko_t=data.t(ko_ind);
    [beta,W,result.errs(1),summary]=stammFitStateModel(result.ind0,result.tko_g,result.tko_t,modelId,reps,false,[],[],tuning);
    fprintf('RSS = %f\n',result.errs(1));

    ind=[result.ind0; gene_other];
    for i=1:length(gene_other)
        % Fixed fit remaining genes
        fprintf('Fitting %s (%d remaining): ',data.g_names{gene_other(i)}, ...
                length(gene_other)-i);
        [beta,W,result.errs(i+1),summary]=stammFitStateModel(ind(1:(cn+i)),result.tko_g,result.tko_t, ...
                                                   modelId,reps,true,W,beta,tuning);
        fprintf('RSS = %f\n',result.errs(i+1));
    end
    
    outfile=[outdir '/' model '-t' num2str(j,'%02d') '.mat'];
    
    % Save results.
    result.tuning = tuning;
    % Make beta column order match Markov-chain order.
    result.beta = fliplr(beta);
    result.W = W;
    result.totalErr = sum(result.errs);
    result.intacc=data.g_accint(ind);
    result.j = j;
    result.model = model;
    result.ind = ind;
    save(outfile,'result');
end

