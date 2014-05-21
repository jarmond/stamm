function stammPerturbation(data,model,reps,clusterFile,preps,outdir,perturbFactor,tuning)
% STAMMPERTURBATION Fit genes to perturbed data
%
%    STAMMPERTURBATION(DATA,MODEL,REPS,CLUSTERFILE,PREPS,OUTDIR,PERTURBFACTOR)
%    Fits genes indexed by DATA.GENE_IND using MODEL and returns RESULT
%    structure. Number of random starting points for fitting specified by
%    integer REPS. Representative genes to use for initial fitting are taken
%    from CLUSTERFILE. Fitting is repeated PREPS times, each with independent
%    Gaussian perturbation on the data with s.d. PERTURBFACTOR. RESULTS
%    structure is saved with MODEL name prefix in OUTDIR.
%
%    STAMMPERTURBATION(DATA,MODEL,REPS,CLUSTERFILE,OUTDIR,PERTURBFACTOR,TUNING)
%    TUNING parameter applies a prior to the beta fit, penalizing excessively
%    large beta's.


if nargin<7
    perturbFactor=0.1;
end
if nargin<8
    tuning=0;
end

% Initial cluster representative genes for first fitting.
c=load(clusterFile);
result.ind0=c.ind;
cn=length(result.ind0);

% Genes for sequential fitting.
gene_other=setdiff(data.gene_ind,result.ind0);
modelId=stammModelToId(model);

% Loop of number of perturbations.
for j=1:preps
    % Full fit clustered genes.
    fprintf('Perturbation %d of %d:\n',j,preps);
    fprintf('Fitting %d genes from clustering (%d remaining):\n',cn, ...
            length(gene_other));
    fprintf('%s ',data.g_names{result.ind0});
    fprintf('\n');
    perturbG=randn(size(data.g))*perturbFactor+data.g;
    result.errs=zeros(length(data.gene_ind),1);
    [beta,W,result.errs(1),summary]=stammFitStateModel(result.ind0,perturbG,data.t,modelId,reps,false,[],[],tuning);
    fprintf('RSS = %f\n',result.errs(1));

    ind=[result.ind0; gene_other];
    for i=1:length(gene_other)
        % Fixed fit remaining genes
        fprintf('Fitting %s (%d remaining): ',data.g_names{gene_other(i)}, ...
                length(gene_other)-i);
        [beta,W,result.errs(i+1),summary]=stammFitStateModel(ind(1:(cn+i)),...
                                          perturbG,data.t,modelId,reps,true,W,beta,tuning);
        fprintf('RSS = %f\n',result.errs(i+1));
    end
    
    outfile=[outdir '/' model '-p' num2str(j,'%02d') '.mat'];
    
    % Save results.
    result.perturbG = perturbG;
    result.tuning = tuning;
    result.beta = beta;
    result.W = W;
    result.totalErr = sum(result.errs);
    result.intacc=data.g_accint(ind);
    result.model = model;
    result.ind = ind;
    save(outfile,'result');
end
