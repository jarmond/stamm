function result=stammFitCluster(data,model,reps,clusterFile,outdir,tuning)
% STAMMFITCLUSTER Fit genes in given set starting with those from clustering
%
%    RESULT = STAMMFITCLUSTER(DATA,MODEL,REPS,CLUSTERFILE,OUTDIR) Fits genes
%    indexed by DATA.GENE_IND using MODEL and returns RESULT structure. Number
%    of random starting points for fitting specified by integer
%    REPS. Representative genes to use for initial fitting are taken from
%    CLUSTERFILE. RESULTS structure is saved with MODEL name prefix in
%    OUTDIR.
%
%    RESULT = STAMMFITCLUSTER(DATA,MODEL,REPS,CLUSTERFILE,OUTDIR,TUNING)
%    TUNING parameter applies a prior to the beta fit, penalizing excessively
%    large beta's.
%
%    Models included with STAMM are for linear Markov Chains with forward-only
%    transistions for 2 to 5 states. If MODEL is empty that the 3 state model
%    is used. The names of these models are:
%
%          stammIps2StateFwd
%          stammIps3StateFwd
%          stammIps4StateFwd
%          stammIps5StateFwd

if isempty(model)
    model='stammIps3StateFwd';
end

% Initial cluster representative genes for first fitting.
c=load(clusterFile);
result.ind0=c.ind;
cn=length(result.ind0);

% Genes for sequential fitting.
gene_other=setdiff(data.gene_ind,result.ind0);
modelId=stammModelToId(model);

% Full fit clustered genes.
result.errs=zeros(length(gene_other)+1,1);
fprintf('Fitting %d genes from clustering (%d remaining):\n',cn, ...
        length(gene_other));
fprintf('%s ',data.g_names{result.ind0});
fprintf('\n');
[beta,W,result.errs(1),summary]=stammFitStateModel(result.ind0,data.g,data.t, modelId, ...
                                      reps,0,[],[],tuning);
fprintf('RSS = %f\n',result.errs(1));

result.ind=[result.ind0; gene_other];
for i=1:length(gene_other)
    % Fixed fit remaining genes
    fprintf('Fitting %s (%d remaining): ',data.g_names{gene_other(i)}, ...
            length(gene_other)-i);
    [beta,W,result.errs(i+1),summary]=stammFitStateModel(result.ind(1:(cn+i)), ...
                                                   data.g,data.t, modelId, ...
                                                   reps,1,W,beta,tuning);
    fprintf('RSS = %f\n',result.errs(i+1));
end

outfile=[outdir '/' model '.mat'];
err=sum(result.errs);

fprintf('Total RSS = %f\n',err);

% Save results.
result.model = model;
result.tuning = tuning;
% Make beta column order match Markov-chain order.
result.beta = fliplr(beta);
result.W = W;
result.totalErr = err;
result.intacc=data.g_accint(result.ind);
save(outfile,'result');
