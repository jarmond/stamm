function stammOptimalPerturbation(data,model,reps,clusterFile,preps,perturbFactors,outdir,tuning)
% STAMMOPTIMALPERTURBATION Fit with perturbed data for range of factors
%
%    STAMMOPTIMALPERTURBATION(DATA,MODEL,REPS,CLUSTERFILE,PREPS,FACTORS,OUTDIR,TUNING)
%    Parameters are identical as for STAMMPERTURBATION except PERTURBFACTORS
%    is a vector of factors.

for i=1:length(perturbFactors)
    fprintf('Perturbation factor: %f\n',perturbFactors(i));
    
    perdir=[outdir '/f' num2str(i,'%02d')];
    mkdir(perdir);
    stammPerturbation(data,model,reps,clusterFile,preps,perdir,perturbFactors(i),tuning);
end

