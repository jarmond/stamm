% Script to reproduce results for "Dissecting stochastic transitions during
% cellular reprograming", Armond et al.

% jobs : 1 - fit S-T core genes
%        2 - fit S-T genome
%        3 - fit S-T with timepoint knockouts
%        4 - fit S-T with Gaussian perturbations
%        5 - fit M core genes
%        6 - analyse Fluidigm data (not currently used)
%        7 - find S-T marker pairs
%        8 - analyse Tang data
%        9 - optimal MAP penalty by cross validation
%       10 - optimal perturbation noise by correlation
%       11 - GO analysis

% g_s: a state gene expression signature set
% ind: index of rows of g_s in original dataset (data.g)
function REPRODUCE(jobs)
if nargin<1
    jobs=setdiff(1:10,6);
end

close all;
tic;
seed=472295; % Random integer from randi(1000000). Used as seed for paper.
stammRandomSeed(seed);

outdir='results';
mkdir(outdir);
fitdir=[outdir '/fits'];
mkdir(fitdir);
genomedir=[outdir '/genome'];
mkdir(genomedir);

% Load data
stammLoadDataMats; % Data mats can be generated from raw data by stammPrepareDataMats.m

clusterFile=[outdir '/cluster.mat'];
reps=100; % Take lowest error of 100 nonlinear least squares with random
          % starting points

mapPenalty=0.1;

if ismember(1,jobs)
    % SAMAVARCHI-TEHRANI et al.

    % Cluster core genes from Table S1 to select k_1 fitting genes.
    k1=7;
    stammClusterGenes(data,k1,clusterFile);

    % Fit models to S-T data.
    models={'stammIps2StateFwd','stammIps3StateFwd',...
            'stammIps4StateFwd','stammIps5StateFwd'};
    stammRunAllModels(models,'normal',data,reps,clusterFile,fitdir,mapPenalty);

    % Fit model with lambda=0 as comparison
    lambda0dir=[outdir '/lambda0'];
    mkdir(lambda0dir);
    stammFitCluster(data,'stammIps4StateFwd',reps,clusterFile,lambda0dir,0);

    % Generate plots.
    % State and expression plots.
    load([fitdir '/stammIps4StateFwd.mat']);
    stammPlotFit(data,[fitdir '/stammIps2StateFwd.mat'],fitdir);
    stammPlotFit(data,[fitdir '/stammIps3StateFwd.mat'],fitdir);
    stammPlotFit(data,[fitdir '/stammIps4StateFwd.mat'],fitdir);
    stammPlotFit(data,[fitdir '/stammIps5StateFwd.mat'],fitdir);
    stammPlotErrs(models,fitdir,fitdir,mapPenalty); % RSS error by model.
    stammPlotCondNumber(models,fitdir,fitdir); % Log reciprocal cond number by model.
    % Selection of genes from Table S1 for Fig. 3f
    sel=[59 55 67 4 64 62 8 7 2 36 68 75 39 51 19 20 17 3 24 30 35 50 79 77 47 81 38 54 82];
    stammBetaBars(data,[fitdir '/stammIps4StateFwd.mat'],fitdir,0,sel);
    stammBetaBars(data,[fitdir '/stammIps5StateFwd.mat'],fitdir,-1);
    % Plot Nanog detection threshold model.
    nanog=stammGeneToIndex('nanog',data.g_names);
    nanog_idx=find(result.ind==nanog);
    mean_nanog=mean(result.beta(nanog_idx,:));
    stammPlotThreshold(data,nanog_idx,[fitdir '/stammIps4StateFwd.mat'],...
                     [mean_nanog/2 mean_nanog 2*mean_nanog], fitdir, 1);

end

if ismember(2,jobs)
    % Fit filtered genome for four-state (takes some time).
    core_ind=data.gene_ind; % Keep core genes.
    % Load raw data for filtering by mean and variance.
    rawdata=stammPrepareSamavarchiTehraniData('../ips_data/wrana2010/mmc3.csv',0);
    genome_ind=stammFilterGenes(rawdata.g);
    data.gene_ind=genome_ind;
    stammFitCluster(data,'stammIps4StateFwd',reps,clusterFile,genomedir,mapPenalty);
    data.gene_ind=core_ind; % Restore core genes.

    load([genomedir '/stammIps4StateFwd.mat']);
    stammGenomePlot(result.beta,[genomedir '/heatmap.pdf']); % Genome-wide heatmap.
    close all;
    % Dump genome results.
    stammBetaToCsv(data,result.ind,result.beta,[genomedir '/TableS3.csv']);
    col=2; % column 2 in beta is 3rd state
    stammRankGenes(data,result.ind,result.beta,col,length(result.ind),'descend',...
                  [genomedir '/TableS2.csv']);

    % Ranking.
    upreg=[1334 360 507 3088 793]; % index of downregulated genes to plot as bars
    downreg=[745 3507 3251]; % index of downregulated genes to plot as bars
    stammPlotRankBars(rawdata,result.ind,result.beta,upreg,[genomedir '/upreg_ranks.pdf']);
    stammPlotRankBars(rawdata,result.ind,result.beta,downreg,[genomedir '/downreg_ranks.pdf']);

    clear rawdata;
end

if ismember(3,jobs)
    % Robustness checking - timepoint knockout.
    models={'stammIps4StateFwd'};
    %models={'stammIps2StateFwd','stammIps3StateFwd',...
    %        'stammIps4StateFwd','stammIps5StateFwd'};
    tkodir=[outdir '/tko'];
    mkdir(tkodir);
    stammRunAllModels(models,'timepointko',data,reps,clusterFile,tkodir,mapPenalty);
    %stammPlotConf(data,length(data.t),'stammIps2StateFwd-t',tkodir);
    %stammPlotConf(data,length(data.t),'stammIps3StateFwd-t',tkodir);
    stammPlotConf(data,length(data.t),'stammIps4StateFwd-t',tkodir);
    %stammPlotConf(data,length(data.t),'stammIps5StateFwd-t',tkodir);
end

if ismember(4,jobs)
    % Robustness checking - perturbation by Gaussian noise.
    %models={'stammIps2StateFwd','stammIps3StateFwd',...
    %        'stammIps4StateFwd','stammIps5StateFwd'};
    models={ 'stammIps4StateFwd'};
    numPerturbations=10;
    perturbdir=[outdir '/perturb'];
    mkdir(perturbdir);
    perturbFactor=0.2;
    stammRunAllModels(models,'perturb',data,reps,clusterFile,numPerturbations,perturbdir,perturbFactor,mapPenalty);
    %stammPlotConf(data,numPerturbations,'stammIps2StateFwd-p',perturbdir);
    %stammPlotConf(data,numPerturbations,'stammIps3StateFwd-p',perturbdir);
    stammPlotConf(data,numPerturbations,'stammIps4StateFwd-p',perturbdir);
    %stammPlotConf(data,numPerturbations,'stammIps5StateFwd-p',perturbdir);
    stammPlotPerturbations(data,numPerturbations,data.gene_ind(1:9), ...
                          'stammIps4StateFwd-p',perturbdir);
end

if ismember(5,jobs)
    % MIKKELSEN et al.
    mikdir=[outdir '/mik_fits'];
    mkdir(mikdir);
    mapFile=[mikdir '/wrana_mik_map.mat'];
    if exist(mapFile)==0
        [mik.map,mik.rmap]=stammGenerateMap(data,mik,mapFile);
    else
        load(mapFile);
        mik.map=map;
        mik.rmap=rmap;
    end
    % Cluster core genes from Table S1 to select k_1 fitting genes.
    matched=intersect(data.gene_ind,mik.rmap);
    mik.gene_ind=mik.map(matched);
    k1=7;
    mikClusterFile=[mikdir '/mik_cluster.mat'];
    stammClusterGenes(mik,k1,mikClusterFile);

    % Fit models to Mikkelsen.
    models={'stammIps2StateFwd','stammIps3StateFwd',...
            'stammIps4StateFwd'}; % Too few timepoints for 5-state model.
    stammRunAllModels(models,'normal',mik,reps,mikClusterFile,mikdir,mapPenalty);

    % Generate plots.
    stammPlotFit(mik,[mikdir '/stammIps2StateFwd.mat'],mikdir);
    stammPlotFit(mik,[mikdir '/stammIps3StateFwd.mat'],mikdir);
    stammPlotFit(mik,[mikdir '/stammIps4StateFwd.mat'],mikdir);
    stammPlotErrs(models,mikdir,mikdir);
    stammPlotCondNumber(models,mikdir,mikdir);
    clear mdata;
end

% $$$ if ismember(6,jobs)
% $$$     % Fluidigm data.
% $$$     figmdir=[outdir '/figm'];
% $$$     mkdir(figmdir);
% $$$     mapFile=[figmdir '/wrana_figm_map.mat'];
% $$$     if exist(mapFile)==0
% $$$         % map genes between S-T and Fluidigm
% $$$         [figm.map,figm.rmap]=generateSTFigmMap(data,figm,mapFile);
% $$$     else
% $$$         load(mapFile);
% $$$         figm.map=map;
% $$$         figm.rmap=rmap;
% $$$     end
% $$$     figm.flatg=flattenFigmData(figm.g); % Lump cells from different days together
% $$$
% $$$     % Rank pairs of genes under model genome fit.
% $$$     [err,ind,w,g_s]=ipsLoad([genomedir '/stammIps4StateFwd.mat']);
% $$$     % choose only those genes that are fitted in genome fit
% $$$     [figm.ind,IA,IB]=intersect(figm.rmap,ind);
% $$$     figm.markers=findMarkers(g_s(IB,:),ind(IB));
% $$$     figm.markers.figmpairs=figm.map(figm.markers.pairsind); % index pairs into
% $$$                                                             % Fluidigm data
% $$$     plotMarkers(figm.markers.disc,figm.markers.figmpairs,figm.g,figm.names,figmdir);
% $$$     % count cells in quadrants
% $$$     [figm.counts,figm.states,figm.cellquad]=countMarkers(figm.markers.figmpairs, ...
% $$$                                                       figm.flatg);
% $$$     pvalueIters=1000000;
% $$$     pvalueFirstVals=50; % number of pairs to use to calculate top number of
% $$$                         % states
% $$$     figm.pvalue=pvalueMarkers(figm.markers.disc,figm.markers.figmpairs,...
% $$$                               figm.flatg,pvalueFirstVals,pvalueIters);
% $$$
% $$$     maWindow=50;
% $$$     randomRepeats=10000;
% $$$     statesMarkers(figm.markers.disc,figm.markers.figmpairs,figm.flatg,...
% $$$                   maWindow,randomRepeats,[figmdir '/figm_avgStates.pdf']);
% $$$     save([figmdir '/figm_results.mat'],'figm');
% $$$     clear figm;
% $$$ end

if ismember(7,jobs)
    % Discriminatory marker pairs
    markerdir=[outdir '/markerdir'];
    mkdir(markerdir);
    load([genomedir '/stammIps4StateFwd.mat']);
    markers=stammFindMarkers(result.beta,result.ind);
    save([markerdir '/markers.mat'],'markers','result');
    stammMarkersToCsv(data,markers,result.beta,100000,[markerdir '/markers.csv']);
    clear markers;
end

if ismember(8,jobs)
    % Tang data.
    tangdir=[outdir '/tang'];
    mkdir(tangdir);
    mapFile=[tangdir '/wrana_tang_map.mat'];
    if exist(mapFile)==0
        % Map genes between S-T and Tang
        [tang.map,tang.rmap]=stammGenerateMap(data,tang,mapFile);
    else
        load(mapFile);
        tang.map=map;
        tang.rmap=rmap;
    end
    % Rank pairs of core genes under model fit.
    load([fitdir '/stammIps4StateFwd.mat']);
    % Choose only those genes that are present in model fit.
    [tang.ind,IA,IB]=intersect(tang.rmap,result.ind);
    tang.markers=stammFindMarkers(result.beta(IB,:),result.ind(IB));
    tang.markers.tangpairs=tang.map(tang.markers.pairsind); % Index into Tang data.

    [tang.counts,tang.states,tang.cellquad]=...
        stammCountMarkers(tang.markers.tangpairs,tang.all);

    % Calculate enrichment for number of states p-value against random ranking
    pvalueIters=1000000;
    pvalueFirstVals=50; % number of pairs to use to calculate top number of states
    tang.pvalue=stammPvalueMarkers(tang.markers.tangpairs,tang.all,...
                              pvalueFirstVals,pvalueIters);
    % P-value comes out as zero, hence p-value < 1/pvalueIters

    % Plot markers under ranking and randomized with discriminatory score.
    maWindow=50;
    randomRepeats=10000;
    stammStatesMarkers(tang.markers.disc,tang.markers.tangpairs,tang.all,...
                      maWindow,randomRepeats,[tangdir '/tang_avgStates.pdf']);

    % Produce bars for some example genes
    %goodSel=[1 2 9 11 14 20 25 28];
    %badSel=[3560 3563 3565 3570];
    %stammBarMarkerScores(tang.markers.disc,goodSel,badSel,[tangdir '/tang_scores.pdf']);

    save([tangdir '/tang_results.mat'],'tang');
end

if ismember(9,jobs)
    % Optimal MAP penalty cross-validation
    crossvaldir=[outdir '/cross_val'];
    mkdir(crossvaldir);

    load([fitdir '/stammIps4StateFwd.mat']);
    lambdas=[0.0 0.01 0.05 0.1 0.2 0.3];
    %validationTimepoint=4;
    core_ind=data.gene_ind;
    data.gene_ind=result.ind(1:7);
    %priorCrossValidation(data,'stammIps4StateFwd',reps,clusterFile,crossvaldir,lambdas,...
    %                     validationTimepoint);
    stammPriorSimpleCrossValidation(data,'stammIps4StateFwd',reps,clusterFile,...
                               crossvaldir,lambdas);
    data.gene_ind=core_ind;
end

if ismember(10,jobs)
    % Optimal perturbation noise
    optperturbdir=[outdir '/opt_perturb'];
    mkdir(optperturbdir);

    load([fitdir '/stammIps4StateFwd.mat']);
    core_ind=data.gene_ind;
    data.gene_ind=result.ind(1:7);
    numPerturbations=10;
    perturbFactors=[0.05 0.1 0.2 0.3 0.4];
    stammOptimalPerturbation(data,'stammIps4StateFwd',reps,clusterFile, ...
                            numPerturbations, perturbFactors,optperturbdir, ...
                            mapPenalty);

    stammPlotOptimalPerturbations(result.beta(1:7,:),numPerturbations, ...
                                'stammIps4StateFwd-p', optperturbdir,1, ...
                                perturbFactors);

    data.gene_ind=core_ind;
end

if ismember(11,jobs)
    % GO analysis
    godir=[outdir '/go'];
    go=geneont('File',[godir '/gene_ontology.obo']);
    mkdir(godir);

    transform=2; % -log10(pval)
    pvalThreshold=[0.01, 1e-6];
    distThreshold=1e-4;

    % UPREGULATED
    % Terms chosen by Kris
    stammGoAnalyseBingoData(godir,[godir '/kris'],'pvalThreshold',pvalThreshold,'filter',1,...
                           'goObject',go,'transform',transform);

    % Unsupervised term selection
    sorting=1; % kmeans
    termThreshold=[200 800];
    stammGoAnalyseBingoData(godir,[godir '/unsuper'],'pvalThreshold',pvalThreshold,'sort',...
                       sorting,'transform',transform,'termThreshold',...
                       termThreshold,'distThreshold',distThreshold,...
                       'goObject',go);

    % Categorised unsupervised
    stammGoAnalyseBingoData(godir,[godir '/unsuper_cat'],'pvalThreshold',pvalThreshold,'filter',...
                       4,'transform',transform,'goObject',go);


    % Categorised unsupervised, discard other category
    stammGoAnalyseBingoData(godir,[godir '/unsuper_cat_noother'],'pvalThreshold',pvalThreshold,...
                       'filter',5,'transform',transform,'goObject',go);


    % Top 50 by state
    n=50;
    stammGoAnalyseBingoData(godir,[godir '/top50'],'pvalThreshold',pvalThreshold,...
                       'topInState',n,'distThreshold',distThreshold,...
                       'goObject',go);


    % DOWNREGULATED
    % Unsupervised term selection
    sorting=1; % kmeans
    termThreshold=[200 800];
    stammGoAnalyseBingoData(godir,[godir '/unsuper_down'],'pvalThreshold',pvalThreshold,'sort',...
                       sorting,'transform',transform,'termThreshold',...
                       termThreshold,'distThreshold',distThreshold,...
                       'goObject',go,'readDownData',1);

    % Categorised unsupervised
    stammGoAnalyseBingoData(godir,[godir '/unsuper_cat_down'],'pvalThreshold',pvalThreshold,'filter',...
                       4,'transform',transform,'goObject',go,'readDownData',1);


    % Categorised unsupervised, discard other category
    stammGoAnalyseBingoData(godir,[godir '/unsuper_cat_noother_down'],'pvalThreshold',pvalThreshold,...
                       'filter',5,'transform',transform,'goObject',go,'readDownData',1);


    % Top 50 by state
    n=50;
    stammGoAnalyseBingoData(godir,[godir '/top50_down'],'pvalThreshold',pvalThreshold,...
                       'topInState',n,'distThreshold',distThreshold,...
                       'goObject',go,'readDownData',1);

end

if ismember(12,jobs)
    % Switch and pulse.
    rankdir=[outdir '/ranks'];
    mkdir(rankdir);

    S3switch={'Mad2l1'; 'Wdr5'; 'Jarid2'; 'Cenpa'; 'Ngfrap1';...
                    'Bex1'; 'Ezh2'; 'Nanog'; 'Gdf3'; 'Sall4'};
    S2switch={'Taf5'; 'Uchl5'; 'Cnnm3'; 'Tmpo'; 'Tert'; 'Rqcd1';...
                    'Pole2'; 'Ccnt1'; 'Klf8'; 'Kdm6a'};
    S2pulse={'Tgfbr2'; 'Tgfbi'; 'Twist1'; 'Tcf4'; 'Grk5'; 'Mmp3'; ...
                   'Jun'; 'Krt14'; 'Dkk2'; 'Tbx4'};
    S3pulse={'Dgka'; 'Ocln'; 'Cldn3'; 'Pkp3'; 'Limk2'; 'Epha1'; ...
                   'Efna4'; 'Mxi1'; 'Sox11'; 'Bmp8b'};
    S4pulse={'Mcm2'; 'Rad51'; 'Klf5'; 'Eed'; 'Gnl3'; 'Bub1'; ...
                    'Mycn'; 'Zfp42'; 'Dppa4'; 'Esrrb'; 'Utf1'; 'Sox2'};

    selNames={S2switch, S3switch, S2pulse, S3pulse, S4pulse};
    % Raw data for calculating coeff. of var., differential expr. and dynamic range.
    rawdata=stammPrepareSamavarchiTehraniData(['../ips_data/wrana2010/' ...
                        'mmc3.csv'],0);
    load([outdir '/genome/stammIps4StateFwd.mat']);
    stammSwitchAndPulse(rawdata,result.beta,result.ind,...
                                 rankdir,selNames);
end

toc;
