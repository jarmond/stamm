Model-based Analysis of State Transitions (MAST) v1.0
-------------------------------------------

The MAST software enables the identification cellular state transitions in
timecourse genomic data, and accompanies the article Armond et al., "A
stochastic model dissects cellular states and heterogeneity in transition
processes."

To use the software data should be provided in the following format:
  A structure DATA with fields:
    DATA.G         nxm matrix where each row is a gene and each column is a 
                   timepoint
    DATA.T         m-length row vector of timepoints
    DATA.G_ACC     gene accession numbers as cell array of strings
    DATA.G_ACCINT  gene accession numbers as integers
    DATA.G_NAMES   gene names as cell array of strings
    DATA.GENE_IND  index of rows in DATA.G of core genes to target initial 
                   fitting
See mastPrepareSamavarchiTehraniData for an example.

A basic analysis procedure uses the following steps:
1) mastClusterGenes to identify core expression profiles for initial fitting.
2) mastFitCluster to fit the core genes to data and then use those genes to fit
   the remaining genes.
3) mastPlotFit to display the fitted model against data.

For parameter information see 'help functionname'.

To run mastClusterGenes, the number of desired clusters must be specified. For
the dataset MAST was created to analysed we used 7 clusters. As is typical with
clustering, it is possible to provide a number of clusters suitable for all
datasets. The suggested procedure is use a set of numbers, e.g. n=4,5,6,7,8,9
and plot the sum of squared residuals (RSS) against n, and to choose n to give
a small RSS but avoiding diminishing returns.

Four models are provided with MAST and solve the Markov Chain for linear
forward-transition models with 2 to 5 states, named as follows:
  mastIps2StateFwd
  mastIps3StateFwd
  mastIps4StateFwd
  mastIps5StateFwd

New models can be developed, and for k states, n genes and m timepoints, must
simply take a W transition matrix, a BETA (n x k) matrix and a T timepoint
m-length vector, and return S, an (n x m) matrix of mean expression levels for
each gene, and P, a (k x m) matrix of state probabilities. The following
functions must also be updated to inform other functions about the new model:
  mastModelToId   Given a model name returns a unique model id.
  mastIdToModel   Given a model id returns a model name.
  mastPackParams  Packs W and BETA matrices into a vector.
  mastUnpack      Unpacks vector into W and BETA matrices.


To reproduce the results discussed in the article, a script called REPRODUCE is
provided. This also serves to illustrate the usage of the software.


A list of MAST functions, with brief description follows. More help is available
by typing 'help functionname', where functionname should be replaced with the
name of the function of interest.

mastAccToIndex                 - Identify index for a gene accession number.
mastAccessionToInteger         - Convert accession number to integer
mastBarMarkerScores            - Plot a bar chart of a selection of marker pair scores.
mastBetaBars                   - Show beta matrix as series of horizontal bar plots
mastBetaToCsv                  - Dumps beta to CSV file
mastCompareModel               - Make plots of model compared to data. Called by mastPlotFit.
mastCompareModelConf           - Compare model and data, with confidence intervals
mastComputeRanks               - Computes the rank order of values.
mastCountMarkers               - Count cells in each quadrant (Fluidigm or Tang data)
mastDumpFit                    - Dump model fit to textfiles for external plotting
mastFilterGenes                - Filter gene data by average level, then by variance
mastFindMarkers                - Test each gene pairing for suitability to identify states
mastFitCluster                 - Fit genes in given set starting with those from clustering
mastFitStateModel              - Fit model to a subset of gene expression timecourse
mastGeneSelection              - Return accession numbers of sets of core genes.
mastGeneToIndex                - Identify dataset index for a gene name
mastGenerateMap                - Generate a gene index map from data1 to data2
mastGenomePlot                 - Plot genome beta fit using clustergram
mastGoAnalyseBingoData         - Analyse GO results from Bingo in context of MAST results
mastGoFilterBingoData          - Filter Bingo data on distance threshold
mastGoKrisGenes                - List of interesting GO categories identified by Kris
mastGoMarioGenes               - List of interesting categories identified by Mario.
mastGoPlotBingoData            - Plot p-values from Bingo analysis as a heatmap
mastGoSortBingoData            - Sort data from Bingo analysis
mastIdToModel                  - Convert model ID to model name string.
mastIntAccToIndex              - Convert integer accession number to dataset index
mastIps2StateFwd               - iPS reprogramming 2-state model, forward transitions only
mastIps3StateFwd               - iPS reprogramming linear 3-state model, forward only
mastIps4StateFwd               - iPS reprogramming 4-state linear model, fwd transitions only
mastIps4StateFwdAna            - iPS reprogramming 4-state linear model, forward transitions only
mastIps5StateFwd               - iPS reprogramming 5-state linear model, fwd transistions only
mastMarkersToCsv               - Dumps beta to CSV file
mastMinGenes                   - Calculate minimum number of genes for k-state linear model
mastModelSummary               - Output summary results of a model fit
mastModelToId                  - Convert model string to numerical ID
mastNumSP                      - Number of parameters and states for given model
mastOptimalPerturbation        - Fit with perturbed data for range of factors
mastPackParams                 - Pack parameters into one dimensional array
mastPerturbation               - Fit genes to perturbed data
mastPlotBetaBars               - Plot beta from fitting as bars
mastPlotCondNumber             - Plots reciprocal condition number of beta matrix of a set of models.
mastPlotConf                   - Plots fitted models with confidence intervals against experimental data
mastPlotErrs                   - Plots the sum squared error of a set of models.
mastPlotFit                    - Plots fitted model against experimental data
mastPlotGenes                  - Plot gene expression data
mastPlotOptimalPerturbations   - Scatter standard betas against robustness checking betas
mastPlotPerturbations          - Plots perturbations to data
mastPlotThreshold              - Plot fraction of cells exceeding a given expression threshold
mastPrepareBingoData           - Load BINGO results from hardcoded paths for article
mastPrepareDataMats            - Convert data sets into convenient mat files.
mastPrepareMikkelsenData       - Load Mikkelsen data into environment.
mastPrepareSamavarchiTehraniData - Load Samavarchi-Tehrani data into environment.
mastPrepareTangData            - Prepare Tang 2010 Cell Stem Cell data
mastPrintBeta                  - Print out fitted beta table
mastPriorSimpleCrossValidation - Cross-validation to choose MAP penalization tuning parameter lambda
mastPvalueMarkers              - Generate random rankings to calculate p-value for ranking
mastRandomSeed                 - Set default stream to clock based seed, or specified integer.
mastRankFoldChanges            - Report fold change rank for top state genes
mastRankGenes                  - Find top n genes, ranked by relative expression in state k
mastReadBingoData              - Reads BINGO output file
mastStartPointRandom           - Generate random start point for optimization
mastStatesMarkers              - Plot number of states against marker rank
mastSummary                    - Output summary data from fitting
mastSwitchAndPulse             - Computes switch and pulse scores for a selection of genes
mastTimepointKnockout          - Fit genes to perturbed data
mastUnpack                     - Unpack parameters packed by mastPackParams into W and beta
