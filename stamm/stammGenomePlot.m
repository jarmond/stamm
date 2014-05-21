function stammGenomePlot(beta,outfile)
% STAMMGENOMEPLOT Plot genome beta fit using clustergram
%
%    STAMMGENOMEPLOT(BETA,OUTFILE) Plots genome BETA fit as a heatmap and
%    saves as PDF in OUTFILE.

cgobj=clustergram(beta,'Cluster','column');
plot(cgobj);
saveas(gcf,outfile);
