function [map,rmap]=stammGenerateMap(data1,data2,outfile)
% STAMMGENERATEMAP Generate a gene index map from data1 to data2
%
%   [MAP,RMAP] = STAMMGENERATEMAP(DATA1,DATA2,OUTFILE) Maps gene indexes
%   between DATA1 and DATA2 using accession numbers. MAP gives the indexes of
%   each gene from DATA1 in DATA2 and RMAP gives the indexes of each gene
%   from DATA2 in DATA1. MAP and RMAP are saved in OUTFILE

[tf,map]=ismember(data1.g_accint,data2.g_accint);
[tf,rmap]=ismember(data2.g_accint,data1.g_accint);
save(outfile,'map','rmap');
