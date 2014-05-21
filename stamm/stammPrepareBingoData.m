function bingo=stammPrepareBingoData(readDownData,godir)
% STAMMPREPAREBINGODATA Load BINGO results from hardcoded paths for article
%
%   BINGO = STAMMPREPAREBINGODATA(READDOWNDATA) Reads up regulated data in
%   READDOWNDATA is 0, otherwise reads down regulated data. 

if nargin<1
    readDownData=0;
end

if readDownData==0
    files={'S1_bgo.csv';
           'S2_bgo.csv';
           'S3_bgo.csv';
           'S4_bgo.csv'};
    headers=[20; 23; 25; 22];
else
    files={'S1_bgo_down.csv';
           'S2_bgo_down.csv';
           'S3_bgo_down.csv';
           'S4_bgo_down.csv'};
    headers=[32; 24; 24; 28];
end    

n=length(headers);

s=cell(n,1);
for i=1:n
    s{i}=stammReadBingoData(fullfile(godir,files{i}),headers(i));
end

bingo=cell(n,1);
for i=1:n
    bingo{i}.id=s{i}{1};
    bingo{i}.pval=s{i}{3};
    bingo{i}.names=s{i}{8};
end
