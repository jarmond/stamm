function stammPrepareDataMats(outdir)
% STAMMPREPAREDATAMATS Convert data sets into convenient mat files.

if nargin<1
    outdir='.';
end

fprintf('Samavarchi-Tehrani data...\n');
wrana=stammPrepareSamavarchiTehraniData('../ips_data/wrana2010/mmc3.csv',1,5);
save([outdir '/wrana_data.mat'],'-struct','wrana');

fprintf('Mikkelsen data...\n');
mik=stammPrepareMikkelsenData('../ips_data/mikkelsen2008/nature07056-s2.csv',1,5);
save([outdir '/mik_data.mat'],'-struct','mik');

fprintf('Tang data...\n');
tang=stammPrepareTangData('../ips_data/tang2010/tang2010.csv');
save([outdir '/tang_data.mat'],'-struct','tang');

fprintf('done.\n');
