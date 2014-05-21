function data=stammPrepareTangData(csvfile)
% STAMMPREPARETANGDATA Prepare Tang 2010 Cell Stem Cell data
    
fid=fopen(csvfile);

fieldfmt=['%q %q ' repmat('%d ',1,32)];
rdata=textscan(fid,fieldfmt,'Delimiter',',','HeaderLines',1);
fclose(fid);

data.g_acc=arrayfun(@strtrim,rdata{1});
data.g_names=arrayfun(@strtrim,rdata{2});
data.g_accint=stammAccessionToInteger(data.g_acc);

i=3;
% esc (embryonic stem cells)
n=13;
data.esc=cell2mat(rdata(i:i+n-1))';

i=i+n;
% day 3 outgrowth Oct4+
n=2;
data.day3=cell2mat(rdata(i:i+n-1))';

i=i+n;
% day 5 outgrowth Oct4+
n=5;
data.day5=cell2mat(rdata(i:i+n-1))';

i=i+n;
% e4-5 epiblast
n=3;
data.e45=cell2mat(rdata(i:i+n-1))';

i=i+n;
% icm (inner cell mass)
n=9;
data.icm=cell2mat(rdata(i:i+n-1))';

% collate data together
data.all=[data.esc; data.day3; data.day5; data.icm];

