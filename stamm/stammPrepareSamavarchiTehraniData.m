function data=stammPrepareSamavarchiTehraniData(csvfile,normal,selection)
% STAMMPREPARESAMAVARCHITEHRANIDATA Load Samavarchi-Tehrani data into environment.
%
%   DATA = STAMMPREPARESAMAVARCHITEHRANIDATA(NORMAL,SELECTION) Loads data from
%   Samavarchi et al., Cell Stem Cell 7, 64-67 (2010) and normalizes
%   according to NORMAL:
%       0 : no normalization
%       1 : per gene normalization (default)
%       2 : global normalalization
%   SELECTION specifies a core set of genes to use for initial
%   clustering. See stammGeneSelection for provided sets, defaults to 5.

if nargin<2
    normal=1;
end
if nargin<3
    selection=5;
end

% Read raw data from CSV.
fieldfmt='%d %d %q %q %q %f %f %f %f %f %f %f %f %f %f %*q %f %f %f %d %f %f %f %f %f %f %f %f %d %d %d %d %d %d %d %d %d %q';

fid=fopen(csvfile);
head1=fgetl(fid);
head2=fgetl(fid);
rawdata=textscan(fid,fieldfmt,'Delimiter',',');
fclose(fid);

% Times (days) are hardcoded from Samavarchi-Tehrani paper
data.t=[0 2 5 8 11 16 21 30];

% Columns of time course
data.g=[rawdata{7} rawdata{8} rawdata{9} rawdata{10} rawdata{11} rawdata{12} ...
        rawdata{13} rawdata{14}];

% Gene accessions
data.g_acc=arrayfun(@strtrim,rawdata{4});

% Gene names
data.g_names=arrayfun(@strtrim,rawdata{3});

% Primary ips
data.pips=rawdata{15};

% Primary MEF
data.pmef=rawdata{6};

% Store variances
data.g_var=var(data.g,0,2);

if normal==1
    fprintf('Per gene normalization\n');
    std_g=std(data.g,0,2);
    mean_g=mean(data.g,2);
    for i=1:length(std_g)
        data.g(i,:)=((data.g(i,:)-mean_g(i))/std_g(i)); % Normalise by mean and std.
    end
elseif normal==2
    fprintf('Global normalization\n');
    std_g=std(data.g(:));
    mean_g=mean(data.g(:));
    data.g=(data.g-mean_g)/std_g;
    
    data.pips=(data.pips-mean(data.pips))/std(data.pips);
    data.pmef=(data.pmef-mean(data.pmef))/std(data.pmef);
end

if selection>0
    data.genes=stammGeneSelection(selection);
    data.gene_ind=stammAccToIndex(data.genes,data.g_acc);
else
    data.genes={};
    data.gene_ind=(1:length(data.g_names))';
end

data.g_accint=stammAccessionToInteger(data.g_acc);
data.normal=normal;
fprintf('Gene set size: %d\n',length(data.gene_ind));
