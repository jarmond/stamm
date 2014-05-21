function data=stammPrepareMikkelsenData(csvfile,normal,selection)
% STAMMPREPAREMIKKELSENDATA Load Mikkelsen data into environment.
%
%   DATA = STAMMPREPARESAMAVARCHITEHRANIDATA(NORMAL,SELECTION) Loads data from
%   Mikkelsen et al., Nature 454, 49-55 (2008) and normalizes
%   according to NORMAL:
%       0 : no normalization
%       1 : per gene normalization (default)
%       2 : global normalalization
%   SELECTION specifies a core set of genes to use for initial
%   clustering. See stammGeneSelection for provided sets, defaults to 5.

% Read data from CSV.
fid=fopen(csvfile);

fieldfmt=['%q %q %q %q ' repmat('%f ',1,18)];
rdata=textscan(fid,fieldfmt,'Delimiter',',','HeaderLines',1);
fclose(fid);

% Times (days).
data.t=[0 4 8 12 16];

% Gene data.
data.g=log2([rdata{7} rdata{8} rdata{9} rdata{10} rdata{11}]);

% Gene accessions.
data.g_acc=arrayfun(@strtrim,rdata{2});

% Gene names.
data.g_names=arrayfun(@strtrim,rdata{3});

% Normalise.
if normal
    std_g=std(data.g,0,2);
    mean_g=mean(data.g,2);
        for i=1:length(std_g)
        data.g(i,:)=((data.g(i,:)-mean_g(i,:))/std_g(i)); % normalise by mean and std
    end
end

% Some genes have std=0 because they have <20 expression throughout
% timecourse. Filter these out.
nanind=find(~isnan(data.g(:,1)));
data.g=data.g(nanind,:);
data.g_acc={data.g_acc{nanind}};
data.g_names={data.g_names{nanind}};
data.g_accint=stammAccessionToInteger(data.g_acc);
