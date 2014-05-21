function stammBetaBars(data,filename,outdir,sorting,sel,colorThresh)
% STAMMBETABARS Show beta matrix as series of horizontal bar plots
%
%   STAMMBETABARS(DATA,FILENAME) Plots beta matrix in FILENAME using gene
%   names from DATA.
%
%   STAMMBETABARS(DATA,FILENAME,OUTDIR) Saves PDF in OUTDIR.
%
%   STAMMBETABARS(DATA,FILENAME,OUTDIR,SORTING,SEL,COLORTHRESH) Optional
%   parameters to control plot. SORTING set to positive integer sorts
%   ascending by that column, set to any negative integer uses k-means
%   clustering to group rows of bars, and zero disables sorting. SEL allows a
%   subset of the fitted genes to be shown (indexed into
%   RESULT.IND). COLORTHRESH sets the value to switch between a green and red bar.

if nargin<3
    outdir=[];
end
if nargin<4
    sorting=0;
end
if nargin<5
    sel=[];
end
if nargin<6
    colorThresh=1; % std devs
end

load(filename);
[k,nw]=stammNumSP(stammModelToId(result.model));

% Select genes to display.
if ~isempty(sel)
    result.beta=result.beta(sel,:);
    result.ind=result.ind(sel);
end

% Rescale data about mean.
mbeta=mean(result.beta);
for i=1:size(result.beta,2)
     result.beta(:,i)=result.beta(:,i)/mbeta(i);
end

if sorting>0
    % Sort by column 'sorting'.
    [y,seq]=sort(result.beta(:,sorting),1,'ascend');
    result.beta=result.beta(seq,:);
elseif sorting<0
    % K-means clustering.
    idx=kmeans(result.beta,k,'replicates',1000,'distance','sqEuclidean');
    [y,seq]=sort(idx,1,'ascend');
    result.beta=result.beta(seq,:);
else
    seq=1:size(result.beta,1);
end
names={data.g_names{result.ind(seq)}};
n=size(result.beta,1);
x=1:n;

figure(1);
clf;
for i=1:k
    % Plot bars for column i of beta.
    subplot(1,k,i);
    hold on;
    idx=result.beta(:,i)<colorThresh;
    if sum(idx)>0
        barh(x(idx),result.beta(idx,i),'r');
    end
    idx=result.beta(:,i)>=colorThresh;
    if sum(idx)>0
        barh(x(idx),result.beta(idx,i),'g');
    end
    plot([1 1],[0 1000],':k');
    hold off;
    bnds=[0 max(result.beta(:))];
    xlim(bnds);
    set(gca,'ylim',[0.5 length(names)+.5],'ytick',[]);
    xlabel(['State ' num2str(i)]);
    if i==1
        xl=get(gca,'xlim');
        tx=xl(1)-0.05*(xl(2)-xl(1));
        for i=1:length(names)
            ht=text(tx,i,names{i},'horizontalalignment','right','fontsize',5);
        end
    end
    if i==k
        xl=get(gca,'xlim');
        tx=xl(2)+0.05*(xl(2)-xl(1));
        for i=1:length(names)
            ht=text(tx,i,names{i},'horizontalalignment','left','fontsize',5);
        end
    end
end

if ~isempty(outdir)
    saveas(gcf,[outdir '/' result.model '_bars.pdf']);
end
