function stammGoPlotBingoData(pval,names,outpref)
% STAMMGOPLOTBINGODATA Plot p-values from Bingo analysis as a heatmap
%
%   STAMMGOPLOTBINGODATA(PVAL,NAMES,OUTPREF) P-values in PVAL are represented as
%   a heatmap, with NAMES aside each row. Figure saved as PDF with prefix
%   OUTPREF.

[n k]=size(pval);

pvallims(1,:)=nanmax(pval);
pvallims(2,:)=nanmin(pval);

% Set nan's to small value so that they appear nicely on (bottom of) colour scale
pval(isnan(pval))=0;

plotNumber=120; % # of categories per plot
i=0;
a=1;
b=plotNumber;
while a<=n
    if b>n
        b=n;
    end
    i=i+1;
    figure(i);
    subplot(1,2,1);
    % Plot heatmap.
    imagesc(pval(a:b,:));
    hold on;
    for j=1:k
        nans=find(pval(a:b,j)==0); % nb: nans have been set to 0
        % Plot white cross in not significant categories.
        if ~isempty(nans)
            plot(j*ones(size(nans)),nans,'xw','markersize',3);
        end
    end
    hold off;
    % Write gene names.
    xl=get(gca,'xlim');
    tx=xl(2)+0.01*(xl(2)-xl(1));
    for j=a:b
        label=names{j};
        ht=text(tx,j-a+1,label,'horizontalalignment','left','fontsize',3);
    end
    set(gcf,'paperunits','centimeters');
    set(gcf,'papersize',[21 29.7]);
    set(gcf,'paperposition',[0.5 0.5 20 28]);
    set(gca,'fontsize',6);
    
    colorbar('location','WestOutside');
    yl=get(gca,'ylim');
    yw=yl(2)-yl(1);
    tx=xl(1);
    tw=xl(2)-xl(1);

    % Write min/max values at bottom/top of each column.
    for j=1:k
        ht=text(tx+((j-1.0)/k)*tw,yl(1)-0.05*yw,num2str(pvallims(1,j)),'fontsize',6);
        ht=text(tx+((j-1.0)/k)*tw,yl(2)+0.05*yw,num2str(pvallims(2,j)),'fontsize',6);
    end

    saveas(gcf,[outpref num2str(i,'_%02d') '.pdf']);
    a=a+plotNumber;
    b=b+plotNumber;
end
