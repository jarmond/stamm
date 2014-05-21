function spScores=stammSwitchAndPulse(rawdata,beta,ind,outdir,selNames,cutoffLine,switches,pulses)
% STAMMSWITCHANDPULSE Computes switch and pulse scores for a selection of genes
%
%   SPSCORES = STAMMSWITCHANDPULSE(RAWDATA,BETA,IND) Calculates switch and pulse
%   scores for genes indexed by IND with BETA matrix and also calculates
%   coefficient of variation, differential expression and dynamic range
%   using RAWDATA. Specifically, by default, it computes S2, S3 switch and
%   S2, S3, S4 pulse.
%
%   STAMMSWITCHANDPULSE(RAWDATA,BETA,IND,OUTDIR) Saves data in CSV format in
%   OUTDIR with name switchandpulse.csv.
%
%   STAMMSWITCHANDPULSE(RAWDATA,BETA,IND,OUTDIR,SELNAMES) 
%   2D cell array SELNAMES specifies genes to scatter plot. Columns should
%   match up with computed switches followed by computed pulses. Each column 
%   is given a different symbol.
%
%   STAMMSWITCHANDPULSE(...,CUTOFFLINE) Draws an line on plot to emphasise
%   a cutoff. By default CUTOFFLINE=200.
%
%   STAMMSWITCHANDPULSE(...,SWITCHES,PULSES) The states for which switch and
%   pulse are calculated can be controlled by logical vectors SWITCHES and
%   PULSES. If BETA is NxK, then each vector must be K length. If SWITCHES is
%   specified, PULSE must be also.
    
if nargin<4
    outdir=[];
end
    
if nargin<5
    selNames=[];
end

if nargin<6
    cutoffLine=200;
end

if nargin<8
    switches=[0 1 1 0];
    pulses=[0 1 1 1];
end

[n k]=size(beta);
pulseScore=zeros(n,k);
switchScore=zeros(n,k);
coeffvar=zeros(n,1);
de=zeros(n,1);
dr=zeros(n,1);

% Compute scores.
for i=1:n
    % Current row.
    b=beta(i,:);
    mb=mean(b);
    for j=1:k
        if switches(j)
            s=[];
            for jj=j:k
                s=[s b(jj)./b(1:j-1)];
            end
            switchScore(i,j)=min(s);
            %a=min([b(2), b(3), b(4)]./b(1));
            %b=min([b(3)/b(1), b(3)/b(2), b(4)/b(1), b(4)/b(2)]);
        end
       
        if pulses(j)
            pulseScore(i,j)=b(j)/mb;
        end
    end
    
    g=rawdata.g(ind(i),:);
    coeffvar(i)=std(g)/mean(g);
    de(i)=g(end)-g(1);
    dr(i)=max(g)-min(g);
end

% Calculate switch and pulse ranks for requested states.
switchRank=zeros(size(switchScore));
pulseRank=zeros(size(pulseScore));
for j=1:k
    if switches(j)
        switchRank(:,j)=stammComputeRank(switchScore(:,j));
    end
    if pulses(j)
        pulseRank(:,j)=stammComputeRank(pulseScore(:,j));
    end
end

% Calculate cov, de, dr ranks.
coeffvarRank=stammComputeRank(coeffvar);
deRank=min([stammComputeRank(de) stammComputeRank(-de)],[],2);
drRank=stammComputeRank(dr);

% Write to CSV.
if ~isempty(outdir)
    % Write to csv.
    f=fopen(fullfile(outdir,'switchandpulse.csv'),'wt');
    % Write header.
    fprintf(f,'Name,');
    for i=1:k
        fprintf(f,num2str(i,'b_%d,'));
    end
    for i=1:k
        if switches(i)
            fprintf(f,num2str(i,'S%d switch rank,'));
        end
    end
    for i=1:k
        if pulses(i)
            fprintf(f,num2str(i,'S%d pulse rank,'));
        end
    end
    fprintf(f,'Cov rank,DE rank,DR rank,');
    for i=1:k
        if switches(i)
            fprintf(f,num2str(i,'S%d switch score,'));
        end
    end
    for i=1:k
        if pulses(i)
            fprintf(f,num2str(i,'S%d pulse score,'));
        end
    end
    fprintf(f,'Cov,DE,DR\n');
    
    % Write values.
    for i=1:n
        % Name.
        fprintf(f,'%s,',rawdata.g_names{ind(i)});
        % Beta.
        for j=1:k
            fprintf(f,'%.15f,',beta(i,j));
        end
        % Ranks.
        for j=1:k
            if switches(j)
                fprintf(f,'%d,',switchRank(i,j));
            end
        end
        for j=1:k
            if pulses(j)
                fprintf(f,'%d,',pulseRank(i,j));
            end
        end
        fprintf(f,'%d,%d,%d,',coeffvarRank(i),deRank(i),drRank(i));
        % Scores.
        for j=1:k
            if switches(j)
                fprintf(f,'%.15f,',switchScore(i,j));
            end
        end
        for j=1:k
            if pulses(j)
                fprintf(f,'%.15f,',pulseScore(i,j));
            end
        end
        fprintf(f,'%.15f,%.15f,%.15f\n',coeffvar(i),de(i),dr(i));
    end
    fclose(f);
end
    
if ~isempty(outdir) && ~isempty(selNames)
    % Make plots of rank vs coefficient of variation.
    l=[cutoffLine 0; cutoffLine n];
    % Find indices.
    ncols=size(selNames,2);
    idx=cell(ncols,1);
    for i=1:ncols
        idx{i}=stammGeneToBetaIndex(selNames{i},rawdata.g_names,ind);
    end
    
    makePlot(switchRank,pulseRank,coeffvarRank,idx,'Coefficient of variation rank');
    saveas(gcf,fullfile(outdir,'rank_vs_cov.pdf'));
    makePlot(switchRank,pulseRank,deRank,idx,'Differential expression rank');
    saveas(gcf,fullfile(outdir,'rank_vs_de.pdf'));
    makePlot(switchRank,pulseRank,drRank,idx,'Dynamic range rank');
    saveas(gcf,fullfile(outdir,'rank_vs_dr.pdf'));
end

% $$$ % Counts above line.
% $$$ covAbove=sum(coeffvarRank<=cutoffLine);
% $$$ deAbove=sum(deRank<=cutoffLine);
% $$$ drAbove=sum(drRank<=cutoffLine);
% $$$ fprintf('Counts above %d out of %d:\n',cutoffLine,n);
% $$$ fprintf('COV: %d DE: %d DR: %d\n',covAbove,deAbove,drAbove);
% $$$ % Switches.
% $$$ for j=1:k
% $$$     if switches(j)
% $$$         stammAbove=sum(switchRank(:,j)<=cutoffLine);
% $$$         fprintf('S%d switch: %d\n',j,stammAbove);
% $$$     end
% $$$ end
% $$$ % Pulses.
% $$$ for j=1:k
% $$$     if pulses(j)
% $$$         stammAbove=sum(pulseRank(:,j)<=cutoffLine);
% $$$         fprintf('S%d pulse: %d\n',j,stammAbove);
% $$$     end
% $$$ end

% Package results.
spScores.switchScore=switchScore;
spScores.pulseScore=pulseScore;
spScores.cov=coeffvar;
spScores.de=de;
spScores.dr=dr;
spScores.switchRank=switchRank;
spScores.pulseRank=pulseRank;
spScores.covRank=coeffvarRank;
spScores.deRank=deRank;
spScores.drRank=drRank;
save(fullfile(outdir,'switchandpulse.mat'),'spScores');

function makePlot(sr,pr,vsrank,idx,ylab)
    figure();

    % Color sequence.
    colors=['krbgcmy'];
    hold on;
    idxCounter=1;
    ncols=length(idx);
    % Plot switches.
    for j=1:k
        if switches(j)
            idxCol=idx{idxCounter};
            plot(sr(idxCol,j),vsrank(idxCol),[colors(j) 'o']);
            idxCounter=idxCounter+1;
        end
    end
    % Plot pulses.
    for j=1:k
        if pulses(j)
            idxCol=idx{idxCounter};
            plot(pr(idxCol,j),vsrank(idxCol),[colors(j) '+']);
            idxCounter=idxCounter+1;
        end
    end
    % Plot cut off line.
    plot(l(:,1), l(:,2),'k:');
    hold off;
    
    xlim([1 length(vsrank)]);
    ylim([1 length(vsrank)]);
    xlabel('Stamm rank');
    ylabel(ylab);
    names={};
    for j=1:k
        if switches(j)
            names=[names num2str(j,'S%d switch')];
        end
    end
    for j=1:k
        if pulses(j)
            names=[names num2str(j,'S%d pulse')];
        end
    end
    legend([names 'Cutoff']);

end

end
