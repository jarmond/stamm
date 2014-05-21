function stammStatesMarkers(disc,pairs,g,s,average,outfile)
% STAMMSTATESMARKERS Plot number of states against marker rank
%
%    STAMMSTATESMARKERS(DISC,PAIRS,G) Plots the number of states in single-cell
%    data G with moving averaging window 50 against marker rank, together with a
%    randomly permuted ranking from PAIRS and discriminatory score DISC.
%
%    STAMMSTATESMARKERS(DISC,PAIRS,G,S) Use moving average window S.
%
%    STAMMSTATESMARKERS(DISC,PAIRS,G,S,AVERAGE) Repeat random permutations
%    AVERAGE times.
%
%    STAMMSTATESMARKERS(DISC,PAIRS,G,S,AVERAGE,OUTFILE) Save plot as PDF to
%    OUTFILE.

if nargin<6
    outfile=[];
end
if nargin<5
    average=0;
end
if nargin<4
    s=50;
end

n=size(pairs,1);
[counts,states]=stammCountMarkers(pairs,g);

% random permutation
if average>0
    randstates=zeros(n,average);
    for i=1:average
        randstates(:,i)=tsmovavg(states(randperm(n))','s',s)';
    end
    stddevRandstatesmovavg=std(randstates,0,2);
    randstatesmovavg=mean(randstates,2);
else
    randstates=states(randperm(length(states)));
    randstatesmovavg=tsmovavg(randstates','s',s);
end

statesmovavg=tsmovavg(states','s',s);
discmovavg=tsmovavg(disc','s',s);

figure(1);
clf;
xaxis=1:(n-s+1);
[ax,h1,h2]=plotyy(xaxis,statesmovavg(s:end),xaxis,discmovavg(s:end));
hold on;
plot(xaxis,randstatesmovavg(s:end),'k');
if average>0
    plot(xaxis,randstatesmovavg(s:end)+stddevRandstatesmovavg(s:end),'r',...
         xaxis,randstatesmovavg(s:end)-stddevRandstatesmovavg(s:end),'r');
end
hold off;
xlim(ax(1),[1 n-s]);
xlim(ax(2),[1 n-s]);
xlabel('Rank');
ylabel(ax(1),'avg # clusters (window size 50)');
ylabel(ax(2),'Distance score');

if ~isempty(outfile)
    saveas(gcf,outfile);
end
