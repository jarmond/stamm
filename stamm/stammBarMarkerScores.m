function stammBarMarkerScores(scores,good,bad,outfile)
% STAMMBARMARKERSCORES Plot a bar chart of a selection of marker pair scores.
%
%    STAMMBARMARKERSCORES(SCORES,GOOD,BAD) Plots bar chart of marker pair SCORES. GOOD pairs are
%    shown in blue and and BAD pairs are shown in red.
%
%    STAMMBARMARKERSCORES(SCORES,GOOD,BAD,OUTFILE) Saves PDF to OUTFILE.

if nargin<4
    outfile=[];
end

ng=length(good);
nb=length(bad);

figure(1);
clf;
if ng>0
    barh(1:ng,scores(good),'b');
end
hold on;
if nb>0
    barh(ng+1:ng+nb,scores(bad),'r');
end
hold off;
set(gca,'YDir','reverse');

if ~isempty(outfile)
    saveas(gcf,outfile);
end
