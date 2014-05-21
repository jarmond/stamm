function stammPlotBetaBars(data,beta,ind,outname,relative)
% STAMMPLOTBETABARS Plot beta from fitting as bars
%
%   STAMMPLOTBETABARS(DATA,BETA,IND) Plots rows from BETA indexed by IND as
%   bars. NB: Reverses state naming order so 1st state is labeled 1 (e.g. 4->1,
%   3->2, 2->3, 1->4).
%
%   STAMMPLOTBETABARS(DATA,BETA,IND,OUTNAME) Saves PDF to OUTNAME.
%
%   STAMMPLOTBETABARS(DATA,BETA,IND,OUTNAME,RELATIVE) Plots bars relative to
%   mean if RELATIVE is 1.

if nargin<4
  outname=[];
end
if nargin<5
    relative=0;
end

n=length(ind);
k=size(beta,2);
fig_n=ceil(sqrt(n));
fig_m=ceil(n/fig_n);


figure(1);
clf;
for i=1:n
    b=beta(i,:);
    if relative(i)
        b=b/mean(beta);
    end
    subplot(fig_m,fig_n,i);
    bar(k:-1:1,beta);
    title(data.g_names(ind(i)));
end

if ~isempty(outname)
    saveas(gcf,outname);
end
