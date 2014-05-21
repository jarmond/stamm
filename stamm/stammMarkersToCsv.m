function stammMarkersToCsv(data,markers,beta,n,outfile)
% STAMMMARKERSTOCSV Dumps beta to CSV file
%
%   STAMMMARKERSTOCSV(DATA,MARKERS,BETA,N,OUTFILE) Outputs top N marker pairs
%   in MARKERS, together with BETA for each gene, with names from DATA. Saves
%   as CSV in OUTFILE.

% Dlmwrite doesn't support cell arrays so to output strings for gene names we
% use fprintf to manually write csv.
f=fopen(outfile,'w+t');

[n k]=size(beta);
% Header.
fprintf(f,'Rank,Gene 1,Gene 2,Score,Gene 1 - ');
for i=1:k
    fprintf(f,['b_' num2str(i) ',']);
end
fprintf(f,'Gene 2 - ');
for i=1:k
    if i==k
        fprintf(f,['b_' num2str(i) '\n']);
    else
        fprintf(f,['b_' num2str(i) ',']);
    end
end

for i=1:n
    fprintf(f,'%d,%s,%s,%.3f,',i,data.g_names{markers.pairsind(i,1)}, ...
            data.g_names{markers.pairsind(i,2)},markers.disc(i));
    b = beta(markers.pairs(i,1),:);
    fprintf(f,'%.2f,',b);
    b = beta(markers.pairs(i,2),:);
    fprintf(f,'%.2f,',b(1:end-1));
    fprintf(f,'%.2f\n',b(end));
end

fclose(f);
