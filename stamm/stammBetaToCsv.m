function stammBetaToCsv(data,ind,beta,outfile)
% STAMMBETATOCSV Dumps beta to CSV file
%
%   STAMMBETATOCSV(DATA,IND,BETA,OUTFILE) Outputs rows of BETA indexed by IND to CSV
%   file named OUTFILE.

% Dlmwrite doesn't support cell arrays so to output strings for gene names we
% use fprintf to manually write csv.
f=fopen(outfile,'w+t');

[n k]=size(beta);
% Header.
fprintf(f,'Gene,');
for i=1:k
    if i==k
        fprintf(f,['b_' num2str(i) '\n']);
    else
        fprintf(f,['b_' num2str(i) ',']);
    end
end

for i=1:n
    fprintf(f,'%s,',data.g_names{ind(i)});
    b = beta(i,:);
    fprintf(f,'%.4f,',b(1:end-1));
    fprintf(f,'%.4f\n',b(end));
end

fclose(f);
