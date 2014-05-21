function stammPrintBeta(data,filename,head,latex,sorting,rescale,outfile)
% STAMMPRINTBETA Print out fitted beta table
%
%   STAMMPRINTBETA(DATA,FILENAME) Prints out the betas in FILENAME
%   corresponding to genes in DATA.
%
%   STAMMPRINTBETA(DATA,FILENAME,HEAD,LATEX,SORTING,RESCALE) Optional
%   parameters may be supplied to modify output. HEAD prints a header line.
%   LATEX formats as a LaTeX table. SORTING specifies while column to sort
%   descending by. RESCALE divides each column by the mean of the column.
%
%   STAMMPRINTBETA(DATA,FILENAME,HEAD,LATEX,SORTING,RESCALE,OUTFILE) Saves
%   output into text file OUTFILE.
    
if nargin<3
    head=true;
end
if nargin<4
    latex=false;
end
if nargin<5
    sorting=0;
end
if nargin<6
    rescale=false;
end
if nargin<7
    outfile=[];
end

load(filename);

if rescale
    % Rescale by column means.
    mbeta=mean(beta);
    for i=1:size(beta,2)
        beta(:,i)=beta(:,i)/mbeta(i);
    end
end

% Print beta.
printout(1);

if ~isempty(outfile)
    % Open output file.
    fid=fopen(outfile,'w');
    printout(fid);
    fclose(fid);
end

%% NESTED FUNCTIONS
   function printout(target)

   if head
       % Print header.
       fprintf(target,'%15s ','Gene');
       fprintf(target,'%7d ',1:size(beta,2));
       fprintf(target,'\n');
   end

   if sorting<1
       % No sorting.
       seq=1:size(beta,1);
   else
       % Sort by column sorting.
       [y,seq]=sort(beta(:,sorting),1,'descend');
   end
    
   % Print each row.
   for j=1:length(seq)
       i=seq(j);
       g_str=data.g_names{ind(i)};
       fprintf(target,'%15s ',g_str);

       if latex
           fprintf(target,'& %7.2f ',beta(i,:));
       else
           fprintf(target,'%7.2f ',beta(i,:));
       end
       
       if latex
           fprintf(target,'\\\\\n');
       else
           fprintf(target,'\n');
       end
   end


   end

end
