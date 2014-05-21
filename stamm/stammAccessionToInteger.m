function accint=stammAccessionToInteger(acc)
% STAMMACCESSIONTOINTEGER Convert accession number to integer

if ~iscell(acc)
    acc={acc};
end

accint=zeros(length(acc),1);
for i=1:length(acc)
    % match the numbers
    r=regexp(acc{i},'\d+','match');
    if isempty(r)
        accint(i)=0;
    else
        accint(i)=str2double(r(1));
    end
end
