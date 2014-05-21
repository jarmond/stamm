function data=stammReadBingoData(filename,skiplines)
% STAMMREADBINGODATA Reads BINGO output file

fid=fopen(filename);
fieldfmt='%d %f %f %d %d %d %d %q %q';
data=textscan(fid,fieldfmt,'Delimiter','\t','HeaderLines',skiplines,'BufSize',16384);
fclose(fid);
