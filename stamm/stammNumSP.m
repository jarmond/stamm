function [k,nw]=stammNumSP(modelId)
% STAMMNUMSP Number of parameters and states for given model
%
%    [K,NW] = STAMMNUMSP(MODELID) Returns number of states K and number of
%    transition parameters NW for model identified by MODELID.

switch modelId
  case 1
    k=2;
    nw=2;
  case {2,3}
    k=3;
    nw=4;
  case 4
    k=4;
    nw=4;
  case {5,6,7}
    k=4;
    nw=6;
  case 8
    k=5;
    nw=8;
  case 9
    k=4;
    nw=8;
  case 10
    k=2;
    nw=1;
  case {11,15,16,17}
    k=4;
    nw=3;
  case 12
    k=5;
    nw=4;
  case {13,14}
    k=3;
    nw=2;
  otherwise
    error('Unknown model');
end
