function [W,beta]=stammUnpack(modelId,x,n)
% STAMMUNPACK Unpack parameters packed by stammPackParams into W and beta

[k,nw]=stammNumSP(modelId);

switch modelId
  case 1
    W=x(1:nw);
  case 2
  W=zeros(3);
  W(1,2)=x(1);
  W(2,1)=x(2);
  W(2,3)=x(3);
  W(3,2)=x(4);
  case 3
  W=zeros(3);
  W(1,3)=x(1);
  W(2,3)=x(2);
  W(3,1)=x(3);
  W(3,2)=x(4);
  case 4
  W=zeros(4);
  W(1,2)=x(1);
  W(2,1)=x(2);
  W(2,4)=x(3);
  W(3,4)=x(4);
  case 5
  W=zeros(4);
  W(1,2)=x(1);
  W(2,1)=x(2);
  W(2,4)=x(3);
  W(3,4)=x(4);
  W(4,2)=x(5);
  W(4,3)=x(6);
  case 6
  W=zeros(4);
  W(1,2)=x(1);
  W(2,1)=x(2);
  W(2,3)=x(3);
  W(3,2)=x(4);
  W(3,4)=x(5);
  W(4,3)=x(6);
  case 7
  W=zeros(4);
    W(1,2)=x(1);
    W(2,1)=x(2);
    W(2,3)=x(3);
    W(2,4)=x(4);
    W(3,2)=x(5);
    W(4,2)=x(6);
  case 8
  W=zeros(4);
    W(1,2)=x(1);
    W(2,1)=x(2);
    W(2,3)=x(3);
    W(3,2)=x(4);
    W(3,4)=x(5);
    W(4,3)=x(6);
    W(4,5)=x(7);
    W(5,4)=x(8);
  case 9
  W=zeros(4);
    W(1,2)=x(1);
    W(1,3)=x(2);
    W(2,1)=x(3);
    W(2,3)=x(4);
    W(3,1)=x(5);
    W(3,2)=x(6);
    W(3,4)=x(7);
    W(4,3)=x(8);
  case 10
    W=x(1);
  case {11,17}
  W=zeros(4);
    W(1,2)=x(1);
    W(2,3)=x(2);
    W(3,4)=x(3);
  case 12
    W=zeros(5);
    W(1,2)=x(1);
    W(2,3)=x(2);
    W(3,4)=x(3);
    W(4,5)=x(4);
  case 13
    W=zeros(3);
    W(1,2)=x(1);
    W(2,3)=x(2);
  case 14
    W=zeros(3);
    W(1,3)=x(1);
    W(2,3)=x(2);
  case 15
    W=zeros(4);
    W(1,2)=x(1);
    W(2,4)=x(2);
    W(3,2)=x(3);
  case 16
    W=zeros(4);
    W(1,2)=x(1);
    W(2,4)=x(2);
    W(3,4)=x(3);
  otherwise
  error('Unknown model');
end

if nargout>1
    beta=reshape(x(nw+1:end),n,k);
end
