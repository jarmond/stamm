function x=stammPackParams(g_s,w,mid)
% STAMMPACKPARAMS Pack parameters into one dimensional array

% Extract w's from w matrix as per model id.
switch mid
  case 1
    wx=w;
  case 2
    wx=zeros(4,1);
    wx(1)=w(1,2);
    wx(2)=w(2,1);
    wx(3)=w(2,3);
    wx(4)=w(3,2);
  case 3
    wx=zeros(4,1);
    wx(1)=w(1,3);
    wx(2)=w(2,3);
    wx(3)=w(3,1);
    wx(4)=w(3,2);
  case 4
    wx=zeros(4,1);
    wx(1)=w(1,2);
    wx(2)=w(2,1);
    wx(3)=w(2,4);
    wx(4)=w(3,4);
  case 5
    wx=zeros(6,1);
    wx(1)=w(1,2);
    wx(2)=w(2,1);
    wx(3)=w(2,4);
    wx(4)=w(3,4);
    wx(5)=w(4,2);
    wx(6)=w(4,3);
  case 6
    wx=zeros(6,1);
    wx(1)=w(1,2);
    wx(2)=w(2,1);
    wx(3)=w(2,3);
    wx(4)=w(3,2);
    wx(5)=w(3,4);
    wx(6)=w(4,3);
  case 7
    wx=zeros(6,1);
    wx(1)=w(1,2);
    wx(2)=w(2,1);
    wx(3)=w(2,3);
    wx(4)=w(2,4);
    wx(5)=w(3,2);
    wx(6)=w(4,2);
  case 8
    wx=zeros(8,1);
    wx(1)=w(1,2);
    wx(2)=w(2,1);
    wx(3)=w(2,3);
    wx(4)=w(3,2);
    wx(5)=w(3,4);
    wx(6)=w(4,3);
    wx(7)=w(4,5);
    wx(8)=w(5,4);
  case 9
    wx=zeros(8,1);
    wx(1)=w(1,2);
    wx(2)=w(1,3);
    wx(3)=w(2,1);
    wx(4)=w(2,3);
    wx(5)=w(3,1);
    wx(6)=w(3,2);
    wx(7)=w(3,4);
    wx(8)=w(4,3);
  case 10
    wx=w;
  case {11,17}
    wx=zeros(3,1);
    wx(1)=w(1,2);
    wx(2)=w(2,3);
    wx(3)=w(3,4);
  case 12
    wx=zeros(4,1);
    wx(1)=w(1,2);
    wx(2)=w(2,3);
    wx(3)=w(3,4);
    wx(4)=w(4,5);
  case 13
    wx=zeros(2,1);
    wx(1)=w(1,2);
    wx(2)=w(2,3);
  case 14
    wx=zeros(2,1);
    wx(1)=w(1,3);
    wx(2)=w(2,3);
  case 15
    wx=zeros(3,1);
    wx(1)=w(1,2);
    wx(2)=w(2,4);
    wx(3)=w(3,2);
  case 16
    wx=zeros(3,1);
    wx(1)=w(1,2);
    wx(2)=w(2,4);
    wx(3)=w(3,4);
end
    
ngs=prod(size(g_s));
% Reshape into 1-d array.
x=[wx; reshape(g_s,[ngs 1])];
