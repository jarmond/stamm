function ig=invgamrnd(a,b,m,n)
% INVGAMRND Generate random numbers from inverse gamma distribution

ig=1.0./gamrnd(a,1.0./b,m,n);
