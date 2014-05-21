function stammRandomSeed(seed)
% STAMMRANDOMSEED Set default stream to clock based seed, or specified integer.

if nargin<1
    seed=sum(100*clock);
end
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
