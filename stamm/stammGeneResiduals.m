function rss=stammGeneResiduals(data,model,resultfile,sel)
% STAMMGENERSS Calculate residuals of indexed gene

load(resultfile);
[S,P]=eval([model '(result.W,result.beta,data.t)']);
S=log2(S);

n=length(sel);
rss=zeros(n,length(data.t));
for i=1:n
    rss(i,:)=data.g(result.ind(sel(i)),:)-S(sel(i));
end

