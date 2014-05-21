function stammPlotOptimalPerturbations(beta,m,filepref,path,m1,p)
% STAMMPLOTOPTIMALPERTURBATIONS Scatter standard betas against robustness checking betas
%
%   STAMMPLOTOPTIMALPERTURBATIONS(BETA,M,FILEPREF,PATH,M1,P) BETA should
%   correspond directly to those with prefix FILEPREF in the files in PATH. The
%    loaded file is the concatenation of FILEPREF and a number from 1 to M with
%    two digits. PDF of plot is saved in PATH. P is the number perturbation
%    factors to check.

if nargin<5
    m1=1;
end

n=m-m1+1;

pc=zeros(length(p),4);
for j=1:length(p)
    c=zeros(n,1);
    for i=1:n
        % Read each perturbation file.
        tfile=[path '/' num2str(j,'f%02d') '/' filepref num2str(i+m1-1,'%02d')];
        %tfile=[path '/' filepref num2str(i+m1-1,'%02d') '-a02.mat'];
        fprintf('Reading: %s\n',tfile);
        load(tfile);
        X=corr(beta,result.beta);
        c(i)=mean(diag(X));
    end
    % Find average for this perturbation factor.
    pc(j,:)=[mean(c) std(c) min(c) max(c)];
end

% Plot mean and std dev.
figure(1);
clf;
hold on
errorbar(p,pc(:,1),pc(:,2));
plot(p,pc(:,3),'r',p,pc(:,4),'r');
ylim([0 1]);
hold off
xlabel('Std dev of noise');
ylabel('Mean correlation');

saveas(gcf,[path '/' filepref '_optimalperturbation.pdf']);
