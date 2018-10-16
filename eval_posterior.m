function [LogPosterior,phi,Loglik,LogPrior] = eval_posterior(dataset,Cstar,C,lSigma,SitesCalib)
% Evaluate posterior distribution (proportional to the likelihood) for a given parameter set
for i=1:15
    phi(i)=exp(-C(i)/Cstar); % probability of non-detection
    tmp=dataset.(['S',num2str(i)]); % load data
    CountZeros=sum(tmp==0); tmp(tmp==0)=[]; tmp(isnan(tmp))=[];
    NonZeros=tmp; CountNonZeros=length(tmp); % count zeros and non-zeros
    if ismember(i,SitesCalib) 
        ContrLoglik(i,1)=CountZeros*log(phi(i)) + CountNonZeros*log(1-phi(i)) + sum(log(normpdf(log(NonZeros),log(C(i)),lSigma)));
    else
        ContrLoglik(i,1)=0;
    end
end
Loglik=sum(ContrLoglik);

LogPosterior = Loglik;

end

