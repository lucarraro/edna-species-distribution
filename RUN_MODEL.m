clear all; close all; clc; rng('shuffle');

%% Load data
load data_wigger.mat; % load Wigger network data
load eDNA_data.mat; % load measured eDNA concentrations 

%% Initialization
ChooseSpecies='Fs'; % if equal to Fs, the model runs for F. sultana; if equal to Tb, the model runs for T. bryosalmonae
if strcmp(ChooseSpecies,'Fs')
    dataset=Fs;
elseif strcmp(ChooseSpecies,'Tb')
    dataset=Tb;
end

% coordinates and reach index for measuring sites
station_coord=[634537.17	240447.56   1;   
               637020.21	238761.43   6;    
               638470.74	236680.43   13;   
               639156.59	233638.54   26;   
               639020.20	233457.90   25;   
               640976.71	229859.55   41;   
               640591.56	226354.76   64;   
               641367.48	225948.74   65;   
               642750.58	219109.92   120;  
               642422.08	218351.57   121;  
               639970.33	217267.33   123;  
               642745.99	222923.29   85;  
               643725.46	222763.63   86;   
               642080.78	224764.12   72;   
               642153.28	225231.59   71];  
N_sites=length(station_coord);
reach_ind=station_coord(:,3);

SitesCalib=1:N_sites; % Here you can exclude some sites from calibration

geometry=v2struct(X,Y,XX,YY,Xc,Yc,Xb,Yb,subcatch,N_reach,AD_pixel,nnodes,outlet,area_upstream,reach); % structure containing river network data used to produce maps
Q0=median(Q_ZOF); % median water discharge from measured time series in Zofingen (site #3) by Swiss Federal Office for the Environment

%% Covariates
Cov_labels=['LocElv'; % 1
            'UpCAr '; % 5
            'GeoMrn'; % 7
            'GeoPea'; % 8
            'GeoWat'; % 9
             ];
N_cov=size(Cov_labels,1);
[covariates,all_covariates]=cov_matrix(Cov_labels,N_cov,reach_ind,mean_subcatchment_altitude,area_upstream,reach_geology_class); % creates covariate matrix
R0=corrcoef(covariates); % matrix of correlation coefficients of the covariates
VIF=diag(inv(R0))'; % calculate Variance Inflation Factors

other=v2struct(all_covariates,reach_upstream,length_downstream,area_local,area_upstream,celerity,length_reach,reach_width,Q0); % structure containing river network data used for model evaluation 

%% Initialize Markov Chain
N_run=2e5; N_steps=(N_cov+4)*N_run/100;

% preallocation of vectors and matrices
p0_vec=zeros(N_steps,1); beta_mat=zeros(N_steps,N_cov); p_mat=zeros(N_steps,N_reach);
tau_vec=zeros(N_steps,1); lSigma_vec=zeros(N_steps,1); C_mat=zeros(N_steps,N_sites); C_all_mat=zeros(N_steps,N_reach);
LogPosterior_vec=zeros(N_steps,1); case_vec=zeros(N_steps,1); Cstar_vec=zeros(N_steps,1); ProbNonDet_mat=zeros(N_steps,N_sites);

% standard deviations for jump distributions
sigma_p0=1e-18; sigma_beta=ones(1,N_cov); sigma_tau=4*3600; sigma_lSigma=0.25; sigma_Cstar=0.2;

% initial parameter values 
% (note that an arbitrary initial parameter set would theoretically work anyways, but it
% might yield LogPosterior = -Inf, making it hard for the chain to start)
p0=1e-18; beta=[0 1 1.5 0 0]; tau=15000; lSigma=0.75; Cstar=1e-16;
for i=1:N_reach
    p(1,i)=p0*exp(all_covariates(i,:)*beta');
end
p_mat(1,:)=p; p0_vec(1)=p0; beta_mat(1,:)=beta; tau_vec(1)=tau; lSigma_vec(1)=lSigma; Cstar_vec(1)=Cstar;
C=eval_concups(p',tau,reach_ind,other); C_mat(1,:)=C';
[LogPosterior,ProbNonDet] = eval_posterior(dataset,Cstar,C,lSigma,SitesCalib);
LogPosterior_vec(1)=LogPosterior; ProbNonDet_mat(1,:)=ProbNonDet;

%% Metropolis-within-Gibbs algorithm
k=1; tic
for ind_run=1:N_run
    for ind_case=1:N_cov+4
        switch ind_case
            case num2cell(1:N_cov)
                % update beta
                betaNew=beta;
                betaNew(1,ind_case)=TruncNormRnd(beta(ind_case),sigma_beta(ind_case),-10,10,1);
                for i=1:N_reach
                    pNew(1,i)=p0*exp(all_covariates(i,:)*betaNew');
                end
                CNew=eval_concups(pNew',tau,reach_ind,other);
                [LogPosteriorNew,ProbNonDetNew] = eval_posterior(dataset,Cstar,CNew,lSigma,SitesCalib);
                tauNew=tau; lSigmaNew=lSigma; p0New=p0; CstarNew=Cstar;
            case N_cov+1
                % update p0
                p0New=TruncNormRnd(p0,sigma_p0,0,1e-15,1);
                for i=1:N_reach
                    pNew(1,i)=p0New*exp(all_covariates(i,:)*beta');
                end
                CNew=eval_concups(pNew',tau,reach_ind,other);
                [LogPosteriorNew,ProbNonDetNew] = eval_posterior(dataset,Cstar,CNew,lSigma,SitesCalib);
                tauNew=tau; lSigmaNew=lSigma; betaNew=beta; CstarNew=Cstar;
            case N_cov+2
                % update tau
                tauNew=TruncNormRnd(tau,sigma_tau,0,40*3600,1);
                CNew=eval_concups(p',tauNew,reach_ind,other);
                [LogPosteriorNew,ProbNonDetNew] = eval_posterior(dataset,Cstar,CNew,lSigma,SitesCalib);
                p0New=p0; betaNew=beta; pNew=p; lSigmaNew=lSigma; CstarNew=Cstar;
            case N_cov+3
                % update lSigma
                lSigmaNew=TruncNormRnd(lSigma,sigma_lSigma,0,2,1);
                [LogPosteriorNew,ProbNonDetNew] = eval_posterior(dataset,Cstar,C,lSigmaNew,SitesCalib);
                p0New=p0; betaNew=beta; pNew=p; tauNew=tau; CstarNew=Cstar;
            case N_cov+4
                % update Cstar
                CstarNew=exp(TruncNormRnd(log(Cstar),sigma_Cstar,-60,0,1));
                [LogPosteriorNew,ProbNonDetNew] = eval_posterior(dataset,CstarNew,C,lSigma,SitesCalib);
                p0New=p0; betaNew=beta; pNew=p; tauNew=tau; lSigmaNew=lSigma; CNew=C;
        end
        if LogPosteriorNew > LogPosterior || rand < exp(LogPosteriorNew-LogPosterior) 
            k=k+1; case_vec(k)=ind_case;
            beta_mat(k,:)=betaNew; p0_vec(k)=p0New; p_mat(k,:)=pNew; tau_vec(k)=tauNew; lSigma_vec(k)=lSigmaNew; Cstar_vec(k)=CstarNew;
            beta=betaNew; p0=p0New; p=pNew; tau=tauNew; lSigma=lSigmaNew; Cstar=CstarNew;
            LogPosterior_vec(k)=LogPosteriorNew; LogPosterior=LogPosteriorNew;
            ProbNonDet_mat(k,:)=ProbNonDetNew; ProbNonDet=ProbNonDetNew;
            C_mat(k,:)=CNew'; C=CNew; 
            C_all_mat(k,:)=eval_concups(pNew',tauNew,1:N_reach,other);
            disp(sprintf('run: %d  -  case: %d  -  time: %.1f  -  accepted: %d / %d  -  LogPosterior: %.1f  -  tau: %.1f',....
                ind_run,ind_case,toc,k,(N_reach+N_cov+3)*(ind_run-1)+ind_case,LogPosterior,tauNew/3600))
        end
    end
    % save results after every 10000 runs
    if mod(ind_run,10000)==0
        if k>10000
        quant_p=quantile(p_mat(10001:k,:),[0.025 0.25 0.5 0.75 0.975]);
        quant_C=quantile(C_all_mat(10001:k,:),[0.025 0.25 0.5 0.75 0.975]);
        else
        quant_p=[]; quant_C=[];
        end
        save(['results_',ChooseSpecies,'.mat'],...
            'LogPosterior_vec','p0_vec','beta_mat','tau_vec','lSigma_vec','case_vec','C_mat','Cstar_vec','ProbNonDet_mat','quant_p','quant_C')
        disp(sprintf(' ')); disp(sprintf(' '));
        disp(sprintf('run: %d',ind_run))
        disp(sprintf(' ')); disp(sprintf(' '));
    end
end
disp(sprintf('End of process'))


%% cut burn-in and end
thr_burnin=10000;

beta_mat(k+1:end,:)=[]; tau_vec(k+1:end)=[]; LogPosterior_vec(k+1:end)=[]; case_vec(k+1:end)=[]; ProbNonDet_mat(k+1:end,:)=[];
p_mat(k+1:end,:)=[]; p0_vec(k+1:end)=[]; lSigma_vec(k+1:end)=[]; C_mat(k+1:end,:)=[]; Cstar_vec(k+1:end)=[];

beta_mat(1:thr_burnin,:)=[]; tau_vec(1:thr_burnin)=[]; LogPosterior_vec(1:thr_burnin)=[]; case_vec(1:thr_burnin)=[]; ProbNonDet_mat(1:thr_burnin,:)=[];
p0_vec(1:thr_burnin)=[]; lSigma_vec(1:thr_burnin)=[]; C_mat(1:thr_burnin,:)=[]; Cstar_vec(1:thr_burnin)=[];

quant_p=quantile(p_mat,[0.025 0.25 0.5 0.75 0.975]);
quant_C=quantile(C_all_mat,[0.025 0.25 0.5 0.75 0.975]);

save(['results_',ChooseSpecies,'.mat'],'LogPosterior_vec','p0_vec','beta_mat','tau_vec','lSigma_vec','case_vec','C_mat','Cstar_vec','ProbNonDet_mat','quant_p','quant_C')

