clear all; close all; clc

%% Load data
load data_wigger.mat; % Wigger network data
load eDNA_data.mat; % measured eDNA concentrations 

%% Initialization

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
other=v2struct(all_covariates,reach_upstream,length_downstream,area_local,area_upstream,celerity,length_reach,reach_width,Q0); % structure containing river network data used for model evaluation 

%% read data
load results_Fs.mat
Fs_beta_mat=beta_mat; Fs_Cstar_vec=Cstar_vec; Fs_LogPosterior_vec=LogPosterior_vec; Fs_tau_vec=tau_vec; Fs_lSigma_vec=lSigma_vec;
Fs_p0_vec=p0_vec; Fs_C_mat=C_mat; Fs_ProbNonDet_mat=ProbNonDet_mat; Fs_quant_p=quant_p; Fs_quant_C=quant_C;

load results_Tb.mat
Tb_beta_mat=beta_mat; Tb_Cstar_vec=Cstar_vec; Tb_LogPosterior_vec=LogPosterior_vec; Tb_tau_vec=tau_vec; Tb_lSigma_vec=lSigma_vec;
Tb_p0_vec=p0_vec; Tb_C_mat=C_mat;  Tb_ProbNonDet_mat=ProbNonDet_mat; Tb_quant_p=quant_p; Tb_quant_C=quant_C;

clear beta_mat p0_vec Cstar_vec LogPosterior_vec lSigma_vec ProbNonDet_mat case_vec  tau_vec C_mat 

%% posterior distributions
dim=60; % number of bins in the x-axis for posterior distribution plots

% beta values
figure('units','centimeters','position',[0 0 18 13])
subplot(2,3,1); histogram(Fs_beta_mat(:,1),[-5:10/dim:5],...
    'Normalization','probability','FaceColor',[100/255 34/255 102/255],'EdgeColor','n');
box off; set(gca,'tickdir','out','xlim',[-5 5],'ylim',[0 0.3],'ytick',[0:0.1:0.3])
hold on; histogram(Tb_beta_mat(:,1),[-5:10/dim:5],...
    'Normalization','probability','FaceColor',[101/255 155/255 65/255],'EdgeColor','n');
xlabel(Cov_labels(1,:)); ylabel('Frequency [-]')
subplot(2,3,2); histogram(Fs_beta_mat(:,3),[-5:10/dim:5],...
    'Normalization','probability','FaceColor',[100/255 34/255 102/255],'EdgeColor','n');
box off; set(gca,'tickdir','out','xlim',[-5 5],'ylim',[0 0.3],'ytick',[0:0.1:0.3])
hold on; histogram(Tb_beta_mat(:,3),[-5:10/dim:5],...
    'Normalization','probability','FaceColor',[101/255 155/255 65/255],'EdgeColor','n');
xlabel(Cov_labels(3,:))
subplot(2,3,3); histogram(Fs_beta_mat(:,5),[-5:10/dim:5],...
    'Normalization','probability','FaceColor',[100/255 34/255 102/255],'EdgeColor','n');
box off; set(gca,'tickdir','out','xlim',[-5 5],'ylim',[0 0.3],'ytick',[0:0.1:0.3])
hold on; histogram(Tb_beta_mat(:,5),[-5:10/dim:5],...
    'Normalization','probability','FaceColor',[101/255 155/255 65/255],'EdgeColor','n');
xlabel(Cov_labels(5,:))
subplot(2,3,4); histogram(Fs_beta_mat(:,2),[-5:10/dim:5],...
    'Normalization','probability','FaceColor',[100/255 34/255 102/255],'EdgeColor','n');
box off; set(gca,'tickdir','out','xlim',[-5 5],'ylim',[0 0.3],'ytick',[0:0.1:0.3])
hold on; histogram(Tb_beta_mat(:,2),[-5:10/dim:5],...
    'Normalization','probability','FaceColor',[101/255 155/255 65/255],'EdgeColor','n');
xlabel(Cov_labels(2,:)); ylabel('Frequency [-]')
subplot(2,3,5); histogram(Fs_beta_mat(:,4),[-5:10/dim:5],...
    'Normalization','probability','FaceColor',[100/255 34/255 102/255],'EdgeColor','n');
box off; set(gca,'tickdir','out','xlim',[-5 5],'ylim',[0 0.3],'ytick',[0:0.1:0.3])
hold on; histogram(Tb_beta_mat(:,4),[-5:10/dim:5],...
    'Normalization','probability','FaceColor',[101/255 155/255 65/255],'EdgeColor','n');
xlabel(Cov_labels(4,:))

% other parameters
figure('units','centimeters','position',[0 0 18 13])
subplot(2,3,1); histogram(Fs_tau_vec/3600,[0:40/dim:40],...
    'Normalization','Probability','FaceColor',[100/255 34/255 102/255],'EdgeColor','n');
box off; set(gca,'tickdir','out','xlim',[0 40],'ylim',[0 0.3],'ytick',[0:0.1:0.3])
hold on; histogram(Tb_tau_vec/3600,[0:40/dim:40],...
    'Normalization','Probability','FaceColor',[101/255 155/255 65/255],'EdgeColor','n');
title('\tau [h]'); ylabel('Frequency [-]')
subplot(2,3,2); histogram(Fs_lSigma_vec,[0:1.5/dim:1.5],...
    'Normalization','Probability','FaceColor',[100/255 34/255 102/255],'EdgeColor','n');
box off; set(gca,'tickdir','out','xlim',[0 1.5],'ylim',[0 0.3],'ytick',[0:0.1:0.3],'xtick',[0:0.75:1.5]); hold on
histogram(Tb_lSigma_vec,[0:1.5/dim:1.5],...
    'Normalization','Probability','FaceColor',[101/255 155/255 65/255],'EdgeColor','n');
title('\sigma')
subplot(2,3,3); histogram(log10(Fs_Cstar_vec),[-17:2/dim:-15],...
'Normalization','Probability','FaceColor',[100/255 34/255 102/255],'EdgeColor','n');
hold on; histogram(log10(Tb_Cstar_vec),[-17:2/dim:-15],...
    'Normalization','Probability','FaceColor',[101/255 155/255 65/255],'EdgeColor','n');
box off; set(gca,'tickdir','out','xlim',[-17 -15],'ylim',[0 0.3],'ytick',[0:0.1:0.3])
title('C* [mol/L]')
subplot(2,3,4); histogram(Fs_p0_vec,[0:2e-17/dim:2e-17],...
    'Normalization','Probability','FaceColor',[100/255 34/255 102/255],'EdgeColor','n');
box off; set(gca,'tickdir','out','xlim',[0 2e-17],'ylim',[0 0.3],'ytick',[0:0.1:0.3],'xtick',[0 1e-17 2e-17])
hold on; histogram(Tb_p0_vec,[0:2e-17/dim:2e-17],...
    'Normalization','Probability','FaceColor',[101/255 155/255 65/255],'EdgeColor','n');
title('p_0 [mol m^{-2} s^{-1}]'); ylabel('Frequency [-]')

%% plot non-detection probability distributions for all sites
figure('units','centimeters','position',[0 0 20 15])
for i=1:15
    subplot(3,5,i)
    histogram(Fs_ProbNonDet_mat(:,i),[0:0.025:1],'FaceColor',[100/255 34/255 102/255],'EdgeColor','n','Normalization','Probability');
    box off; set(gca,'xlim',[0 1],'tickdir','out','ylim',[0 0.3]); hold on;
    tmp=Fs.(['S',num2str(i)]);
    CountZeros=sum(tmp==0); tmp(tmp==0)=[]; tmp(isnan(tmp))=[];
    NonZeros=tmp; CountNonZeros=length(tmp);
    NonDet(i)=CountZeros/(CountZeros+CountNonZeros);
    plot([NonDet(i) NonDet(i)],[0 0.3],'color',[100/255 34/255 102/255],'LineWidth',1.5)
    histogram(Tb_ProbNonDet_mat(:,i),[0:0.025:1],'FaceColor',[101/255 155/255 65/255],'EdgeColor','n','Normalization','Probability');
    tmp=Tb.(['S',num2str(i)]);
    CountZeros=sum(tmp==0); tmp(tmp==0)=[]; tmp(isnan(tmp))=[];
    NonZeros=tmp; CountNonZeros=length(tmp);
    NonDet(i)=CountZeros/(CountZeros+CountNonZeros);
    plot([NonDet(i) NonDet(i)],[0 0.3],'--','color',[101/255 155/255 65/255],'LineWidth',1.5)
    title(['#',num2str(i)]);
    if mod(i,5)==1; ylabel('Frequency [-]'); end
end

%% plot cumulative distribution functions for modelled eDNA concentration (C) at all sites
figure('units','centimeters','position',[0 0 20 15])
for i=1:15
    % F. sultana
    % build empirical distribution function
    clear edf
    max_lim=nanmax([Fs.(['S',num2str(i)]);  Tb.(['S',num2str(i)])]);
    max_lim=ceil(0.1*max_lim*1e17)/1e16;
    vec=[unique(Fs_C_mat(:,i)); max_lim];
    tmp=Fs.(['S',num2str(i)]); tmp(tmp==0)=[]; tmp(isnan(tmp))=[];
    Fs_NonZeros=sort(tmp); CountNonZeros=length(tmp);
    pos_vec=(1:length(vec))/length(vec);
    pos_data=(1:CountNonZeros)/CountNonZeros;
    vec_edf=zeros(2*CountNonZeros+2,1); pos_edf=zeros(2*CountNonZeros+2,1);
    for j=1:CountNonZeros
        vec_edf(j*2:j*2+1)=Fs_NonZeros(j);
        pos_edf(j*2+1:j*2+2)=j/CountNonZeros;
    end
    vec_edf(end)=max_lim;
    for j=1:length(vec)
        ind=find(Fs_NonZeros<vec(j),1,'last');
        if isempty(ind)
            edf(1,j)=0;
        else
            edf(1,j)=pos_data(ind);
        end
    end
    % KS test
    if max(abs(pos_vec-edf)) < 1.36*sqrt(1/CountNonZeros)
        Fs_alpha=0.05;
    elseif max(abs(pos_vec-edf)) < 1.63*sqrt(1/CountNonZeros)
        Fs_alpha=0.01;
    elseif max(abs(pos_vec-edf)) < 1.73*sqrt(1/CountNonZeros)
        Fs_alpha=0.005;
    elseif max(abs(pos_vec-edf)) < 1.97*sqrt(1/CountNonZeros)
        Fs_alpha=0.001;
    else
        Fs_alpha=0;
    end
    % plot
    subplot(3,5,i);  hold on; plot(vec_edf,pos_edf,'color',[100/255 34/255 102/255]);
    plot(Fs_NonZeros,pos_data,'o','markerfacecolor',[100/255 34/255 102/255],'markeredgecolor','n','markersize',5); box off
    l1=plot(vec,pos_vec,'color',[100/255 34/255 102/255],'linewidth',1.5);
    % T. bryosalmonae
    % build empirical distribution function
    clear edf
    vec=[unique(Tb_C_mat(:,i)); max_lim];
    tmp=Tb.(['S',num2str(i)]); CountZeros=sum(tmp==0); tmp(tmp==0)=[]; tmp(isnan(tmp))=[];
    Tb_NonZeros=sort(tmp); CountNonZeros=length(tmp);
    pos_vec=(1:length(vec))/length(vec);
    pos_data=(1:CountNonZeros)/CountNonZeros;
    vec_edf=zeros(2*CountNonZeros+2,1); pos_edf=zeros(2*CountNonZeros+2,1);
    for j=1:CountNonZeros
        vec_edf(j*2:j*2+1)=Tb_NonZeros(j);
        pos_edf(j*2+1:j*2+2)=j/CountNonZeros;
    end
    vec_edf(end)=max_lim;
    for j=1:length(vec)
        ind=find(Tb_NonZeros<vec(j),1,'last');
        if isempty(ind)
            edf(1,j)=0;
        else
            edf(1,j)=pos_data(ind);
        end
    end
    % KS test
    if max(abs(pos_vec-edf)) < 1.36*sqrt(1/CountNonZeros)
        Tb_alpha=0.05;
    elseif max(abs(pos_vec-edf)) < 1.63*sqrt(1/CountNonZeros)
        Tb_alpha=0.01;
    elseif max(abs(pos_vec-edf)) < 1.73*sqrt(1/CountNonZeros)
        Tb_alpha=0.005;
    elseif max(abs(pos_vec-edf)) < 1.97*sqrt(1/CountNonZeros)
        Tb_alpha=0.001;
    else
        Tb_alpha=0;
    end
    % plot
    plot(vec_edf,pos_edf,'color',[101/255 155/255 65/255]);
    plot(Tb_NonZeros,pos_data,'o','markerfacecolor',[101/255 155/255 65/255],'markeredgecolor','n','markersize',5)
    l2=plot(vec,pos_vec,'color',[101/255 155/255 65/255],'linewidth',1.5);
    set(gca,'tickdir','out','xlim',[0 max_lim],'xtick',max_lim*[0 0.5 1])
    legend([l1 l2],['\alpha: ',num2str(Fs_alpha)],['\alpha: ',num2str(Tb_alpha)],'location','southeast')
end

%% Plot maps
DrawWiggerMap('netw',Fs_quant_p(3,:),2e-18,0,'Fs p [mol m^{-2} s^{-1}]','GER',geometry,1,1,1) % Fs production
DrawWiggerMap('netw',Tb_quant_p(3,:),2e-18,0,'Tb p [mol m^{-2} s^{-1}]','GER',geometry,1,1,1) % Tb production
DrawWiggerMap('netw',Fs_quant_C(3,:),3e-16,0,'Fs C [mol/L]','GER',geometry,1,1,1) % Fs concentration
DrawWiggerMap('netw',Tb_quant_C(3,:),3e-16,0,'Tb C [mol/L]','GER',geometry,1,1,1) % Tb concentration

%% correlation plots of eDNA production vs. covariates
figure('units','centimeters','position',[0 0 25 11])
% F. sultana
for i=1:5
subplot(2,5,i)
plot(all_covariates(:,i),Fs_quant_p(3,:),'.k'); hold on; box off
b=polyfit(all_covariates(:,i),Fs_quant_p(3,:)',1);
l1=plot([-1 1],[-b(1)+b(2) b(1)+b(2)],'r');
K=corrcoef(all_covariates(:,i),Fs_quant_p(3,:)'); K=K(2,1);
legend(l1,num2str(K));
set(gca,'tickdir','out','ylim',[0 1.5e-18])
if i==1; ylabel('Fs production [mol m^{-2} s^{-1}]'); end
end
% T. bryosalmonae
for i=1:5
subplot(2,5,i+5)
plot(all_covariates(:,i),Tb_quant_p(3,:),'.k'); hold on; box off
b=polyfit(all_covariates(:,i),Tb_quant_p(3,:)',1);
l1=plot([-1 1],[-b(1)+b(2) b(1)+b(2)],'r');
K=corrcoef(all_covariates(:,i),Tb_quant_p(3,:)'); K=K(2,1);
legend(l1,num2str(K));
set(gca,'tickdir','out','ylim',[0 6e-18])
xlabel(Cov_labels(i,:))
if i==1; ylabel('Tb production [mol m^{-2} s^{-1}]'); end
end

% correlation plot of Fs vs Tb eDNA productions
figure;
plot(Fs_quant_p(3,:),Tb_quant_p(3,:),'.k'); hold on; box off
b=polyfit(Fs_quant_p(3,:),Tb_quant_p(3,:),1);
l1=plot([0 1.5e-19],[b(2) 1.5e-19*b(1)+b(2)],'r');
K=corrcoef(Fs_quant_p(3,:),Tb_quant_p(3,:)); K=K(2,1);
legend(l1,num2str(K));
set(gca,'tickdir','out','ylim',[0 6e-18],'xlim',[0 1.5e-18],'ytick',1e-18*[0 2 4 6])

%% plot maps of covariate values
for i=1:5
    DrawWiggerMap('netw',all_covariates(:,i),1,-1,Cov_labels(i,:),'GER',geometry,1,1,1);
end

%% scatter plot of observed vs measured eDNA concentrations (for both F. sultana and T. bryosalmonae) 
for i=1:15
    tmp=Fs.(['S',num2str(i)]); tmp(tmp==0)=[]; tmp(isnan(tmp))=[];
    Fs_obs_median(i,1)=median(tmp);
    tmp=Tb.(['S',num2str(i)]); tmp(tmp==0)=[]; tmp(isnan(tmp))=[];
    Tb_obs_median(i,1)=median(tmp);
end

NS_Fs=1-sum((Fs_obs_median'-Fs_quant_C(3,reach_ind)).^2)/sum((Fs_obs_median-mean(Fs_obs_median)).^2);
NS_Tb=1-nansum((Tb_obs_median'-Tb_quant_C(3,reach_ind)).^2)/nansum((Tb_obs_median-nanmean(Tb_obs_median)).^2);

figure; hold on; plot([0 1],[0 1],'color',[.5 .5 .5]); box off; grid on
l1=plot(Fs_obs_median',Fs_quant_C(3,reach_ind),'o','markerfacecolor',[100/255 34/255 102/255],'markeredgecolor','n');
l2=plot(Tb_obs_median',Tb_quant_C(3,reach_ind),'o','markerfacecolor',[101/255 155/255 65/255],'markeredgecolor','n');
set(gca,'xlim',[0 1.2e-16],'ylim',[0 1.2e-16],'tickdir','out','xtick',1e-17*[0:3:12],'ytick',1e-17*[0:3:12]);
legend([l1 l2],num2str(NS_Fs),num2str(NS_Tb),'location','northwest')



