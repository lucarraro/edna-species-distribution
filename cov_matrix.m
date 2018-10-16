function [covariates,all_covariates]=cov_matrix(Cov_labels,N_cov,reach_ind,mean_subcatchment_altitude,area_upstream,reach_geology_class)
% function that builds the matrix of covariates

% normalization constants for covariates
minAlt=min(mean_subcatchment_altitude); maxAlt=max(mean_subcatchment_altitude); avgAlt=0.5*(minAlt+maxAlt);
minAreaUp=min(area_upstream); maxAreaUp=max(area_upstream); avgAreaUp=0.5*(minAreaUp+maxAreaUp);
minGeoMor=min(reach_geology_class(:,3)); maxGeoMor=max(reach_geology_class(:,3)); avgGeoMor=0.5*(minGeoMor+maxGeoMor); 
minGeoPea=min(reach_geology_class(:,4)); maxGeoPea=max(reach_geology_class(:,4)); avgGeoPea=0.5*(minGeoPea+maxGeoPea);  
minGeoWat=min(reach_geology_class(:,5)); maxGeoWat=max(reach_geology_class(:,5)); avgGeoWat=0.5*(minGeoWat+maxGeoWat); 

for i=1:N_cov
    if strcmp(Cov_labels(i,:),'LocElv')
        covariates(:,i)=(mean_subcatchment_altitude(reach_ind)-avgAlt)/(maxAlt-avgAlt);
        all_covariates(:,i)=(mean_subcatchment_altitude(:)-avgAlt)/(maxAlt-avgAlt);
    elseif strcmp(Cov_labels(i,:),'UpCAr ')
        covariates(:,i)=(area_upstream(reach_ind)-avgAreaUp)/(maxAreaUp-avgAreaUp);
        all_covariates(:,i)=(area_upstream(:)-avgAreaUp)/(maxAreaUp-avgAreaUp);
    elseif strcmp(Cov_labels(i,:),'GeoMrn')
        covariates(:,i)=(reach_geology_class(reach_ind,3)-avgGeoMor)/(maxGeoMor-avgGeoMor);
        all_covariates(:,i)=(reach_geology_class(:,3)-avgGeoMor)/(maxGeoMor-avgGeoMor);
    elseif strcmp(Cov_labels(i,:),'GeoPea')
        covariates(:,i)=(reach_geology_class(reach_ind,4)-avgGeoPea)/(maxGeoPea-avgGeoPea);
        all_covariates(:,i)=(reach_geology_class(:,4)-avgGeoPea)/(maxGeoPea-avgGeoPea);
    elseif strcmp(Cov_labels(i,:),'GeoWat')
        covariates(:,i)=(reach_geology_class(reach_ind,5)-avgGeoWat)/(maxGeoWat-avgGeoWat);
        all_covariates(:,i)=(reach_geology_class(:,5)-avgGeoWat)/(maxGeoWat-avgGeoWat);
        
    end
end


end