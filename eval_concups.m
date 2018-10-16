function ConcUps = eval_concups_T(ConcLoc,DecayTime,reach_ind,other)

[all_covariates,reach_upstream,length_downstream,area_local,area_upstream,celerity,length_reach,reach_width,Q0]=v2struct(other);

for ind=1:numel(reach_ind)
    j=reach_ind(ind);
    subset=reach_upstream(reach_upstream(:,j)>0,j); % find subset of reaches upstream of j
    Qbar(j)=area_upstream(j)/area_upstream(13)*Q0; % median value of discharge; '13' is the position of the stretch where Q0 was measured
    ConcUps(ind,1)=1e-3./Qbar(j).*...
        sum(exp(-(length_downstream(subset,j)-length_downstream(j,j))./(celerity(subset,j)*DecayTime)).*...
    length_reach(subset).*reach_width(subset).*ConcLoc(subset));        
    % note that the factor 1e-3 corrects the discrepancy between units
    % (from L at the denominator of the l.h.s. to m^3 at the r.h.s.)
end


end

