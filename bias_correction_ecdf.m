    function [gcm_obs,gcm_f_t]=bias_correction_ecdf...
        (obs,gcm_h,gcm_f);

    %Preset output array
    temp_out = NaN(size(gcm_h));
    temp_out_f = NaN(size(gcm_f));

    %==============================================
    %Get empirical cdf of gcm historical data
    %==============================================
    %get empirical cdf of obs
    [cdf_obs,xobs]=ecdf(obs);
    [cdf_gcm,xgcm]=ecdf(gcm_h);
    [cdf_gcm_f,xgcm_f]=ecdf(gcm_f);
    
    %make cdf estimates distinct for the cdf obs lookup table
    [xobs_t,index_num]=unique(xobs);
    cdf_obs_t =cdf_obs(index_num);
    %modify distribution so goes from [0,1]
    N=length(cdf_obs)-1; %ecdf always duplicates 1/N twice
    cdf_obs_t = (N*cdf_obs_t - 1)/(N-1);
    if(cdf_obs_t(1)~=0)
        cdf_obs_t(1) = 0.0;
    end

    
    %make cdf estimates distinct
    [xgcm_t,index_num]=unique(xgcm);
    cdf_gcm_t =cdf_gcm(index_num);
    %modify distribution so goes from [0,1]
    N=length(cdf_gcm)-1; %ecdf duplicates 1/N twice
    cdf_gcm = (N*cdf_gcm_t - 1)/(N-1);
    if(cdf_gcm_t(1)~=0)
        cdf_gcm_t(1) = 0.0;
    end
    cdf_gcm =interp1  (xgcm_t,cdf_gcm_t,gcm_h);
    
    %make cdf estimates distinct
    [xgcm_f_t,index_num]=unique(xgcm_f);
    cdf_gcm_f_t =cdf_gcm_f(index_num);
    %modify distribution so goes from [0,1]
    N=length(cdf_gcm_f)-1; %ecdf duplicates 1/N twice
    cdf_gcm_f = (N*cdf_gcm_f_t - 1)/(N-1);
    %make cdf estimates distinct
    if(cdf_gcm_f_t(1)~=0)
        cdf_gcm_f_t(1) = 0.0;
    end
    cdf_gcm_f =interp1  (xgcm_f_t,cdf_gcm_f_t,gcm_f);
   

    %map gcm values onto obs
        gcm_obs=interp1(cdf_obs_t,xobs_t,cdf_gcm);
        gcm_f_obs=interp1(cdf_obs_t,xobs_t,cdf_gcm_f);
   
    %==============================================
    %PRESERVE QUANTILE DIFFERENCES FUTURE TO HISTORICAL (EPOCH ADJUSTMENT)
    %==============================================
    %interpolate the cdf values of future onto the historical cdf
   
        xgcm_f_t_gcm=interp1  (cdf_gcm_t,xgcm_t,cdf_gcm_f_t);
  
    %find difference between the gcm values  (gcm - gcm_f)
    future_offset = xgcm_f_t- xgcm_f_t_gcm;    %corresponding y value is cdf_gcm_f


    cdf_offset=cdf_gcm_f_t;   %already unique
    factor_fut=interp1(cdf_offset,future_offset,cdf_gcm_f);
    gcm_f_t =  gcm_f_obs + factor_fut;
    [a,b]=sort(gcm_f);
    [a2,b2]=sort(gcm_f_t);
    gcm_f_t(b)=gcm_f_t(b2);

