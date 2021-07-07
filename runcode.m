[obsdata]=csvread('observed_climate_fire.csv',1,0);
Z=obsdata(:,5);

BA=obsdata(:,2);

% truncate to 1984-2020 when BA data is available
yr=1984-1949:2020-1949;
BA=BA(yr);
Z=Z(yr);
%area of forest extent from LANDFIRE ESP
totalarea=833330;

% reference fuel memory (years)
nmem=30;

% number of CMIP6 models
nmod=30;

% annual max burned fraction
capannual=0.075;

[wus]=csvread('westernusforest_extension.csv',1,0);
wus=wus(:,2);
[cmip6]=csvread('cmip6_model_cwd.csv',1,0);
% years 1950-2100 only
cmip6=cmip6(101:251,2:31);

% bias correct output

[cmip6]=biascorrect(Z,cmip6,nmod);



% loop over cmip models
for model=1:nmod
    % loop over fire-fuel feedback models
    % option 1= static
    % option 2= weak-constant
    % option 3= moderate-constant
    % option 4= strong-constant
    % option 5= weak-fading
    % option 6= moderate-fading
    % option 7= strong-fading

for option=1:7
    % spin up historical fuel estimates based on fire activity
[firehistory,rat1(:,option)]=spinupfirehistory(wus,option-1,totalarea,nmem);
fire1950=firehistory(:,1950-1916);
firehistory=firehistory(:,end-36);
[r,p,FFA,ra,addr,raold,permold]=ba_feedbackmodel(BA,Z,totalarea,cmip6(:,model),option-1,capannual,nmem,firehistory,fire1950);
if model==1 
   rat1(69:105,option)=raold+permold;
end
ratg(:,:,option,model)=ra;
rval(:,:,option,model)=r;
burned(:,:,option,model)=FFA;
arat(:,:,option,model)=addr;
oldrat(option,:,:)=raold;
end
end


