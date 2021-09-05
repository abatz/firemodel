function [cmip6]=biascorrect(Z,cmip6,nmod);

% first remove 31-yr moving running mean of data, we'll add this back at
% the end
t=movmean(cmip6,31,1);
cmip6=cmip6-t;

% loop over models: run cdf for each 30-40 yr interval for safe measure
for model=1:nmod
    [Z1(:,model),ZZ1(1:41,model)]=bias_correction_ecdf...
        (Z(end-29:end),cmip6(42:71,model),cmip6(1:41,model));
    [Z1(:,model),ZZ1(42:71,model)]=bias_correction_ecdf...
        (Z(end-29:end),cmip6(42:71,model),cmip6(42:71,model));
    [Z1(:,model),ZZ1(72:101,model)]=bias_correction_ecdf...
        (Z(end-29:end),cmip6(42:71,model),cmip6(72:101,model));
    [Z1(:,model),ZZ1(102:151,model)]=bias_correction_ecdf...
        (Z(end-29:end),cmip6(42:71,model),cmip6(102:151,model));
end;
ZZ1(1:41,:)=ZZ1(1:41,:)-repmat(mean(ZZ1(1:41,:),1),[41 1]);
ZZ1(42:71,:)=ZZ1(42:71,:)-repmat(mean(ZZ1(42:71,:),1),[30 1]);
ZZ1(72:101,:)=ZZ1(72:101,:)-repmat(mean(ZZ1(72:101,:),1),[30 1]);
ZZ1(102:150,:)=ZZ1(102:150,:)-repmat(mean(ZZ1(102:150,:),1),[49 1]);


ZZ1=ZZ1+t;

%finally standardize to ensure than year=42-71 for each model have mean of 0
ZZ1=ZZ1-repmat(mean(ZZ1(42:71,:),1),[150 1]);
cmip6=ZZ1;
