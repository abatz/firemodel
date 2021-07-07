function [r,p,BA,ra,permrat,ratold,permratold]=ba_feedbackmodel(BA,Z,totalarea,futureZ,option,capannual,nmem,firehistory,firehistory1950);

% experimental function that builds in a negative feedback to burned area
nensemble=1000;
Z(6:end+5)=Z;
Z5=movmean(Z,[5 0]);
Z5=Z5(6:end);
Z=Z(6:end);


if nargin==7 % if you do not have fire history layer
% SPINUP: assume some fraction of previous say 30-yrs have burnhistory
% default it 0.1% of lands burn
firehistory=0.001*ones(nmem,1);

for i=2:length(BA)
    firehistory(1:nmem-1,i)=firehistory(2:nmem,i-1);
    firehistory(nmem,i)=BA(i-1)/totalarea;
end
end

% reduce available area to burn as a function of fire history
addrat=0;
for i=1:length(BA)
 [rat,addrat,firehistory]=firehistorybufferadv(totalarea,firehistory,option,nmem,Z(i),Z5(i),addrat,0.1);
 firehistory(end+1)=BA(i)/totalarea;
 firehistory=firehistory(2:end);
 rrr(i)=rat+addrat;
 permoff(i)=addrat;
end
rat=rrr;
Z2=Z;


% statistics about the model fit
a2=fitlm((Z2'),log10(BA./(1-rat)'));
mm=predict(a2,Z2);mm=(10.^mm)./(1-rat)';
[r(1),p(1)]=corr(log10(BA./(1-rat)'),Z2);
[r(2),p(2)]=corr(detrend(log10(BA./(1-rat)')),detrend(Z2));
[r(3),p(3)]=corr(log10(BA(1:16))./(1-rat(1:16)'),Z2(1:16));
[r(4),p(4)]=corr(log10(BA(17:37))./(1-rat(17:37)'),Z2(17:37));

% residuals and constrained residuals
residual_log=(log10(BA)-log10(mm));
r1=prctile(residual_log,[15 85]);
residual_log(residual_log<r1(1))=r1(1);
residual_log(residual_log>r1(2))=r1(2);

% save permoff firehistory and rrr
permratold=permoff;
ratold=rat;


reset1950=0;
if reset1950
% if you want to reset burned area starting at 1950
firehistory=firehistory1950;
clear BA rat
startyr=1;
else
    startyr=72;
end

if option==0 startyr=1;end  % this is for no-feedback model

% starting in 1950 iterative forward
    Z5(6:156)=movmean(futureZ(1:151),[5 0]);
    addrat=0;
    
    initfhist=firehistory;
    for kk=1:nensemble
        firehistory=initfhist;addrat=0;
for i=startyr:151
    [rat,addrat,firehistory]=firehistorybufferadv(totalarea,firehistory,option,nmem,futureZ(i),Z5(i),addrat,.1);
     BA(i,kk)=10.^(predict(a2,(futureZ(i)))+residual_log(ceil(rand(1,1)*37))).*(1-rat);
    
    % limit such that the burned area annually can not
    % exceed more than X% of existing land
    BA(i,kk)=min(BA(i,kk),(1-rat)*totalarea*capannual);
    ra(i,kk)=rat;
    permrat(i,kk)=addrat;
    firehistory(end+1)=BA(i,kk)/totalarea;
    firehistory=firehistory(2:end);

end
end   