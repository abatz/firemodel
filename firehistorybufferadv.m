function [rat,permrat,firehistory]=firehistorybuffer(totalarea,firehistory,option,nmem,Z,Z5,permrat,max_norecovery);


% Z5 is the mean aridity from the past 5-years, this is used to immediately
% pull a percent of burned area out of rotation (forest --> non-forest)

% max_norecovery: what is the maximum extent of burned area from 6-yrs ago that does not

% scales linear function of aridity in past 5-yrs
%percent_norecovery=(cos(max(min(1,Z5-1)-1,-1)*90*pi/180)+1)/1;
percent_norecovery(Z5<1)=0;
percent_norecovery(Z5>2)=1;
percent_norecovery(Z5>=1 & Z5<=2)=Z5(Z5>=1 & Z5<=2)-1;
removedforest=max_norecovery*percent_norecovery.*firehistory(end-5,:);

% burned area from 6 years ago is modified
addrat=removedforest;
firehistory(end-5,:)=firehistory(end-5,:)-removedforest;

% build this back in at the end
fh=firehistory;
permrat=permrat+addrat;

% Z is aridity in current year, we use this as a modifier to
% reduce the time period of negative feedbacks in L following Parks et al.
% 2018 where Z values of 2 shorten feedbacks by a third

shorten(Z>2)=round(nmem/3);
shorten(Z>0 & Z<=2)=round(nmem/3*Z(Z>0 & Z<=2)/2);
shorten(Z<=0)=0;
%shorten=(cos(max(min(2,Z)-2,-2)*90*pi/180)+1)/2;
firehistory=firehistory(1+shorten:end,:);
nmem=nmem-shorten;

switch option
    case 0, rat=zeros(size(firehistory,2),1);
    case 1  % this makes all land that burned in previous 30-yrs out of bounds
        rat=sum(firehistory,1)'*.5;
    case 2
        % like case 1, but further reduce land available by a factor to
        % account for buffering effects of fire (e.g., fire perimeters
        % stopping/slowing advancing fire fronts)
        rat=sum(firehistory,1)'*1;
    case 3 
        % like case 1, but further reduce land available by a factor to
        % account for buffering effects of fire (e.g., fire perimeters
        % stopping/slowing advancing fire fronts)
        rat=sum(firehistory,1)'*1.5;
    case 4
        % all land burned in last 5 years are off limits, then use a
        % sinusodial function over next 25 years to reintroduce land
        nmem=nmem(1);
        b=(1-cos((1:(nmem-5))*pi/(nmem-5)))/2;b(nmem-4:nmem)=1;
        rat=firehistory'*b'*.5;
    case 5 
        % like case 3, but further reduce land available by a factor to
        % account for buffering effects of fire (e.g., fire perimeters
        % stopping/slowing advancing fire fronts)
        nmem=nmem(1);
        b=(1-cos((1:(nmem-5))*pi/(nmem-5)))/2;b(nmem-4:nmem)=1;
        rat=firehistory'*b'*1;
    case 6 
        % like case 3, but further reduce land available by a factor to
        % account for buffering effects of fire (e.g., fire perimeters
        % stopping/slowing advancing fire fronts)
        nmem=nmem(1);
        b=(1-cos((1:(nmem-5))*pi/(nmem-5)))/2;b(nmem-4:nmem)=1;
        rat=firehistory'*b'*1.5;
    case 7  % linear function of inelegible land based on years since it burned
        b=nmem:-1:1;
        b=1-b/nmem;
        rat=firehistory'*b;
end

firehistory=fh;
