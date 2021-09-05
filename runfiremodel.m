function [BA,rat,permrat,firehistory]=runfiremodel(totalarea,firehistory,option,nmem,futureZ,Z5,addrat,a2,residual,mcap,alphaval,capannual);

[rat,addrat,firehistory]=firehistorybufferadv(totalarea,firehistory,option,nmem,futureZ,Z5,addrat,alphaval,nyeardrought);

        
    BA=10.^(predict(a2,futureZ.*(1-rat)));

    % could also be some limit such that the burned area annually can not
    % exceed more than X% of existing land
    BA=min(BA,(1-rat)*totalarea*mcap);
    permrat=addrat;
    firehistory(end+1)=BA/totalarea;
    firehistory=firehistory(2:end);
