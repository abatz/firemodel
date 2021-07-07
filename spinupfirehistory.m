function [firehistory,rat]=spinupfirehistory(BA,option,totalarea,nmem);

% MFI averaged for western US forest is ~150yr
firehistory=1/150*ones(nmem,1);

for i=2:length(BA)
    firehistory(1:nmem-1,i)=firehistory(2:nmem,i-1);
    firehistory(nmem,i)=BA(i-1)/totalarea;
end


% next use some function to reduce available area to burn as a function of
% fire history

[rat,addrat,firehistory]=firehistorybufferadv(totalarea,firehistory,option,nmem,zeros(1,length(BA)),zeros(1,length(BA)),0,0);
