function [temp,temptopp]=In_initial_temperature(ncomp)

%Initial temperature profile 
%IN
%ncomp=number of compartiments;
%OUT
%tempp:initial temperature for all compartiments
%temptopp= soil surface temperatute at time t=0
%--------------------------------------------------------------------------
%M.SALL 26/11/08
%J. Vansteenkiste 18/09/09

j=3;

if j==1
temptopp=20;                 %soil surface temperature
tempbotp= 7;                 %bottom temperature
for i=1:ncomp
    temp(i)=i*(tempbotp-temptopp)/ncomp +temptopp;
end

elseif j==2
    temp(1:ncomp)=7;
    temptopp=0;
    
elseif j ==3
    firstcomp = 7;
    lastcomp = 7;
    for i=1:ncomp
    temp(i)=i*(lastcomp-firstcomp)/ncomp +firstcomp;
end
    temptopp=0;
end