function  [rtex,EPRA] = RER(t,ph,simplant,dx,units,plant_date,harvest_date);

%RER ::  Calculation of the extraction profile
%
%IN:
%   t: current time
%   ph: ph distribution
%   simplant: plant simulation (logical)
%   dx: spatial discretization
%OUT:
%   rtex: root extraction rate (1/day)
%   EPRA: total root extraction (1/day)
%CALLS:
%   In_ETsplit.m
%CALLED BY:
%   state_var.m
%--------------------
% F. Hupet 15/10/2000
% Javaux, 2006
%Modified by M. Sall 21/11/2009

[AKc,EPA,esa,DRZ,ARER,BRER] = In_ETsplit(t,simplant,units,plant_date,harvest_date);
ncs=length(ph);%nber of comp. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate maximal root extraction rate (Smax)  


if simplant & t>=plant_date & t<=harvest_date;  
irz=ceil(DRZ/dx);
Smax(1:ncs,1)=0;
i=1;
    while(i<=irz)
        heigth=dx*i;  
        Smax(i,1)=ARER-(BRER.*heigth);           %%%% (day-1)
        i=i+1;
    end
    if irz>DRZ/dx;
        Smax(irz,1)=0; 
    end
  
%IN WAVE.FOR
%Calculate the sum of maximal root extraction rate
    sum_smax=sum(Smax);
%Calculate fraction of root distribution for each compartment
    fraction=Smax/sum_smax;
% recalculate smax to avoid the depassment of EPA
    Smax=fraction*EPA/dx;     %%%ESSAI A ENLEVER
 
    %if irz>DRZ/dx;
    %    Smax(irz,1)=fraction*EPA/(irz*dx-DRZ);
    %end

    i=1;
    while(i<=ncs)
        if (Smax(i,1)<0);
            Smax(i,1)=0;
        else
            Smax(i,1)=Smax(i,1);
            i=i+1;
        end
    end

%Reduction of potential root water uptake
    alpha_h(1:ncs,1)=0;
    rtex(1:ncs,1)=0;

%Determination of alpha_h for each compartiment
    i=1;
    while (i<=ncs)
        if (abs(ph(i))<50)
            alpha_h(i)=0;
        elseif (abs(ph(i))>=50)&(abs(ph(i))<=100);
            alpha_h(i)= 1*((abs(ph(i))-50)/(100-50));
        elseif (abs(ph(i))>=100)&(abs(ph(i))<=600)
            alpha_h(i)=1;
        elseif (abs(ph(i))>=600)&(abs(ph(i))<=15000)
            alpha_h(i)=1*((15000-abs(ph(i)))/(15000-600));
        else
            alpha_h(i)=0;
        end 
   i=i+1;
    end
alpha_h(1:ncs)=1;   %%%ESSAI A ENLEVER
%Determination of rtex for each compartiment

    rtex=alpha_h.*Smax;   %%%%%%%%  (t-1)
    EPRA=0;
    for i=1:ncs
        EPRA=EPRA+rtex(i);
    end

elseif simplant==0 |t<plant_date |t>plant_date;
    rtex=zeros(ncs,1);
    EPRA=0;
else 
    rtex=zeros(ncs,1);
    EPRA=0;

end

%PRECEDENT CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EPRA=0;
%i=1;
%while(i<=ncs)
%   EPRA=EPRA+rtex(i,1);  
%   if (EPRA<=EPA)               %%%Remark: EPRA must be divided by dx before comparing withe EPA with respect to dimension
 %     rtex(i,1)=rtex(i,1);
 %  elseif (EPRA>=EPA);
%      rtex(i:1:ncs,1)=0;
%   end
%   i=i+1;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

