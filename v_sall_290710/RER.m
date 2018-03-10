function  [rtex,EPRA,esa] = RER(t,ph,simplant,dx,units,plant_date,harvest_date);

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

[AKc,epa,esa,drz,Smax_param,Smax_given] = In_ETsplit(t,simplant,units,plant_date,harvest_date);
ncs=length(ph);%nber of comp. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p0=Smax_param(1);
p1= Smax_param(2);
p2l= Smax_param(3);
p2h=Smax_param(4);
p3=Smax_param(5); wilting_point=p3;
ilinr= Smax_param(6);

if Smax_given;
     Smax=Smax_param(7:end);
else
    arer= Smax_param(7);
    brer= Smax_param(8);
end

if simplant & t>=plant_date & t<harvest_date;     
irz=ceil(drz/dx);

%calculate maximal root extraction rate (Smax) if Smax is not imput
    if Smax_given==0
        Smax(1:ncs,1)=0;
        i=1;
        while(i<=irz)
            heigth=dx*i;  
            Smax(i,1)=arer-(brer.*heigth);           %%%% (day-1)
            i=i+1;
        end
        if irz>drz/dx;
            Smax(irz,1)=0; 
        end
   end
%IN WAVE.FOR
%Calculate the sum of maximal root extraction rate
    sum_smax=sum(Smax);
%Calculate fraction of root distribution for each compartment
    fraction=Smax/sum_smax;
% recalculate smax to avoid the depassment of epa
    Smax=fraction*epa/dx;     %%%ESSAI A ENLEVER
 
    %if irz>drz/dx;
    %    Smax(irz,1)=fraction*epa/(irz*dx-drz);
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
    if irz < ncs
      p2=p2h;
        if epa < 0.1
        p2=p2l;
        elseif epa < 0.5
            p2=p2h+((0.5 - epa)/0.4)*(p2l-p2h);
            
        else
        end
       
        for i=1:irz
            if ph(i) >= p1
                if p0 >= p1
					alpha_h(i)= max((p0-ph(i))/(p0-p1),0.0);
                else
					alpha_h(i) = 0;
                end
            end
            if ph(i) < p2
                if ph(i) > p3
                    if p3 < p2
                        if ilinr 
                            alpha_h(i)= max ((p3-ph(i))/(p3-p2),0);
                        else
							alpha_h(i)=10^((p2-ph(i))/p3);
                        end
                    else
                        alpha_h(i) = 0.;
                    end
                else
                    alpha_h(i) = 0.;
                end
            else
				alpha_h(i) = 1.;
            end
        end
       %		root extraction rate is also limited to avoid extraction below
%		wilting point (wcwp)
        %for i = 1: irz
		%wcwp(i)=wilting_point;
        %end
        for i = 1: irz
		rtex(i)= max(0., alpha_h(i)*(Smax(i))); %(wc(i)-wcwp(i))/dt) );
        end

%		rtex(i)= 0 for compartments below the root zone
        for i = irz+1: ncs
		rtex(i) = 0.;
        end
        
    end
%Determination of rtex for each compartiment

    %rtex=alpha_h.*Smax;   %%%%%%%%  (t-1)
    EPRA=0;
    for i=1:ncs
        EPRA=EPRA+rtex(i);
    end

elseif simplant==0 |t<plant_date |t>=harvest_date;
    rtex=zeros(ncs,1);
    EPRA=0;
else 
    rtex=zeros(ncs,1);
    EPRA=0;

end

