function  [rtex,EPRA] = RER(t,ph,simplant,dx);

%Calculation of the extraction profile
%
%IN:
%   t: current time
%   ph: ph distribution
%   simplant: plant simulation (logical)
%   dx: spatial discretization
%OUT:
%   rtex: root extraction rate (1/min)
%   EPRA: total root extraction (1/min)
%CALLS:
%   In_ETsplit.m
%CALLED BY:
%   state_var.m
%--------------------
% F. Hupet 15/10/2000
% Javaux, 2006

[AKc,EPA,esa,DRZ,ARER,BRER] = In_ETsplit(t,simplant);
ncs=length(ph);%nber of comp. 

Smax(1:ncs,1)=0;

i=1;
while(i<=ncs)
   heigth=dx*i;
   Smax(i,1)=ARER-(BRER.*heigth);
   i=i+1;
end

i=1;
while(i<=ncs)
   if (Smax(i,1)<0)
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
   elseif (abs(ph(i))>=50)&(abs(ph(i))<=100)
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

%Determination of rtex for each compartiment

rtex=alpha_h.*Smax;
EPRA=0;

i=1;

while(i<=ncs)
   EPRA=EPRA+rtex(i,1);
   if (EPRA<=EPA)
      rtex(i,1)=rtex(i,1);
   elseif (EPRA>=EPA)
      rtex(i:1:ncs,1)=0;
   end
   i=i+1;
end


