function ph = In_initial_pressure_head(compartiments_number,soil_parameters,dx_inter,...
    BOUNDARY_CONDITIONS_MATRIX,initial_gwl)

%Initial pressure head
% if i=
%		1) Initial pres. head profile calculated from top and bottom pres. head values: i=1
%		2) Initial pres. head profile calculated from some pres. head values given for some depths: i=2
%		3) Uniform initial pres. head profile: i=3
%		4) Initial pres. head profile calculated from initial water content given for some depths
%       5) Initial pres. head profile calculated from equilibrum with groundwater table
% IN:
%	compartiments_number,soil_parameters, dx_inter,BOUNDARY_CONDITIONS_MATRIX,initial_gwl 
%   dx_inter
% OUT:
%	ph: initial pressure head 
% CALLS:
%	none


i=3;%Fill it

if i==1;
   phbot = -40.;	
   phtop = -41.;	
   ph=phtop:(phbot-phtop)/(compartiments_number-1):phbot;
  
elseif i==2;
   ph=[-90:-5/29:-95 -95:70/69:-25];
  
elseif i==3;

    ph(1:compartiments_number)=-100;
   
elseif i==4;
   zdata=[0 2.75 5.5 8.25 11 13.75 16.5 19.25 22];		%depths (cm) for which moisture contents are available
   wcdata=[0.22 0.21 0.09 0.09 0.09 0.12 0.22 0.27 0.26];    %Moisture content for each depth
   z=[1:1:compartiments_number]
   WC=spline(zdata,wcdata,z);%MAMADOU
   WCR=soil_parameters(:,1);
   WCS=soil_parameters(:,2);
   ALFA=soil_parameters(:,3);
   N=soil_parameters(:,4);
   M=1-1./N;
   Se=(WC'-WCR)./(WCS-WCR);
   ph=-(Se.^(-1./M)-1).^(1./N)./ALFA;ph=ph'; %origine
   
else  %if Equilibrum with a ground water table)
    ncs=compartiments_number;
    gwl=initial_gwl;
    ph=-abs(gwl)+abs([1:ncs].*dx_inter(1:ncs));
end

