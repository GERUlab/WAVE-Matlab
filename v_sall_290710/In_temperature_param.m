function [ttop,tbot,rocp,lamda]=In_temperature_param(t,WC,wcs,...
    ncs,ttopp,tbotp);

%MODULE TEMPERATURE
%IN_TEMPERATURE_PARAM calculate top and bottom temperature and parameters 
%IN:
%t=time(mn)
%WC=water content(-); wcs=saturated water content

%ttopp= soil surface temperatute at time t-dt
%tbotp= soil bottom temperatute at time t-dt
%ncs=compartiment number
%OUT:
%ttop= soil surface temperatute at time t
%tbot= soil bottom temperatute at time t
%rocp,lamda:  parameters at time t
%CALL:none
%CALLED BY: SOLVE_TEMPERATURE
%--------------------------------------------------------------------------
%M.SALL 26/11/08

%calculate soil surface temperature;
tempmax=30; tempmin=4;
gamma=tempmax-tempmin;
ttop=ttopp+gamma*sin(2*pi*t/1);
tbot=tbotp;
%calculate parameters
for i=1:ncs;
ro_s= 2.7;              %(bulk density of soil solid)
Cs=840;              %(spesific gravimetric heat capacity of the soil solid J/Kg/°C)
ro_w=1;              %(density of water)
C_w=4.2;             % (spesific gravimetric heat capacity of water J/Kg/°C)
rocp(i)=ro_s*Cs+WC(i)*ro_w*C_w;   %volumetric heat capacity of the soil

%volume fraction X of different components
X_quartz=0.54;
X_OM=0.045;
X_other=0.015;
X_water=WC(i);
X_air=1-X_quartz-X_OM-X_OM-WC(i);
X=[X_quartz;X_OM;X_other;X_water;X_air];
% thermal conductivity of the different components
lamda_quartz=20.4;lamda_OM=0.6;lamda_other=7;lamda_air=0.0615;
lamda_water=0.176;
lamda_component=[lamda_quartz;lamda_OM;lamda_other;lamda_air;...
lamda_water];

    if WC(i)>=0.2;
  lamda_v(i)=0.176; lamda_component(4)=lamda_air+lamda_v(i); lamda_0=lamda_water;
    else
% lamda decrease linearly from 0.176 to 0 (between WC=0.2 and oven dry)
   lamda_v(i)=0.176*(WC(i)-0.2)/0.2; lamda_component(4)=lamda_air+lamda_v(i);
   lamda_0=lamda_air;
    end

% value of g factors for the different conponents 
g_quartz=[0.125,0.125,0.75];
g_OM=[0.5,0.5,0];
g_other=[0.125,0.125,0.75];
%g_water is not explicitly defined in the WAVE_MANUAL (cas WC<0.2)
g_water=[1/3,1/3,1/3];
%calculate g_air witch decrease linearly from 0.333 to 0.105 (between wcs 
%and  WC=0.2) and from 0.105 to 0.015 (between WC<0.2 and oven dry)
    if WC(i)>=0.2;
g_air(1)=(WC(i)-wcs(i))*(0.105-0.333)/(0.2-wcs(i)) +0.333;
g_air(2)=g_air(1); g_air(3)=1-g_air(1)-g_air(2);
    else
g_air(1)=(WC(i)-0.2)*(0.015-0.105)/(0-0.2) +0.105;
g_air(2)=g_air(1); g_air(3)=1-g_air(1)-g_air(2);
    end

g=[g_quartz; g_OM; g_other; g_air; g_water];

% calculate the value of k coefficients for the different components

    for I=1:length(lamda_component);
        if I== 1;
K(1)=0; 
        else
    K(I,1)=1/3*sum(1./(1+((lamda_component(I)/lamda_0) -1).*g(I,:)));
        end
        if WC(I)>=0;
    c=1.65;
        else
    c=1;
        end
    end
%CALCULATE global lamda OF THE i_th compartiment of the SOIL   
lamda(i)=c*sum(K.*X.*lamda_component)/sum(K.*X);
end




