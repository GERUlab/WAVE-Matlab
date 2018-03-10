function [allam,beta] = calc_thermal_conduc2(units,ncomp,soil_parameters, solute_param1,WC,wcp)

%Calculates thermal conductivity based on Campbell

%In :
% units: specify which units are used in in_general_data
% ncomp: number of compartments
% soil_parameters(:,2) :saturated water conten
% solute_param1(:,1) : bulk density
% WC: water content of the current time step
% wcp: water content of the previous time step

%OUT:
% allam: thermal conductivity between the nodes, halfway the time step 
%in (J L-1 T-1 C-1)
% beta : heat capacity of the soil (J/C/L^3)

% CALLS:
% In_temperature_param
% calc_temp_k
% calc_temp_g_air
% CALLED BY:
% solve_temperature_JV


%From other modules:
wcs = soil_parameters(:,2); %saturated water content
bd = solute_param1(:,1); %bulk density
bd(ncomp+1) = bd(ncomp);
WC = WC'; %Water content

%Define extra water contents for calculation allam (ncomp +1) 
%Joachim V. (Not included in fortran)
WC(ncomp +1)= WC(ncomp);
wcp(ncomp +1) = wcp(ncomp);
wcs(ncomp+1) = wcs(ncomp);


    %%calculate water content halfway the node, halfway the time step
    i = [2:ncomp+1];
    tw(1) = (WC(1) + wcp(1))/2;   
    tw(i) = (WC(i) + wcp(i) + WC(i-1) + wcp(i-1))/4;
    tw = tw';

%User defined input of temperature parameters
[A,B,C,D,E,bd_water,heatcapsolids_corr, heatcapwater_corr,disp]= In_temperature_param2(units,ncomp);

%Calculate volumetric heat capacity (eq 4-2 Wave manual)
i = [1:ncomp+1];
beta(i) = heatcapsolids_corr*bd(i) + heatcapwater_corr*tw(i)*bd_water; %in J/C/L^3 
beta = beta';
Cw(i)= heatcapwater_corr*tw(i)*bd_water; %in J/C/L^3 
Cw = Cw';

lo_campbell = A + B*tw-(A-D)*exp(-(-C*tw).^E); 

%Conversion of lo
if strcmp(units{1,1},'cm')==1 
    corr_L = 1/100;
elseif strcmp(units{1,1},'m')==1
    corr_L = 1;
elseif strcmp(units{1,1},'mm')==1 
    corr_L = 1/1000;
end

if strcmp(units{1,2},'day')==1 
    corr_T = (60*60*24);
elseif strcmp(units{1,2},'h')==1
    corr_T = (60*60);
elseif strcmp(units{1,2},'min')==1 
    corr_T = 60;
elseif strcmp(units{1,2},'s')==1 
    corr_T = 1;
end

allam = lo_campbell*corr_T*corr_L;
