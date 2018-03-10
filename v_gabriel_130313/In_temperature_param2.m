function [A,B,C,D,E,bd_water,heatcapsolids_corr, heatcapwater_corr,disp]= In_temperature_param2(units,ncomp)

%As by Campbell
solid = 0.59; 
organicmatter = 0; 
disp = 5; %L
quartz = 0.54;
other = 0.045;
clay = 0.015;

heatcapsolids = 840;                    % specific heat capacity solids (J/kg/C)
heatcapwater = 4186.8;                  % specific heat capacity water (J/kg/C)	


A = (0.57 + 1.73*quartz + 0.93*other)/(1-0.74*quartz-0.49*other)-2.8*solid*(1-solid);
B = 2.8*solid;
C = 1+2.6*(clay^(-0.5));
D = 0.03+0.7*solid^2;
E=4;

%%End of inputsection, all input parameters are converted to the correct
%%unit and/or put in a matrix form

%Conversion of heat capacity to J/M/C
if strcmp(units{1,3},'kg')==1 
    corr_M = 1;
elseif strcmp(units{1,3},'g')==1
    corr_M = 1/1000;
elseif strcmp(units{1,3},'mg')==1 
    corr_M = 1/1000000;
end
heatcapsolids_corr = heatcapsolids*corr_M; %in J/M/C
heatcapwater_corr = heatcapwater*corr_M; %in J/M/C

%Conversion bulk density water to M L^-3
if strcmp(units{1,1},'cm')==1 
    corr_L = 1/100^3;
elseif strcmp(units{1,1},'m')==1
    corr_L = 1;
elseif strcmp(units{1,1},'mm')==1 
    corr_L = 1/1000^3;
end
bd_water = 1000*corr_L/corr_M;
%1000 is bulk density in kg/m^3

