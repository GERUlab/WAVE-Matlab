function [b1,b2,b3,bd_water,heatcapsolids_corr, heatcapwater_corr,disp]= In_temperature_param3(units,ncomp)

disp = 0;
%Adopted from the paper by Chung and Horton (units W/m/°C)
%Clay
% b1 = -0.197;
% b2 = -0.962;
% b3 = 2.2521;
%Loam
b1 = 0.243;
b2 = 0.393;
b3 = 1.534;
%Sand
% b1 = 0.228;
% b2 = -2.406;
% b3 = 4.909;

heatcapsolids = 840;                    % specific heat capacity solids (J/kg/C)
heatcapwater = 4186.8;                  % specific heat capacity water (J/kg/C)	

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

%Conversion of b factors
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

b1 = b1*corr_T*corr_L;
b2 = b2*corr_T*corr_L;
b3 = b3*corr_T*corr_L;