function [part,K_i,bd_water,heatcapsolids_corr, heatcapwater_corr,shape_par, shape_par2, K_vapour]= In_temperature_param4(units,ncomp)

%Input parameters for the temperature module as by De Vries

%mainly fixed values
%Parameters are a fixed value,based on literature.
%Components: 
%quartz| organicmatter| other solids| air| water
part = [0.3 0.27 0.03 0 0];          % volume fraction of each component
K_i = [8.8 2.9 0.25 0.025 0.57];    % thermal conductivity of each component W m-1 C-1
shape_par = [1 1 0 999 0];            %Determines the shape of the particles
shape_par2 = [3.5 3.5 3 999 0];       %Determines the location
heatcapsolids = 840;                    % specific heat capacity solids (J/kg/C)
heatcapwater = 4186.8;                  % specific heat capacity water (J/kg/C)
K_vapour = 0.04; %W m-1 C-1


%%End of inputsection, all input parameters are converted to the correct
%%unit and/or put in a matrix form

part = repmat(part,ncomp+1,1);
K_i = repmat(K_i, ncomp+1,1);%thermal conductivity of each component (mCal cm-1 sec-1)

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

%Convert thermal conductivity of W m-1 C-1 to J L-1 T-1 C-1
%Calculate the conversion factor needed for the length
if strcmp(units{1,1},'cm')==1 
    corr_L = 1/100;
elseif strcmp(units{1,1},'m')==1
    corr_L = 1;
elseif strcmp(units{1,1},'mm')==1 
    corr_L = 1/1000;
end

%Calculate the conversion factor needed for the time
if strcmp(units{1,2},'day' )==1
    corr_T = (60*60*24);
elseif strcmp(units{1,2},'hour' )==1
    corr_T = (60*60);
elseif strcmp(units{1,2},'min' )==1 
    corr_T = 60;
elseif strcmp(units{1,2},'s' )==1
    corr_T = 1;    
end

%use the correction factors
K_i = K_i*corr_L*corr_T;
K_vapour = K_vapour*corr_L*corr_T;

