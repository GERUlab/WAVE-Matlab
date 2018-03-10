 function  [output] = calc_esa(t,simplant,plant_date,harvest_date);

 %Calculate the potential evaporation
% Script based on In_ETsplit, but adapted to calculate esa for every
% climatic measurement
%Joachim Vansteenkiste, 11/08/2010

 output = [];
%% Plant information on water uptake
climate = in_climatic_data;
time = climate(:,1); 
ETo = climate(:,2);

for i=1:length(climate)
    t = climate(i,1);
boolean_plant = (t>= plant_date);
boolean_harvest = (t<= harvest_date);
boolean_plant_harvest = boolean_plant + boolean_harvest;
index = find(boolean_plant_harvest == 2);

if isempty(index) == 1
    [drz, rlai, kc, Smax_given,Smax_param,root_length_max] = in_crop_number1(t,plant_date, harvest_date,index); %not used in calculations
elseif index == 1
[drz, rlai, kc, Smax_given,Smax_param,root_length_max] = in_crop_number1(t,plant_date, harvest_date,index);
elseif index == 2
[drz, rlai, kc, Smax_given,Smax_param,root_length_max] = in_crop_number2(t,plant_date, harvest_date,index);
elseif index == 3
[drz, rlai, kc, Smax_given,Smax_param,root_length_max] = in_crop_number3(t,plant_date, harvest_date,index);
else
    error('define new function for the fourth crop')
end

                       
 
%% Searches for the appropriate plant water uptake for the current time t 
i=max(find(climate(:,1)<=t));
ETo=climate(i,2); % L/TET
AKc=kc;

if simplant & (boolean_plant_harvest(index) == 2);  %Simulation with crop
   EToA=AKc*ETo;        %Potential evapotranspiration of the crop (L T-1)
   esa=exp(-0.6*rlai)*EToA;		%Potential evaporation of the soil (L T-1)
   epa=EToA-esa;	%Potential transpiration of the crop (L T-1)
else
   AKc=1;
   esa=ETo;     %Potential evaporation of the soil (L T-1)
   epa=0;     %Potential transpiration of the crop (L T-1)
   rlai=0;
   drz=0;
end

output(end+1) =esa;
end
% EVAPORATION OF INTERCEPTED WATER (in WAVE.FOR (RER.FOR))
%ev_intc=min(stor,epa)
%epa_min_intc= epa-av_intc;
%stor=stor-ev_intc;




