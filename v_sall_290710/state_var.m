function [WC,kh,CH,rtex,EPRA,esa]=state_var(ph,soil_parameters,dt,update_historic,t,simplant,dx,units,...
    plant_date,harvest_date);

%STATE_VAR  Calculates state variables (moisture retention, hydraulic conductivity, differential 
%   moisture capacity and root extraction rates for a given pressure head)
%
%IN:
%   ph,soil_parameters,dt,update_historic,t,simplant,dx
%OUT:
% 	WC = volumetric soil moisture content
%   kh = soil unsaturated hydraulic conductivity (cm/min)
%   CH = differential moisture capacity (1/cm)
% 	rtex = root extraction rate (1/min)
%   EPRA=plant parameter
%CALL:
%	moist_ret,conduct,diff_moist_capac,RER
%CALLED BY:
%   solve_flow.m
%---------------------------
% M. Vanclooster 2/2/2000

WC = moist_ret(ph,soil_parameters,dt,update_historic);
kh=conduct(ph,soil_parameters);
CH=diff_moist_capac(ph,soil_parameters);
%if simplant==1  %%MODIF 24/03/10 MSALL
   [rtex,EPRA,esa] = RER(t,ph,simplant,dx,units,plant_date,harvest_date);
%else             %%MODIF 24/03/10 MSALL
%   rtex=0.*ph';   %%MODIF 24/03/10 MSALL
 %  EPRA=0.*ph';    %%MODIF 24/03/10 MSALL
%end                  %%MODIF 24/03/10 MSALL
rtex=rtex';
