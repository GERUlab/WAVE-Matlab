function [rocpp,tempp,lamdap,ttopp,tbotp,wcs]=In_initial_temperature(ncs,soil_parameters);
%MODULE TEMPERATURE
%IN_INITIAL_TEMPERATURE: initial soil temperature
%IN
%ncs=number of compartiments
%soil parameters=soil parameters used for the water flow module
%OUT
%tp:initial temperature for the different compartiments
%ttop= soil surface temperatute at time t=0
%tbot= soil bottom temperatute at time t=0
%rocp,lamda:  parameters at time t=0
%saturated water content
%--------------------------------------------------------------------------
%M.SALL 26/11/08


%% Initialisation de ROCP,LAMDA et T
rocpp=zeros(ncs,1);       %volumetric heat capacity of the soil   
lamdap=zeros(ncs,1);      %thermal coductivity
ttopp=18;                   %soil surface temperature
tbotp=7;                    %bottom temperature
wcs=soil_parameters(:,2); 
for i=1:ncs
    tempp(i)=i*(7-18)/ncs +18;
end