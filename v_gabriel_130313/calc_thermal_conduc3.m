function [allam, beta] = calc_thermal_conduc3(units,ncomp,soil_parameters, solute_param1,WC,wcp);

%Calculates thermal conductivity based on Chung

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

[b1,b2,b3,bd_water,heatcapsolids_corr, heatcapwater_corr,disp]= In_temperature_param3(units,ncomp);

    %%calculate water content halfway the node, halfway the time step
    i = [2:ncomp+1];
    tw(1) = (WC(1) + wcp(1))/2;   
    tw(i) = (WC(i) + wcp(i) + WC(i-1) + wcp(i-1))/4;
    tw = tw';
    
%Calculate volumetric heat capacity (eq 4-2 Wave manual)
i = [1:ncomp+1];
beta(i) = heatcapsolids_corr*bd(i) + heatcapwater_corr*tw(i)*bd_water; %in J/C/L^3 
beta = beta';

Cw(i)= heatcapwater_corr*tw(i)*bd_water; %in J/C/L^3 
Cw = Cw';

lo_chung = b1 + b2*WC + b3*(tw).^0.5; 
allam = lo_chung;

%l_theta_chung = lo_chung + disp.*Cw.*abs(wat_flxs');

