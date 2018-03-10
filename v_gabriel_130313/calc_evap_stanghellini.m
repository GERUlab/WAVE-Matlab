function ESoilDay = calc_evap_stanghellini(t)
%Calculate evaporation according to Stangehellini (in case that no crop is
%present on the field) ==> needed as input into the soilmodel.

%Read in climate file for simulation
RefClim = dlmread('last_used_climate.txt');

%Find relevant climatic data for that specific day
ESoilDay =0;
for i=1:24
d = find(RefClim(:,1)==round(t)+1 & RefClim(:,2)==(i));

JUL = RefClim(d,1);
IHE = RefClim(d,2);
solrad = RefClim(d,3);
TMPA = RefClim(d,4);
RH = RefClim(d,5);

%VPD calculation in kPa
Ps = 0.61078 * exp((17.2694 * TMPA)/(TMPA+237.3));
%Calculate actual vapor pressure PV
Pv = RH/100 * Ps;
VPD = Ps-Pv;
VPD_nok = VPD*10^3;
%Specific heat of air J/kg K
Cp = 1000;   
%Heat of evaporisation J/kg
lambda = 2502535.259-(2385.76*TMPA);

%AtmPressure in Pa
AtmP = 100000;

%Psychrometric constant along to Moneith p 181 Value=67 Pa/C at 20 C
gamma = (Cp*AtmP)/(lambda*0.622);

%Slope of saturationn curve Pa/C
delta = 41.45*exp(0.06088*TMPA);    

%Air density in kg/m^3
rho.air = 100000/(287*(TMPA+273.16));
%--------------------------------------------------------------------------
%Soil evaporation. Classical Penmann equation (Monteith, Unsworth, 1990) with resistance 
%--------------------------------------------------------------------------
    %Surface conductance m/s
    rs = 1/70;
    albedo = 0.13;
    %ESoil in kg / m^2 s
    ESoil = (delta*solrad*albedo + rho.air*Cp*VPD_nok*rs) / (lambda*(delta + gamma));
    ESoilHour = ESoil*3600;%(1-0/100);
    
    ESoilDay = ESoilDay + ESoilHour;
end