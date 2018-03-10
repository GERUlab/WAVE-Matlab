%Penmann-Monteith implementation by Marnik Vanclooster, Jose Luis Gabriel
%Perez
function [ET0_time ,ET0_cm_per_day] =calc_ET02(climate)
%Read in climate data
time = jday2matlab(climate(:,2), climate(:,1))+ climate(:,3)/24;
Tair = climate(:,5);
RH = climate(:,6);
SlrMJm2 = climate(:,4)*60*60*10^-6;
u = climate(:,8);%Climatic file needs to be updated to also include wind speed. OK

[timeday u2] = convert_to_daily_mean(u,time);
[timeday tave] = convert_to_daily_mean(Tair,time);
[timeday rs] = convert_to_daily_sum(SlrMJm2,time);
[timeday RH] = convert_to_daily_mean(RH,time);

% Saturated vapour pressure, kPa 
esat=0.6108.*exp(17.27.*tave./(tave+237.5)); % Eqn.11 of Allen et al. (1998) 
% Actual vapour pressure, kPa 
eact =RH./100.*esat;
% vapor pressure deficit, kPa 
dsat = esat-eact;

[year month day hr minutes seconds]=datevec(timeday); % convert datenumber to year, month, day, ...
dnr=day-32+floor(275*month/9)+2*floor(3./(month+1))+floor(month/100-mod(year,4)/4+0.975); % Eq. 52; julian day number (1 .. 365 or 366)

%% Calculate reference evapotranspiration according FAO-56
% Reference: http://www.fao.org/docrep/X0490E/X0490E00.htm#Contents
% Parameters
ele= 150.; % elevation above sea level(m)
degreenord=50; %longitude belgium (°)
minutenord=50; %longitude belgium (minutes)

% % Incoming climate data 
% u2=data(:,5); %Windspeed, m/s
% tave=data(:,3); % Mean temperature °C
% dsat=data(:,7)/10; % Vapour pressure deficit, kPa
% rs=data(:,2).*1000./(3600*24);  %Incoming solar radiation, MJ/m2/day 
% 
% 
% % Actual vapour pressure, kPa
% eact=esat-dsat; 

%Extraterestrial radiation (ra), MJ/m2/day
gsc=0.082;
dec_deg=degreenord+(minutenord.*1.)./60.;	%Coordinates Nord Bruxelles (50° 50' N)								
fi=(dec_deg./180.).*pi; % Coordinates in radians
dr=1.+0.033*cos(2.*pi*dnr./365.); % Distance earth sun, Eqn 23 of Allen et al. (1998)
delta=0.409*sin(2.*pi.*dnr./365.-1.39); % Solar declination, Eqn 24 of Allen et al. (1998)
omega=acos(-tan(fi).*tan(delta)); %  Eqn 25 of Allen et al. (1998)
x1=24.*60./pi;
x2=sin(fi);
x3=sin(delta);
x4=cos(fi);
x5=cos(delta);
x6=sin(omega);
ra=max(x1.*gsc.*dr.*(omega.*x2.*x3+x4.*x5.*x6),0);	%  Eqn 21 of Allen et al. (1998)
clear x1 x2 x3 x4 x5 x6;

%Clear sky solar radiation, MJ/m2/day
rso=(0.75+2.E-5*ele)*ra; %  Eqn 37 of Allen et al. (1998) 

% Net short wave radiation (albedo = 0.23), MJ/m2/day
rns=0.77*rs;   % Eqn 38 of Allen et al. (1998)

% Long wave radiation, MJ/m2/day 
stef_bol=4.903E-9; 
x1=((tave+273.).^4)./2.;
x2=0.34-0.14.*sqrt(eact);
x3=1.35*(rs./rso)-0.35;
rnl=max(stef_bol.*x1.*x2.*x3,0);	
clear x1 x2 x3; 

%Net radiation MJ/m2/day
rn=rns-rnl;

% Slope of saturation deficit curve 
x1=tave+237.3; 
x2=17.27.*tave./x1;
p1=4096.;
p2=0.6108;
delta=p1.*p2.*exp(x2)./x1.^2;
clear x1 x2;

% Gamma 
x3=(293.-0.0065.*ele)./293.;
p=101.3.*x3.^5.26; % air pressure (kPa) 
labda=2.45; 
eps=0.622; 
cp=1.013E-3;
gamma= cp.*p./(eps.*labda);
clear x3; 

%et0
x5=0.408.*delta.*rn;
x6=gamma*900./(tave+273.);
x7=u2.*dsat;
x8=delta+gamma.*(1.+0.34.*u2);
et0=(x5+x6.*x7)./x8;
clear x5 x6 x7 x8;

ET0_time = timeday;
ET0_cm_per_day= et0/10;

% figure
% plot(et0(:,1),et0(:,2));
% xlabel('Temps (Jours)');
% ylabel('Evapotranspiration de réference (mm)');
% datetick('x',19)
% xlim([et0(1,1) et0(end,1)])
% 