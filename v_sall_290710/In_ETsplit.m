function  [AKc,epa,esa,drz,Smax_param, Smax_given] = In_ETsplit(t,simplant,units,plant_date,harvest_date);

% Calculation of the different plant water uptake parameters/variables
% IN:
% 	t: current time 
% OUT:
%AKc : crop coefficient
%epa: potential transpiration (mm)
%esa: potential evaporation mm)
%drz: rooth length
% Smax_param ()
    % pressure head below which the roots  start to extract water from the  soil: P0(cm)
    % pressure head below which the roots start to extract water optimally from the soil: P1(cm)
    % pressure head below which the roots can no c   longer extract water optimally
%                               at a high evaporative demand : P2h (cm)       
%                               at a low  evaporative demand: P2l (cm)
    %pressure head at which the water uptake by the roots ceases (wilting point)P3: (cm)
    % ilinr: switch for the relation between the reduction of the root wateruptake and pressure head (linar = true)
    %arer: intercept for water uptake (day-1)
    %brer: slope for water uptake (mm-1 day-1)
%Smax_given: switch for Smax input or Smax calculated (= 1 or 0)

%%CALL: In_Evaporation
%%CALLED BY: RER
%-------------------------------------
% F. Hupet 15/10/2000
%Modified by M. Sall 21/11/2009

time=[0:1:49]';
 Eto = [ones(50,1)*2.5]; %%Eto= In_evaporation(); %%% (mm)
drz= [ones(10,1)*200; ones(20,1)*300; ones(20,1)*400]; %%% (mm)
rlai= [ones(10, 1)*1.3; ones(20, 1)*2.5; ones(20, 1)*3.5];
kc=[ones(10,1)*0.85; ones(20,1)*1.05; ones(20,1)*0.6];


PLANT_WATER_UPTAKE_MATRIX=[time Eto drz rlai kc]; 
    
%% parametres for water uptake reducion
p0= -10;             %%% (cm)
p1= -25;             %%% (cm)
p2h= -800;             %%% (cm)
p2l= -1500;             %%% (cm)
p3= -8000;              %%% (cm)
Smax_given=1;
Smax=[ones(20,1)*0.031; ones(20,1)*0.0000; ones(20,1)*0.00; ones(10,1)*0.;...
    ones(10,1)*0.0; zeros(60,1)]; 
ilinr = 1;   %% used if Smax_given=0;
arer= 0.035;  %%  used if Smax_given=0;
brer= 0.00007;  %% used if Smax_given=0;


% Searches for the appropriate plant water uptake for the current time t 
i=max(find(PLANT_WATER_UPTAKE_MATRIX(:,1)<=t));

ETo=PLANT_WATER_UPTAKE_MATRIX(i,2);%en mm jour-1
drz=PLANT_WATER_UPTAKE_MATRIX(i,3); %%% mm
RLAI=PLANT_WATER_UPTAKE_MATRIX(i,4); %%%m-2 m-2
AKc=PLANT_WATER_UPTAKE_MATRIX(i,5);
%ARER=PLANT_WATER_UPTAKE_MATRIX(i,6); %%%day-1
%BRER=PLANT_WATER_UPTAKE_MATRIX(i,7); %%%mm-1 day-1

if simplant==0 |t<plant_date |t>=harvest_date; %% NO CROP
   AKc=1;
   esa=ETo;
   epa=0;
   RLAI=0;
   drz=0;
   %ARER=0;
   %BRER=0;
elseif simplant & t>=plant_date & t<harvest_date;    %Simulation with crop
   EToA=AKc*ETo;	%Potential evapotranspiration of the crop (mm day-1)
   esa=exp(-0.6*RLAI)*EToA;		%Potential evaporation of the soil (mm day-1)
   epa=EToA-esa;	%Potential transpiration of the crop (mm day-1)
   esa=esa;  %%% mm day-1
   epa=epa;  %%% mm day-1
end

          %% (= 1 if Smax is input else Smax_given=0)

%%Conversion of esa and epa to the units addopted for the simulation
%geometry conversion
if units{1,1}=='mm';
    esa=esa; epa=epa; drz=drz; brer=brer;
elseif units{1,1}=='cm';
    esa=0.1*esa; epa=0.1*epa; drz=0.1*drz; brer=10*brer;
elseif  units{1,1}== 'm';
    esa=0.001*esa; epa=0.001*epa; drz=0.001*drz; brer=1000*brer;
else
    error('the geometry unit must be m; cm or mm (check file In_general_data)')
end

%time conversion

if units{1,2}=='day';
    esa=esa; epa=epa; arer=arer; brer=brer;
elseif units{1,2}=='h';
    esa=esa/24; epa=epa/24; arer=arer/24;brer=brer/24;   
elseif units{1,2}=='min';                       %% units{1,1}= m
    esa=esa/1440; epa=epa/1440; arer=arer/1440;  brer=brer/1440;     
elseif units{1,2}=='s'; 
    esa=esa/86400; epa=epa/86400; arer=arer/86400; brer=brer/86400; 
else 
    error('the time unit must be day; h;; min or s, (check file In_general_data)')
end 

if Smax_given
    Smax_param=[p0; p1; p2l; p2h; p3; ilinr; Smax];
else
    Smax_param = [p0 p1 p2l p2h p3 ilinr arer brer ];
end

% EVAPORATION OF INTERCEPTED WATER (in WAVE.FOR (RER.FOR))
%ev_intc=min(stor,epa)
%epa_min_intc= epa-av_intc;
%stor=stor-ev_intc;



