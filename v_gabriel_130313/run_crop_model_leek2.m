function [crop_result, RootResDay,crop_result2] =... 
                 run_crop_model_leek2(LeekParams,output_crop,simplant, plant_date, harvest_date,t, Nuptake_time,output_cropNup,PLM2,cropinit,crop_climate)

ncrop = find_cropnumber (t,simplant,plant_date, harvest_date); 
JDAY = t-plant_date(ncrop)+1;
	DATE = t;%mod(TIME+NSTART,365);
	DAP = JDAY;
NDAYS = harvest_date(ncrop) - plant_date(ncrop);
NSTART = plant_date(ncrop);

%Read in climate file for simulation
RefClim = crop_climate; %dlmread('last_used_climate.txt');

%%
%Loop parameters
    %per day simulation
    DELT = 1;
    %Per hour simulation
    NFAST = 24;
    DTFAST = 1.0/NFAST;

%development stage parameters
    %base temperature, blow which plant growth ceases
    Tb = LeekParams(10);

%Photosynthesis parameters %from Baumann 0.44-0.8
    XK = LeekParams(11);
    %universal ones:
    XM = LeekParams(12);
    QE = LeekParams(13);
    TAU1 = LeekParams(14); 
    TAU2 = LeekParams(15);

%Temperature effect on photosynthesis PMAX
PGRED = LeekParams(16:24);
TMPG = LeekParams(25:33);

%DVS effect on photosynthesis PMAX and SLA
 PDVS = LeekParams(34:41);
 
 %SLA in cm² per g
 DAPN = LeekParams(50:57);
 DVSArray = LeekParams(42:49);DVSArray(end) = DVSArray(end)+600;
 
 SLA = LeekParams(58:65);
 
%  plot(DAPN(:),SLA(:)); xlabel('TSum'); ylabel('SLA');
%  plot(DVSArray(:),SLA(:)); xlabel('TSum'); ylabel('SLA');
 
%%
%This part needs to be changed into N-supply in the soil
    
    %N supply effect on photosynthesis
dosis = LeekParams(66:70);
NCor =  LeekParams(71:75);

% plot(dosis(:),NCor(:)); xlabel('Dosis N'); ylabel('Correction factor for N dosis on PMax');

FertDosis = LeekParams(4);

    %N-content in leaf and shaft for dosis=3 and methode=4, calculated on 09 data
DAPN = LeekParams(50:57);
PercNShaft = [3.03000000000000 3.03000000000000 2.84000000000000 2.50000000000000 2.44000000000000 2.35000000000000 2.28000000000000 2.28000000000000];
%LeekParams(76:83);
PercNLeaf = LeekParams(84:91);
TSumN = LeekParams(92:99);

TSumN = LeekParams(92:99);

        %par(mfrow=c(2,1))
        %plot(TSumN,PercNLeaf,type='l')
%         %plot(TSumN,PercNShaft,type='l')
%         figure; 
%         subplot(2,1,1);plot(TSumN(:),PercNLeaf(:));
%         subplot(2,1,1);plot(TSumN(:),PercNShaft(:));

%VPD effect on photosynthesis
VPDL = LeekParams(100);
CK = LeekParams(101);

%Plant density in pl/m^2
%PLM2 = LeekParams(5);

%Proportionality factor for coverage calculated on PLM2=18.18 to PLM2 = 15.58
PropFact = LeekParams(6);

%%
if DAP ==1
%Initialisation parameters, 
    %Initial DW leaves in g
    DMLeaf = cropinit.DWLeaf(ncrop);%LeekParams(7);
    %Initial leaf area in m^2
    TotLeafArea = DMLeaf * cropinit.SLA(ncrop)/10000;
    %Initial DW stem in g
    DMStem = cropinit.DWStem(ncrop);
    %Initial nitrogen content
    InitN(1) = cropinit.NLeaf(ncrop) ;%
    InitN(2) = cropinit.NStem(ncrop) ;%
end
%Maintenance respiration parameters
    %Temperature effect parameter
    Q10 = LeekParams(102);
    %Respiration parameters for Leaf and Curd (F)g(CH2O)/g(DW) day
    RMRL = LeekParams(103);
    RMRF = LeekParams(104);

%Growth efficiency in g (DM) / g (CH2O)
    GREF = LeekParams(105);

% calibration frac leaf and shaft 2009 data in relation to T Sum
DVSArray2 = LeekParams(106:113);DVSArray2(end) = DVSArray2(end)+600;
LeafFrac = LeekParams(114:121);
StemFrac = LeekParams(122:129);
RootFrac = LeekParams(130);

        %win.graph()
        %par(mfrow=c(2,1))
        %plot(DVSArray2,LeafFrac,type='l',xlab='T Sum',ylab='Leaf Fraction')
        %plot(DVSArray2,StemFrac,type='l',xlab='T Sum',ylab='Stem Fraction')
        %LeafFrac + StemFrac
%         figure;
%         plot(DVSArray2(:), LeafFrac(:)); xlabel('T Sum'); ylabel('Leaf Fraction');
%         plot(DVSArray2(:), StemFrac(:)); xlabel('T Sum'); ylabel('Stem Fraction');

%Deadleaf fraction in relation to Effective Temperature (Te)
%#Tef = c(0,5,10,15,20,25,30,35,40)
%#Kd = c(0.1,0.2,0.3,0.4,0.5,0.6,0.6,0.6,0.6)
%#par(mfrow=c(1,1))
%#plot(Tef,Kd,type='l',xlab='Te',ylab='Kd dead leaves')

%Leaf coverage in relation to TSum  data 2009
DAPCov = LeekParams(131:139);
DVSAr = LeekParams(140:148);DVSAr(end) = DVSAr(end)+600;
LeafCoverage = LeekParams(149:157);

        %win.graph()
        %plot(DAPCov,LeafCoverage,type='l',xlab='DAP',ylab='Leave Coverage in %')
        %#plot(DVSAr,LeafCoverage,type='l',xlab='T Sum',ylab='Leave Coverage in %')
%         figure;
%         plot(DAPCov(:), LeafCoverage(:)); xlabel('DAP'); ylabel('Leaf Coverage in %');
        
%Leaf coverage in relation to TSum  Preliminary extraction fromdata 2011
%#DVSAr<-c(0,358,383,775,1170,1533,1834,2045,2114,3000)
%#Leaf coverage in percent
%#LeafCoverage<-c(1,4,10,22,41,52,59,64,64,64)
%#plot(DVSAr,LeafCoverage,type='l',xlab='DVS',ylab='Leave Coverage in %')

%Initial LAI
if DAP ==1
LAICov = TotLeafArea*PLM2 / (LeafCoverage(1)/100);
end
%par(mfrow=c(1,1))

%Root DM partitioning function for 11 virtual layers
RDPF = zeros(1, 11);%LeekParams(158:168);

%Virtual root depth
VirRD = zeros(1,11);%LeekParams(169:179);
        %plot(RDPF,-VirRD,type='l',ylab='Virtual Rooting Depth',xlab='Root DM partitioning')
%        figure; plot(RDPF(:),-VirRD(:)); xlabel('Root DM partitioning'), ylabel('Virtual Rooting Depth'); 

%Initialise results file
if DAP==1
results = zeros(1, 34);


%Initialisation of Effective Temperature Sum (TeSum)
TeSum = 0;

%Initialisation of integration
TotDMLeaf = DMLeaf;
TotDMStem = DMStem;
TotDMRoot = 0;
TotDMShoot = DMLeaf + DMStem;
TSum = 0;
TotNperLeaf = 0;
TotNperStem = 0;
CumTranspiration = 0;
CumEvaporation = 0;
TotNperLeaf = 0;
TotNperShaft = 0;
CumTranspiration = 0;
CumEvaporation = 0;
TotNDemandLeaf = 0;
TotNDemandStem = 0;
TotNDemandLeafStem = 0;
TotNDemandRoot = 0;
CumPlantNUptakeActualLeaf = TotDMLeaf*InitN(1)/100 ;
CumPlantNUptakeActualStem =  TotDMStem*InitN(2)/100;
CumPlantNUptakeActualLeafStem = CumPlantNUptakeActualLeaf + CumPlantNUptakeActualStem;
CumPlantNUptakeActualRoot = 0 ;
%Initiatlisation
Holding(1,:) = [0,0,0,0,0,0,0,0,0];

%Initial root depth at DAP=0 in cm
TotRootDepth2 = LeekParams(9);
else
% DATE = output_crop{ncrop}(end,1);
% DAP = output_crop{ncrop}(end,2)
GP = output_crop{ncrop}(end,3);
RMAINT = output_crop{ncrop}(end,4);
MeanDayTemp= output_crop{ncrop}(end,5);
Te= output_crop{ncrop}(end,6);
TeSum= output_crop{ncrop}(end,7);
TotDM= output_crop{ncrop}(end,8);
DMRoot= output_crop{ncrop}(end,9);
DMLeaf= output_crop{ncrop}(end,10);
LeafArea= output_crop{ncrop}(end,11);
DMStem= output_crop{ncrop}(end,12);
TotDMLeaf= output_crop{ncrop}(end,13);
TotDMStem= output_crop{ncrop}(end,14);
TotDMRoot= output_crop{ncrop}(end,15);
TotDMShoot= output_crop{ncrop}(end,16);
TotLeafArea= output_crop{ncrop}(end,17);
LAICov= output_crop{ncrop}(end,18);
LAI= output_crop{ncrop}(end,19);
NLeaf= output_crop{ncrop}(end,20);
NStem= output_crop{ncrop}(end,21);
NperLeaf= output_crop{ncrop}(end,22);
NperStem= output_crop{ncrop}(end,23);
TotNperLeaf = output_crop{ncrop}(end,24);
TotNperStem = output_crop{ncrop}(end,25);
TotRootDepth1 = output_crop{ncrop}(end,26);
RCil = output_crop{ncrop}(end,27);
Transpiration= output_crop{ncrop}(end,28);
Evaporation= output_crop{ncrop}(end,29);
CumTranspiration= output_crop{ncrop}(end,30);
CumEvaporation= output_crop{ncrop}(end,31);
TotNperShoot = output_crop{ncrop}(end,32);
ShaftDi =output_crop{ncrop}(end,33);
ShaftLe = output_crop{ncrop}(end,34);
TotNDemandLeaf = 0;

CumPlantNUptakeActualLeaf =output_cropNup{ncrop}(end,1) ;
CumPlantNUptakeActualStem =  output_cropNup{ncrop}(end,2);
CumPlantNUptakeActualLeafStem = output_cropNup{ncrop}(end,3) ;
CumPlantNUptakeActualRoot = output_cropNup{ncrop}(end,4) ;
end
% PercNShoot = CumPlantNUptakeActualLeafStem/(TotDMLeaf+TotDMStem)*100;
% PercNLeaftab = CumPlantNUptakeActualLeaf/(TotDMLeaf)*100;
% PercNStemtab = CumPlantNUptakeActualStem/(TotDMStem)*100;
% NFac = PercNLeaftab/ interp1(TSumN,PercNLeaf,TeSum);
% NFac = min(NFac,1);
% NFac = max(NFac, 0.05);

%**********************************************************************
%Daily loop
%**********************************************************************
% for JDAY = 1:NDAYS
% 	TIME = (JDAY-1)*DELT;
% 	DATE = mod(TIME+NSTART,365);
% 	DAP = JDAY;
    %Initialize variables that accumulate during the day
      GP = 0;
      RMAINT = 0;
      Transpiration = 0;
      Evaporation = 0;
    %Parameter to calculate average day temperature
        Tsum = 0;
     
     %*********************************************************************
     %Fast loop (if time-step = 1 hour, NFAST = 24)
     %*********************************************************************
    for JF = 1 : NFAST
        
        % TFAST is the hour of the day (h)
        TFAST = (JF-1)* 24/NFAST;
        
        %*********************************************************************
        %Input of hourly climate data --- source(RefClim) ---
        %*********************************************************************
        d = find(RefClim(:,1)==DATE & RefClim(:,2)==(TFAST+1));
            if length(d)>1 
                disp('more than 1 found!!');
            end
            if isempty(d)
                disp(' invalid index of climate data');
            end
		JUL = RefClim(d,1);
		IHE = RefClim(d,2);
		solrad = RefClim(d,3);
		TMPA = RefClim(d,4);
		RH = RefClim(d,5);
		%VPD calculation in kPa
		Ps = 0.61078 * exp((17.2694*TMPA)/(TMPA+237.3));
		%Calculate actual vapor pressure PV
		Pv = RH/100 * Ps;
		VPD = Ps-Pv;
     VPD_nok = VPD*(10^3);
		CO2L = LeekParams(181);
		PAR = solrad * LeekParams(182); 
	    PPFD = PAR * LeekParams(183);
        
            %*********************************************************************
            %Calculation of gross photosynthesis --- source('PHOTO2.R') Acock's model---
            %*********************************************************************
            GPF = 0;
			%Effect of CO2 on Pmax (from Gainesville)
			PMAX = TAU1 * CO2L;
                if CO2L>1500
                    PMAX = (TAU1 * 1500) + (TAU2 * (CO2L-1500));
                end		
			%Reduction of Pmax at extreme temperatures
			PMAX = PMAX * interp1(TMPG,PGRED,TMPA);

			%Reduction of Pmax in relation to DVS, senescence  in LEEK ?????
			PMAX = PMAX * interp1(DVSArray,PDVS,TeSum);
            
           
			%Reduction of PMax in relation to N-dosis
%			PMAX = PMAX * interp1(dosis,NCor,FertDosis);
% if isempty (output_cropNup{1,ncrop})==0
% PMAX = PMAX * output_cropNup{1,ncrop}(end,9);
% end	
% 
if DAP < 30
    PMAX = PMAX/2;
else
    PMAX = PMAX*1.2;
end
			%Reduction of XK with leaf area in m^2(from Baumaan  0.8 to 0.44
			%XK<-ifelse(TotLeafArea < 0.235,0.90 - 3.00*TotLeafArea + 6.38*TotLeafArea^2,XKK)

			if PPFD >= 0.001
			    % reduction of Pmax by VPD
			    if VPD>=VPDL 
                    PMAX = PMAX * exp(CK*(VPD-VPDL));
                end
                % Acock's model
			    TOP = (1-XM) * PMAX + QE * XK * PPFD;
			    BOT = (1-XM) * PMAX + QE * XK * PPFD * exp(-XK * TotLeafArea * PLM2);

			    GPF = (PMAX/XK) * log(TOP/BOT);
  
			    % conversion from CO2 to CH2O (30/44 = 0.682)
			    GPF = GPF * 0.682;
			    % conversion of GPF from µM/m2-s into g/m2-day
                    %µM/m2-s x 0.000044g/µM x 3600s/h x 24h/d = 3.8016 g/m2-day
			    % GPF in g per day per plant
			    GPF = GPF * 3.8016/PLM2;
            end

            %*****************************************************************
            %Calculation of maintenance respiration 
            %*****************************************************************
        	% Effect of temperature on maintenance respiration on hourly basis
     		RMAINTF = (RMRL * (DMLeaf)+RMRF * DMStem) * (Q10^(0.1 * TMPA-2.0));

            %integration of variables on 24 hours
            GP = GP + GPF * DTFAST;
            RMAINT = RMAINT + RMAINTF * DTFAST;
    
    %*****************************************************************
    %Crop transpiration on hourly basis 
    %*****************************************************************
    %Internal resistance in s/m
    ri = LeekParams(184);
    re = LeekParams(185);
    %Specific heat of air J/kg K
    Cp = LeekParams(186);
    %AtmPressure in Pa
    AtmP = LeekParams(187);

    %Short wave radiation calculations W/m^2
    Is = solrad * LeekParams(188);

    %Long wave rad calculations  W/m^2
    Il = solrad * LeekParams(189);

    %Air density in kg/m^3
    rho.air = 100000/(287*(TMPA+273.16));

    %Slope of saturationn curve Pa/C
    delta = 41.45 * exp(0.06088*TMPA);

    %Heat of evaporisation J/kg
    lambda = 2502535.259-(2385.76*TMPA);

    %Psychrometric constant along to Moneith p 181 Value=67 Pa/C at 20 C
    gamma = (Cp * AtmP)/(lambda * 0.622);

    %Canopy temperature day and nigth
    %T.Canopy = ifelse(solrad>5,TMPA+0.006*Is-0.25*(VPD/gamma),TMPA-0.1*(VPD/gamma))
 
    %Extended ri and re calculations, Stanghelini...VERY UNSTABLE
    %l std leaf length in m
    %l=0.3 
    %u is windspeed in m/s
    %u=0
    %re boundary layer resistance on canopy level
    %re=(1174*l^0.5)/((l*abs(T.Canopy-TMPA) + 207*u^2)^0.25)
    %ri internal resistance;
    %ri=82*((Is+4.3)/(Is +0.54))*(1 + 0.023*(T.Canopy-24.5)^2)*((0.004 + exp(-0.54*VPD*1000))^-0.25)

    %Extended Rn calculations
    x = exp(-0.24 * LAICov);
    %As absorption coefficient short wave radiation, smoothing rules Stanghellini
    As = 0.83 * (1-0.70 * x) * (1+0.58 * x * x) * (0.88-x * x + 0.12 * x^(8/3)); 
    %Al absorption coefficient long wave radiation
    Al = 0.83 * (1-0.70 * x) * (1-x^(8/3));
    Rn = As * Is + Al * Il;

    %H = convection loss by canopy if not neglected
    %H=(LAICov*rho.air*Cp*(TMPA-T.Canopy))/re
    H = 0;
 
    %Transpiration heat in W/m^2 covered ground surface
    LE = ((1/(1+(delta/gamma) + (ri/re))*( ((delta/gamma)* (Rn-H)) + (rho.air * Cp *  VPD_nok)/(gamma * re))))*TotLeafArea * PLM2;
    %Transpiration in kg/s m^2 coverage
    E = LE/2/lambda;
    EHour = E * 3600;

        %integration of variables on 24 hours
        Transpiration = Transpiration + EHour;


    %*****************************************************************
    %Soil evaporation. Classical Penmann equation (Monteith, Unsworth, 1990) with resistance 
    %*****************************************************************
    %Surface conductance m/s
    rs = 1/300;%LeekParams(190);
    albedo = LeekParams(191);
    %ESoil in kg / m^2 s
    ESoil = (delta * solrad * albedo + rho.air * Cp *  VPD_nok * rs) / (lambda * (delta + gamma));
    ESoilHour = ESoil * 3600 ;%* (1-(interp1(DVSAr,LeafCoverage,TeSum)/ 100));

    %Integration of variables on 24 hours
        Evaporation = Evaporation + ESoilHour;

    %*****************************************************************
    %Calculate temperature sum for average day temperature calculations
    %*****************************************************************

    %calculate average day temperature
	Tsum = Tsum + TMPA;

    end
    %*****************************************************************
    %End of fast loop
    %*****************************************************************

%GP and RMAINT in g CH2O per plant
    %Transfert to Biomass, calculate Total DM in g per plant per dag
TotDM = (GP - RMAINT) * GREF; 

%Calculate T Sum  as measure for development stage
	MeanDayTemp = Tsum/24;
	Td = MeanDayTemp - Tb;
%Effective temperature set to zero if negative, non-decreasing Tsum
    if Td > 0
        Te = Td;
    else Te = 0;
    end
%Calculate Temperature sum
	TeSum = TeSum + Te;

%Compartimentation over leaf, stem and curd
    %Growth rates g per planr per day
DMLeaf = TotDM * interp1(DVSArray2,LeafFrac,max(TeSum, DVSArray2(end)))/(1 + RootFrac);
DMStem = TotDM * interp1(DVSArray2,StemFrac,max(TeSum, DVSArray2(end)))/(1 + RootFrac);

%Compartimentation to roots, 9 % of shoot and constant
DMShoot = DMLeaf + DMStem; 
DMRoot = DMShoot * RootFrac;

%Calculate leaf area out of SLA-DVS curve
LeafArea = DMLeaf * (interp1(DVSArray,SLA,TeSum)/10000);

%Integrate to cumulative total DM
TotDMLeaf = TotDMLeaf + DMLeaf;
TotDMStem = TotDMStem + DMStem;
TotDMShoot = TotDMShoot + DMShoot;
TotDMRoot = TotDMRoot + DMRoot;
TotLeafArea = TotLeafArea + LeafArea;


%Diameter (mm) en length (cm) calculations
RGRDi = 0.0006775542;   % mm / d C
RGRLe = 0.0004561056;  % cm / d C 
ShaftDi = exp(2.3095946463) * exp(RGRDi * TeSum);
ShaftLe = exp(2.0991201269) * exp(RGRLe * TeSum);


%Direct LAI calculation 
LAI = TotLeafArea * PLM2;
%direct LAI calculation per m^2 coverage
LAICov = TotLeafArea * PLM2 / (interp1(DVSAr,LeafCoverage,TeSum)/ 100);

%Integrate transpiration and evaporation to cumulative values
CumTranspiration = CumTranspiration + Transpiration;
CumEvaporation = CumEvaporation + Evaporation;


    %*****************************************************************
    %Root volume and length calculations
    %*****************************************************************
    %Specific root length in m/g
    SRL = LeekParams(180);
    TLayer =1;% LeekParams(192);
    %RootDepth in cm follows  linear curve in function of TeSum
    TotRootDepth1 = (10 + 0.02088 * TeSum)*1;
%     rootDAP = [0 32 68 102 200];
% z1_param = [7.45 7.45 3.05 2.63 2.63]; 
% z2_param = [1.96 1.96 4.99 5.34 5.34];
% rootdensparam1 =interp1(rootDAP,z1_param,DAP); %Find current parameter1
% rootdensparam2 =interp1(rootDAP,z2_param,DAP); %FInd current parameter2
% depth = 1:TLayer:100;
% rdens = gampdf(depth,rootdensparam1,rootdensparam2); %Calculate rdens over profile
% maxindex = find(rdens == max(rdens)); 
% minindex = find(rdens< 0.0001); %Calculate depths with low root densities
% TotRootDepth1 = minindex(min(find((minindex>maxindex)==1))); %find the firs

    %Calculate number of layers, layer tickness = TLayer in cm

    NumLayer = round(TotRootDepth1/TLayer);
    %Calculate DM per layer, partition along Root DM Partitioning function 
    RDMLa(1:NumLayer) = 0;
    %Root Length per layer
    RLengthLa(1:NumLayer) = 0;
    %Root Volume per layer, using big root model with radius rr in cm
    rr = LeekParams(193);
    RVolumeLa(1:NumLayer) = 0;
    %Equivalent radius of the volume with tickness TLayer, in the most simple model this is a cilinder per layer
    ERLa(1:NumLayer) = 0;
    %Root Length Density per layer in cm/cm^3
    RLDLa(1:NumLayer) = 0;

    %Initialize RootRes
    RootRes = zeros(NumLayer,6);

    %Loop over layers to calculte sum of RDPF
    RDDistProp(1:NumLayer) = 0;
    RDDistSum = 0;
    for j = 1 : NumLayer
        RDDist = VirRD * NumLayer;
        RDDistProp(j) =0;% interp1(RDDist,RDPF,j);
        RDDistSum = RDDistSum + RDDistProp(j);
    end

    %Loop over layers
    for i = 1 : NumLayer
        RDDist = VirRD * NumLayer;
        RDMLa(i) =  0;%TotDMRoot * (interp1(RDDist,RDPF,i)/RDDistSum);
        RLengthLa(i) = 0;%RDMLa(i) * SRL;
        RVolumeLa(i) = 0;%RLengthLa(i) * pi * rr^2;
        ERLa(i) = 0;%sqrt(RVolumeLa(i)/(TLayer*pi));
        %For now we use the actual volume to calculate Root Length Density, could also use cilider with RCil as radius and TLayer as tickness
        if RVolumeLa(i) ~= 0
            RLDLa(i) =0;% RLengthLa(i);
        else RLDLa(i) = 0;
        end
        RootRes(i,:) = [i,RDMLa(i),RLengthLa(i),RVolumeLa(i),ERLa(i),RLDLa(i)];
    end


    %Calculate Radius of enclosing cilinder RCil
        %Root Length Density per layer in cm/cm^3 for enclosing cilinder
    RCil = max(ERLa);
    RLDRCil = RLengthLa;%RLengthLa /((RCil^2) * pi * TLayer);

    RDATE(1:NumLayer) = DATE; RRDATE = RDATE';
    RDAP(1:NumLayer) = DAP; RRDAP = RDAP';
    RLDRCil =RLDRCil';
    RootResDay = [RRDATE,RRDAP,RootRes,RLDRCil];

%    Holding = [Holding;RootResDay];

%*****************************************************************
%Calculate % N in the leaves, stems and curds, along to observed concentrations
%*****************************************************************

%N in leaves
NLeafOpt = interp1(DAPN,PercNLeaf,DAP);
NStemOpt = interp1(DAPN,PercNShaft,DAP);

% %Calculate grams N per plant
NperLeaf = (NLeafOpt/100) * DMLeaf;
NperStem = (NStemOpt/100) * DMStem;

%Integrate Tot N
TotNperLeaf = TotNperLeaf + NperLeaf;
TotNperStem = TotNperStem + NperStem;
TotNperShoot = TotNperLeaf + TotNperStem;


demandoption = 1;
if demandoption ==1
%N demand for optimal growth in leaves in g per plant
LeafNDemand = NLeafOpt * (TotDMLeaf)/100- CumPlantNUptakeActualLeaf;
LeafNDemand = max(0,LeafNDemand);

StemNDemand = NStemOpt * (TotDMStem)/100- CumPlantNUptakeActualStem;
StemNDemand = max(0,StemNDemand);
LeafStemNDemand = StemNDemand + LeafNDemand;
RootNDemand = 0;
elseif demandoption ==2
LeafNDemand = NperLeaf; 
StemNDemand = NperStem ;
LeafStemNDemand = StemNDemand + LeafNDemand;
RootNDemand = 0;
else
    if ncrop == 1
%forcingfunction = xlsread('ForcingNuptakeLeek.xls');
%dlmwrite('forcing_leek.txt', forcingfunction)
forcingfunction = dlmread('forcing_leek.txt');
dailydemand = diff(forcingfunction(:,3));
    else
%forcingfunction = xlsread('ForcingNuptakeLeek2.xls');
%dlmwrite('forcing_leek2.txt', forcingfunction)
forcingfunction = dlmread('forcing_leek2.txt');
dailydemand = diff(forcingfunction(:,3));
    end
demand = dailydemand(DAP);
StemNDemand = demand *0.5;
LeafNDemand = demand*0.4;
RootNDemand = demand*0.1;
LeafStemNDemand = LeafNDemand + StemNDemand;
end

%Calculate tota N demand per day in g per plant
PlantNDemand = LeafNDemand + StemNDemand+ RootNDemand;
if PlantNDemand ~=0
ratioleaf = LeafNDemand/PlantNDemand;
ratiostem = StemNDemand/PlantNDemand;
ratioroot = RootNDemand/PlantNDemand;
else
ratioleaf = 0;
ratiostem = 0;
ratioroot = 0;
    
end

%Integrate Tot N demand  
% TotNDemandLeaf = TotNDemandLeaf + LeafNDemand;
% TotNDemandStem = TotNDemandStem + StemNDemand;
% TotNDemandRoot = TotNDemandRoot + RootNDemand;

%TotNDemandPlant = TotNDemandLeaf + TotNDemandStem + TotNDemandRoot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NSupply van deze dag
if t== plant_date(ncrop)
    NSupply = PlantNDemand;
else
   NSupply = Nuptake_time(end,2); %Input of soil model
%   NSupply = PlantNDemand;
end
%with Npool add actualNpool to Nsupply

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PlantNDemand <= NSupply
%PlantNUptakeActual = PlantNDemand;
PlantNUptakeActual = NSupply;
else
PlantNUptakeActual = NSupply;
end

PlantNUptakeActualLeaf = PlantNUptakeActual*ratioleaf;
PlantNUptakeActualStem = PlantNUptakeActual*ratiostem;
PlantNUptakeActualLeafStem = PlantNUptakeActualLeaf + PlantNUptakeActualStem;
PlantNUptakeActualRoot = PlantNUptakeActual*ratioroot;


CumPlantNUptakeActualLeaf = CumPlantNUptakeActualLeaf+PlantNUptakeActualLeaf ;
CumPlantNUptakeActualStem = CumPlantNUptakeActualStem+PlantNUptakeActualStem ;
CumPlantNUptakeActualLeafStem = CumPlantNUptakeActualLeafStem+PlantNUptakeActualLeafStem  ;
CumPlantNUptakeActualRoot = CumPlantNUptakeActualRoot +PlantNUptakeActualRoot ;

%Two approaches : nitrogen concentration in the freshly formed organs
%(ActPlantNperc) of nitrogen concentration averaged over the shoot
%(PercNShoot)
%ActPlantNPerc = (PlantNUptakeActual/TotDM)*100; % Nitrogen concentration in newly formed organs
PercNShoot = CumPlantNUptakeActualLeafStem/(TotDMLeaf+TotDMStem)*100;
PercNLeaftab = CumPlantNUptakeActualLeaf/(TotDMLeaf)*100;
PercNStemtab = CumPlantNUptakeActualStem/(TotDMStem)*100;
%Calculate NFac. N correction factor on Pmax in PS model 
%NFac = 1 - (interp1(DVS_meas,Shoot_Nc,DVS) - ActPlantNPerc) / (interp1(DVS_meas,Shoot_Nc,DVS) - interp1(DVS_meas,Shoot_Nm,DVS));
%NFac = 1;
%NFac = 1 - (interp1(DVS_meas,Shoot_Nc,DVS) - PercNShoot) / (interp1(DVS_meas,Shoot_Nc,DVS) - interp1(DVS_meas,Shoot_Nm,DVS));
NFac = PercNLeaftab/ interp1(DAPN,PercNLeaf,DAP)*PercNStemtab/ interp1(DAPN,PercNShaft,DAP);
NFac = min(NFac,1);
NFac = max(NFac, 0.05);

crop_result = [DATE,DAP,GP,RMAINT,MeanDayTemp,Te,TeSum,TotDM,DMRoot,DMLeaf,...
    LeafArea,DMStem,TotDMLeaf,TotDMStem,TotDMRoot,TotDMShoot,TotLeafArea,LAICov,LAI,NLeafOpt,...
    NStemOpt,NperLeaf,NperStem,TotNperLeaf,TotNperStem,TotRootDepth1,RCil,Transpiration,Evaporation,CumTranspiration,...
    CumEvaporation,TotNperShoot,ShaftDi,ShaftLe];
crop_result2 = [CumPlantNUptakeActualLeaf,CumPlantNUptakeActualStem,CumPlantNUptakeActualLeafStem, CumPlantNUptakeActualRoot,...
PlantNUptakeActualStem,PlantNUptakeActualLeaf, PlantNUptakeActualLeafStem,...
    PlantNUptakeActualRoot,NFac,PercNShoot,PercNLeaftab, PercNStemtab,PlantNDemand];

end
