function         [crop_result, RootResDay,crop_result2] =... 
                 run_crop_modeles2(CropParams,output_crop,simplant, plant_date, harvest_date,t, Nuptake_time,output_cropNup)
ncrop = find_cropnumber (t,simplant,plant_date, harvest_date); 
JDAY = t-plant_date(ncrop)+1;
	DATE = t;%mod(TIME+NSTART,365);
	DAP = JDAY;
   
NDAYS = harvest_date(ncrop) - plant_date(ncrop);
NSTART = plant_date(ncrop);

%Read in climate file for simulation
RefClim = dlmread('last_used_climate.txt');

%Planting date
NSTART = CropParams(2);
%From plant to end of harvest
NDAYS = CropParams(3);

%Loop parameters
    %per day simulation
    DELT = 1;
    %Per hour simulation
    NFAST = 24;
    DTFAST = 1.0/NFAST;
	TIME = (JDAY-1) * DELT;
%development stage parameters
    %base temperature, blow which plant growth ceases
    Tb = CropParams(11);


%Photosynthesis parameters %from Olesen, 1997
    XKK = CropParams(14);
    %universal ones
    XM = CropParams(15);
    QE = CropParams(16);
    TAU1 = CropParams(17);
    TAU2 = CropParams(18);

%Temperature effect on photosynthesis PMAX
PGRED = CropParams(19:27);
TMPG = CropParams(28:36);

%plot(TMPG,PGRED);
%xlabel('Temperature in C');ylabel('Temperature correction factor on PMax');

%DVS effect on photosynthesis PMAX and SLA
 DVSArray = CropParams(37:45);
 PDVS = CropParams(46:54);
 SLA = CropParams(55:63);
%plot(DVSArray,SLA); xlabel('DVS');ylabel('SLA');
%plot(DVSArray,PDVS); xlabel('DVS');ylabel('Leaf senescence factor on PMax');


%Critical concentrations are necessary for optimal growth and development
%Shoot critical N concentration in % in relation to DVS
%DVS at  biom measures
%DVS.meas <- c(0,0.40,0.64,0.82,1.08,1.69,2)
%Corresponding DAP in 2009, just for calibration purposes, DAP at  biom measures
%DAP.meas  <-  c(0,43,55,64,73,91,100)
%Nc  in percent needed for 100 % shoot production
%Shoot.Nc  <- c(5.5,5.5,5.25,5.0,3.0,2,1)
%Shoot.Nc <- c(4.5,4,3.5,3,2.5,2.25,2)
DVS_meas = CropParams(64:69);
Shoot_Nc = [5.5 5.5 5 4 3.5 3];%CropParams(70:75)*1.25;
%plot(DVS.meas,Shoot.Nc,type='l',xlab='DVS',ylab='Optimal N Conc in Shoot in %')

%Curd critical N concentration in %
Curd_Nc = 3.5;%CropParams(76);

%Root critical N concentration in %
Root_Nc = CropParams(77);

%FertDosis is necessary to subset experimental data of 2009
FertDosis = CropParams(4);

%VPD effect on photosynthesis
VPDL = CropParams(78);
CK = CropParams(79);

%Plant density in pl/m^2
PLM2 = 2.85;%CropParams(5);
    %proportionality factor for coverage calculated on PLM2=4 to PLM2 = 2.8
    PropFact =1;% 2.8 / PLM2;

    if DAP == 1
%initialisation parameters, estimated in ReadDSSATCBT.R
    %initial DW leaves in g
    DMLeaf = CropParams(6);
    %initial leaf area in m^2
    TotLeafArea = DMLeaf * SLA(1)/10000;
    %Initial DW stem in g
    DMStem = CropParams(7);
    %Initial DW curd in g
    DMCurd = CropParams(8);
    %Initial nitrogen content
    InitN = 5.4;%
    end
%Maintenance respiration parameters
%Temperature effect parameter
Q10 = CropParams(80);
%Respiration parameters for Leaf and Curd (F)g(CH2O)/g(DW) day
RMRL = CropParams(81);
RMRF = CropParams(82);

%phenological parameters
    %Effective T sum for end vegetative phase with Tb=4 in Cd
    TeSumCurdIni = CropParams(12);
    %Effective T sum for end curd maturity with Tb=4 in Cd
    TeSumCurdMat = CropParams(13);

%Growth efficiency
GREF = CropParams(83);

%old compartimentation functions
DVSArray2 = CropParams(84:94);
LeafFrac = CropParams(95:105);
StemFrac = CropParams(106:116);
CurdFrac = CropParams(117:127);
RootFrac = CropParams(128);
%Residual leaves staying on the field is 55 % of the shoot
LeafResFrac = CropParams(129);

 %subplot(1,3,1); plot(DVSArray2,LeafFrac); xlabel('DVS');ylabel('Leaf Fraction');
 %subplot(1,3,2); plot(DVSArray2,StemFrac); xlabel('DVS');ylabel('Stem Fraction');
 %subplot(1,3,3); plot(DVSArray2,CurdFrac); xlabel('DVS');ylabel('Curd Fraction');

%Adjust tempsum!!!!!

%#new calibration curd initialisation at DVS=1
%#DVSArray2 = c(0,0.15,0.5,0.6,0.8,1,1.2,1.4,1.6,1.8,2)
%#LeafFrac = c(rep(0.8,2),0.85,0.85,0.90,0.9,0.75,0.65,rep(0.5,3))
%#StemFrac=c(rep(0.2,2),0.15,0.15,0.1,rep(0.1,6))
%#CurdFrac=c(rep(0,6),0.15,0.25,0.40,rep(0.4,2))
%#RootFrac = 0.09
%#win.graph()
%#par(mfrow=c(3,1))
%#plot(DVSArray2,LeafFrac,type='l',xlab='DVS',ylab='Leaf Fraction')
%#plot(DVSArray2,StemFrac,type='l',xlab='DVS',ylab='Stem Fraction')
%#plot(DVSArray2,CurdFrac,type='l',xlab='DVS',ylab='Curd Fraction')
%#LeafFrac + StemFrac + CurdFrac



%Deadleaf fraction in relation to Te
Tef = CropParams(130:138);
Kd = CropParams(139:147);
%plot(Tef,Kd); xlabel('Te'); ylabel('Kd dead leaves');

%Leaf coverage in relation to DVS
DVSAr = CropParams(148:168);
%Leaf coverage in percent
LeafCoverage = CropParams(169:189);%plot(DVSAr,LeafCoverage); xlabel('DVS'); ylabel('Leave Coverage in %');

if DAP == 1
%Initial LAI
LAICov = TotLeafArea * PLM2 / (LeafCoverage(1)/100);
end


%Root DM partitioning function for 11 virtual layers
RDPF = CropParams(190:200);
%Virtual root depth
VirRD = CropParams(201:211);
%plot(RDPF,-VirRD); ylabel('Virtual Rooting Depth'); xlabel('Root DM partitioning');

if DAP == 1
%initialise results file
results = zeros(NDAYS, 38);

%Initialisation of Effective Temperature Sum (TeSum)
TeSum = 0; 
%Initialisation of integration
TotDMLeaf = DMLeaf;
TotDMStem = DMStem;
TotDMCurd = DMCurd;
TotDMRoot = 0;
TotDMShoot = DMStem+DMLeaf+DMCurd;
DVS = CropParams(10);
TotDMGreenLeaf = 0;
TotDMDeadLeaf = 0;
TotNperLeaf = 0;
TotNperStem = 0;
TotNperCurd = 0;
CumTranspiration = 0;
CumEvaporation = 0;
TotNDemandLeafStem = 0;
TotNDemandCurd = 0;
TotNDemandRoot = 0;
CumPlantNUptakeActualCurd = 0;
CumPlantNUptakeActualLeaf = TotDMLeaf*InitN/100 ;
CumPlantNUptakeActualStem =  TotDMStem*InitN/100;
CumPlantNUptakeActualLeafStem = CumPlantNUptakeActualLeaf + CumPlantNUptakeActualStem;
CumPlantNUptakeActualRoot = 0 ;

%Initialisation
 Holding(1,:) = [0,0,0,0,0,0,0,0,0];

%Initial root depth at DAP=0 in cm
TotRootDepth2 = CropParams(9);;

else
% DATE = output_crop{ncrop}(end,1);
% DAP = output_crop{ncrop}(end,2)
GP = output_crop{ncrop}(end,3);
RMAINT = output_crop{ncrop}(end,4);
MeanDayTemp= output_crop{ncrop}(end,5);
Te= output_crop{ncrop}(end,6);
TeSum= output_crop{ncrop}(end,7);
DVS= output_crop{ncrop}(end,8);
TotDM= output_crop{ncrop}(end,9);
DMRoot= output_crop{ncrop}(end,10);
DMLeaf= output_crop{ncrop}(end,11);
LeafArea= output_crop{ncrop}(end,12);
DMStem= output_crop{ncrop}(end,13);
DMCurd= output_crop{ncrop}(end,14);
TotDMLeaf= output_crop{ncrop}(end,15);
TotDMStem= output_crop{ncrop}(end,16);
TotDMCurd= output_crop{ncrop}(end,17);
TotDMRoot= output_crop{ncrop}(end,18);
TotDMShoot= output_crop{ncrop}(end,19);
TotLeafArea= output_crop{ncrop}(end,20);
LAICov= output_crop{ncrop}(end,21);
LAI= output_crop{ncrop}(end,22);
TotDMGreenLeaf= output_crop{ncrop}(end,23);
TotDMDeadLeaf= output_crop{ncrop}(end,24);
TotRootDepth1= output_crop{ncrop}(end,25);
RCil= output_crop{ncrop}(end,26);
Transpiration= output_crop{ncrop}(end,27);
Evaporation= output_crop{ncrop}(end,28);
CumTranspiration= output_crop{ncrop}(end,29);
CumEvaporation= output_crop{ncrop}(end,30);
LeafStemNDemand= output_crop{ncrop}(end,31);
CurdNDemand= output_crop{ncrop}(end,32);
RootNDemand= output_crop{ncrop}(end,33);
PlantNDemand= output_crop{ncrop}(end,34);
TotNDemandLeafStem= output_crop{ncrop}(end,35);
TotNDemandCurd= output_crop{ncrop}(end,36);
TotNDemandRoot= output_crop{ncrop}(end,37);
TotNDemandPlant= output_crop{ncrop}(end,38);

CumPlantNUptakeActualCurd = output_cropNup{ncrop}(end,1);
CumPlantNUptakeActualLeaf =output_cropNup{ncrop}(end,2) ;
CumPlantNUptakeActualStem =  output_cropNup{ncrop}(end,3);
CumPlantNUptakeActualLeafStem = output_cropNup{ncrop}(end,4) ;
CumPlantNUptakeActualRoot = output_cropNup{ncrop}(end,5) ;
end

%for JDAY = 1:NDAYS
% JDAY = t-plant_date(ncrop)+1;
% 	TIME = (JDAY-1) * DELT;
% 	DATE = t;%mod(TIME+NSTART,365);
% 	DAP = JDAY;
    % initialize variables that accumulate during the day
       GP = 0;
       RMAINT = 0;
       Transpiration = 0;
       Evaporation = 0;
    %parameter to calculate average day temperature
       Tsum = 0;
   
       % *********************************************************************
       % fast loop (if time-step = 1 hour, NFAST = 24)
       % ******************************************************************
        for JF = 1:NFAST
            % TFAST is the hour of the day (h)
            TFAST = (JF-1) * 24/NFAST;
            % ***************************************************************
            % input of hourly climate data --- source(RefClim) ---
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
            Ps = 0.61078 * exp((17.2694 * TMPA)/(TMPA+237.3));
            %Calculate actual vapor pressure PV
            Pv = RH/100 * Ps;
            VPD = Ps-Pv;
            VPD_nok = VPD*(10^3);
            CO2L = CropParams(213);
            PAR = solrad * CropParams(214);
            PPFD = PAR * CropParams(215);
        
            % ***************************************************************
            % calculation of gross photosynthesis --- source('PHOTO2.R') ---
            % ***************************************************************
            
		     GPF = 0;
			% effect of CO2 on Pmax (from Gainesville)
			PMAX = TAU1 * CO2L;
			if CO2L>1500
                PMAX = (TAU1 * 1500) + (TAU2 * (CO2L-1500));
            end
            % reduction of Pmax at extreme temperatures
			PMAX = PMAX * interp1(TMPG,PGRED,TMPA);
            
            %reduction of Pmax in relation to DVS, senescence
			PMAX = PMAX * interp1(DVSArray,PDVS,DVS);

			%Reduction of PMax in relation to N-dosis
%			PMAX = PMAX * interp1(dosis,NCor,FertDosis);
if isempty (output_cropNup{1,ncrop})==0
PMAX = PMAX * output_cropNup{1,ncrop}(end,11);
end

			  %reduction of XK with leaf area (from Olesen, 1997)
            if TotLeafArea < 0.235
                XK = 0.90 - 3.00 * TotLeafArea + 6.38 * TotLeafArea^2;
            else XK = XKK;
            end
            
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
			      % conversion of GPF from µM/m2-s into g/m2-day   %µM/m2-s x 0.000044g/µM x 3600s/h x 24h/d = 3.8016 g/m2-day
			      GPF = GPF * 3.8016/PLM2;
            end

	%**************************************************************
    % calculation of maintenance respiration 
   	%**************************************************************
     	% Effect of temperature on maintenance respiration on hourly basis
                RMAINTF = (RMRL * (DMLeaf+DMStem)+RMRF * DMCurd) * (Q10^(0.1 * TMPA-2.0));

                %integration of variables on 24 hours
                GP = GP + GPF * DTFAST;
                RMAINT = RMAINT + RMAINTF * DTFAST;

%*****************************************************************
%Crop transpiration on hourly basis 
%*****************************************************************
    %Internal resistance in s/m
    ri = CropParams(216);
    re = CropParams(217);
    %Specific heat of air J/kg K
    Cp = CropParams(218);   
    %AtmPressure in Pa
    AtmP = CropParams(219);

    %Short wave radiation calculations W/m^2
    Is = solrad * CropParams(220);

    %Long wave rad calculations  W/m^2
    Il = solrad * CropParams(221);

    %Air density in kg/m^3
    rho.air = 100000/(287*(TMPA+273.16));

    %Slope of saturationn curve Pa/C
    delta = 41.45*exp(0.06088*TMPA);    

    %Heat of evaporisation J/kg
    lambda = 2502535.259-(2385.76*TMPA);

    %Psychrometric constant along to Moneith p 181 Value=67 Pa/C at 20 C
    gamma = (Cp*AtmP)/(lambda*0.622);

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
    x = exp(-0.24*LAICov);
    %As absorption coefficient short wave radiation, smoothing rules Stanghellini
    As = 0.83*(1-0.70*x)*(1+0.58*x*x)*(0.88-x*x + 0.12*x^(8/3)) ;
    %Al absorption coefficient long wave radiation
    Al = 0.83*(1-0.70*x)*(1-x^(8/3));
    Rn = As*Is + Al*Il;

    %H = convection loss by canopy if not neglected
    %H=(LAICov*rho.air*Cp*(TMPA-T.Canopy))/re
    H = 0;
 
    %Transpiration heat in W/m^2 covered ground surface multiplied with
    %leaf area and plantdensity
    LE = (1/(1+(delta/gamma) + (ri/re))*( ((delta/gamma)* (Rn-H)) + (rho.air*Cp*VPD_nok)/(gamma*re)))*TotLeafArea * PLM2; 
    %Transpiration in kg/s m^2 coverage
    E = LE/lambda;
    EHour = E * 3600;

    %integration of variables on 24 hours
    Transpiration = Transpiration + EHour;



%--------------------------------------------------------------------------
%Soil evaporation. Classical Penmann equation (Monteith, Unsworth, 1990) with resistance 
%--------------------------------------------------------------------------
    %Surface conductance m/s
    rs = 1/70;%CropParams(222);
    albedo = CropParams(223);
    %ESoil in kg / m^2 s
    ESoil = (delta*solrad*albedo + rho.air*Cp*VPD_nok*rs) / (lambda*(delta + gamma));
    ESoilHour = ESoil*3600;%*(1-(interp1(DVSAr,LeafCoverage,DVS) / 100));

    %Integration of variables on 24 hours
        Evaporation = Evaporation + ESoilHour;    

%--------------------------------------------------------------------------        
%Calculate temperature sum for average day temperature calculations
%--------------------------------------------------------------------------                
    %calculate average day temperature
    Tsum = Tsum + TMPA;
    
               %******************************************************************
        end    %      end of fast loop
               %***********************************************************************


%#GP and RMAINT in g CH2O per plant
    %#Transfert to Biomass, calculate Total DM in g per plant per dag
    TotDM = (GP - RMAINT) * GREF; 

%calculate Development Stage (DVS)
MeanDayTemp = Tsum/24;
Td = MeanDayTemp - Tb;
    %effective temperature set to zero if negative, non-decreasing Tsum
    if Td > 0
        Te = Td;
    else Te = 0;
    end
	%calculate Temperature sum
    if DAP ==1
	TeSum = TeSum + Te;
    else
       	TeSum = output_crop{ncrop}(end,7) + Te;
    end
    %DVS for leaves + stem
    DVS_VegC = TeSum/TeSumCurdIni;
    if DVS_VegC > 1
        DVS_Veg = 1;
    else DVS_Veg =DVS_VegC;
    end
    %DVS for curd
    if DVS_VegC > 1
        TeSumC = TeSum - TeSumCurdIni;
    else TeSumC = TeSum;
    end
    if DVS_VegC > 1
        DVS_Curd = TeSumC/TeSumCurdMat;
    else DVS_Curd = 0;
    end
    %Overall DVS
    DVS = DVS_Veg + DVS_Curd;
    
    DVS = min(DVS,2);
    
%Compartimentation over leaf, stem and curd
%Growth rates g per planr per day
DMLeaf = TotDM * interp1(DVSArray2,LeafFrac,DVS)/(1 + RootFrac);
DMStem = TotDM * interp1(DVSArray2,StemFrac,DVS)/(1 + RootFrac);
DMCurd = TotDM * interp1(DVSArray2,CurdFrac,DVS)/(1 + RootFrac);
    
 %Calculate DM of green and dead leaves
    Kdl = interp1(Tef,Kd,Te);
    if DVS > 1
        DMDeadLeaf = DMLeaf*(exp(Kdl)-1);
    else DMDeadLeaf = 0;
    end
    if DVS > 1
        DMGreenLeaf = DMLeaf-DMDeadLeaf;
    else DMGreenLeaf = DMLeaf;
    end
    
%compartimentation to roots, 9 % of shoot and constant
DMShoot = DMGreenLeaf - DMDeadLeaf + DMStem + DMCurd;
DMRoot = DMShoot * RootFrac;


%calculate leaf area out of SLA-DVS curve
LeafArea = DMGreenLeaf * interp1(DVSArray,SLA,DVS)/10000;

%integrate to total DM
TotDMLeaf = TotDMLeaf + DMLeaf;
TotDMStem = TotDMStem + DMStem;
TotDMCurd = TotDMCurd + DMCurd;
TotDMShoot = TotDMShoot + DMShoot;
TotDMRoot = TotDMRoot + DMRoot;
TotDMGreenLeaf = TotDMGreenLeaf + DMGreenLeaf;
TotDMDeadLeaf = TotDMDeadLeaf + DMDeadLeaf;
TotLeafArea = TotLeafArea + LeafArea;



%direct LAI calculation
LAI = TotLeafArea * PLM2;
%direct LAI calculation per m^2 coverage
LAICov = TotLeafArea*PLM2 / (interp1(DVSAr,LeafCoverage,DVS)/ 100);

%Integrate transpiration and evaporation to cumulative values
%if DAP == 1
CumTranspiration = CumTranspiration + Transpiration;
CumEvaporation = CumEvaporation + Evaporation;
% else
%  CumTranspiration = output_crop{ncrop}(end,32)+ Transpiration;
% CumEvaporation = output_crop{ncrop}(end,33) + Evaporation;   
% end

%###############################################################################
%#Root volume and length calculations
%##############################################################################
    TLayer = CropParams(224);
 %Specific root length in cm/g
  SRL = CropParams(212);
 %methode 1 forcing function of depth in relation to DVS
 %RootDepth follows  sigmoidal curve like LeafCoverage, to be calibrated
 TotRootDepth1 = interp1(DVSAr,LeafCoverage,DVS);

%Results of measurements done on roots cauliflower gamma dist
rootDAP = [0 20 40 68 120];
z1_param = [2.98 2.98 1.77 1.68 1.68]; 
z2_param = [1.36 1.36 5.24 7.57 7.57];
rootdensparam1 =interp1(rootDAP,z1_param,DAP); %Find current parameter1
rootdensparam2 =interp1(rootDAP,z2_param,DAP); %FInd current parameter2
depth = 1:TLayer:100;
rdens = gampdf(depth,rootdensparam1,rootdensparam2); %Calculate rdens over profile
maxindex = find(rdens == max(rdens)); 
minindex = find(rdens< 0.0001); %Calculate depths with low root densities
TotRootDepth1 = minindex(min(find((minindex>maxindex)==1))); %find the first depth after maxindex that has a low root densiti
%==> this depth is defined as the rooting depth;


%Calculate number of layers, layer = TLayer in cm

    NumLayer = round(TotRootDepth1/TLayer);
    %#Calculate DM per layer, partition along Root DM Partitioning function
    RDMLa(1:NumLayer) = 0;
    %Root Length per layer
    RLengthLa(1:NumLayer) = 0;
    %Root Volume per layer, using big root model with radius rr in cm
    rr = CropParams(225);
    RVolumeLa(1:NumLayer) = 0;
    %Equivalent radius of the volume with tickness TLayer, in the most simple model this is a cilinder per layer
    ERLa(1:NumLayer) = 0;
    %Root Length Density per layer in cm/cm^3
    RLDLa(1:NumLayer) = 0;

 %Initialize RootRes
    RootRes = zeros(NumLayer, 6);

%Loop over layers to calculte sum of RDPF
    RDDistProp(1:NumLayer) = 0;
    RDDistSum = 0;
    for j = 1:NumLayer
        RDDist = VirRD * NumLayer;
        RDDistProp(j) = interp1(RDDist,RDPF,j);
        RDDistSum = RDDistSum + RDDistProp(j);
    end

%Loop over layers
    for i = 1:NumLayer
        RDDist = VirRD * NumLayer;
        RDMLa(i) = TotDMRoot * interp1(RDDist,RDPF,i)/RDDistSum;
        RLengthLa(i) = RDMLa(i) * SRL;
        RVolumeLa(i) = RLengthLa(i) * pi * rr^2;
        ERLa(i) = sqrt(RVolumeLa(i)/(TLayer*pi));
        %For now we use the actual volume to calculate Root Length Density, could also use cilider with RCil as radius and TLayer as tickness
        if RVolumeLa(i) ~= 0
            RLDLa(i) = RLengthLa(i);
        else RLDLa(i) = 0;
        end
        RootRes(i,:) = [i,RDMLa(i), RLengthLa(i),RVolumeLa(i),ERLa(i),RLDLa(i)];
    end


%Calculate Radius of enclosing cilinder RCil
        %Root Length Density per layer in cm/cm^3 for enclosing cilinder
        RCil = max(ERLa);
        RLDRCil = RLengthLa /((RCil^2) * pi * TLayer);

        RDATE(1:NumLayer) = DATE; RRDATE = RDATE';
        RDAP(1:NumLayer) = DAP; RRDAP = RDAP';
        RLDRCil =RLDRCil';
        RootResDay = [RRDATE,RRDAP,RootRes,RLDRCil];


%    Holding = [Holding;RootResDay];
    %###############################################################################
    %End of root module
    %##############################################################################
demandoption =1;
if demandoption ==1
%N demand for optimal growth in leaves in g per plant
LeafStemNDemand = interp1(DVS_meas,Shoot_Nc,DVS) * (TotDMGreenLeaf + TotDMDeadLeaf+TotDMStem)/100-...
    CumPlantNUptakeActualLeafStem;
LeafStemNDemand = max(0,LeafStemNDemand);
LeafNDemand = LeafStemNDemand*((TotDMGreenLeaf+TotDMDeadLeaf)/(TotDMGreenLeaf+TotDMDeadLeaf+TotDMStem));
StemNDemand = LeafStemNDemand*((TotDMStem)/(TotDMGreenLeaf+TotDMDeadLeaf+TotDMStem));
%N demand in the curds in g per plant
CurdNDemand = Curd_Nc*(TotDMCurd/100)-CumPlantNUptakeActualCurd;
CurdNDemand = max(CurdNDemand,0);
%N demand in the roots in g per plant
RootNDemand = Root_Nc*TotDMRoot/100-CumPlantNUptakeActualRoot;
RootNDemand = max(RootNDemand,0);
else
%N demand for optimal growth in leaves in g per plant
LeafStemNDemand = interp1(DVS_meas,Shoot_Nc,DVS) * (DMGreenLeaf + DMDeadLeaf + DMStem) / 100;
%N demand in the curds in g per plant
CurdNDemand = Curd_Nc*DMCurd/100;
%N demand in the roots in g per plant
RootNDemand = Root_Nc*DMRoot/100;
end

%Calculate tota N demand per day in g per plant
PlantNDemand = LeafStemNDemand + CurdNDemand + RootNDemand;
if PlantNDemand ~=0
ratiostem = StemNDemand/PlantNDemand;
ratioleaf = LeafNDemand/PlantNDemand;
ratioroot = RootNDemand/PlantNDemand;
ratiocurd = CurdNDemand/PlantNDemand;
else
ratioleaf = 0;
ratiostem = 0;
ratioroot = 0;
ratiocurd = 0;
    
end

%Integrate Tot N demand  
TotNDemandLeafStem = TotNDemandLeafStem + LeafStemNDemand;
TotNDemandCurd = TotNDemandCurd + CurdNDemand;
TotNDemandRoot = TotNDemandRoot + RootNDemand;

TotNDemandPlant = TotNDemandLeafStem + TotNDemandCurd + TotNDemandRoot;

if DAP>70
    disp('drf')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
PlantNUptakeActualCurd = PlantNUptakeActual*ratiocurd;
PlantNUptakeActualLeaf = PlantNUptakeActual*ratioleaf;
PlantNUptakeActualStem = PlantNUptakeActual*ratiostem;
PlantNUptakeActualLeafStem = PlantNUptakeActualLeaf + PlantNUptakeActualStem;
PlantNUptakeActualRoot = PlantNUptakeActual*ratioroot;

CumPlantNUptakeActualCurd = CumPlantNUptakeActualCurd +PlantNUptakeActualCurd ;
CumPlantNUptakeActualLeaf = CumPlantNUptakeActualLeaf+PlantNUptakeActualLeaf ;
CumPlantNUptakeActualStem = CumPlantNUptakeActualStem+PlantNUptakeActualStem ;
CumPlantNUptakeActualLeafStem = CumPlantNUptakeActualLeafStem+PlantNUptakeActualLeafStem  ;
CumPlantNUptakeActualRoot = CumPlantNUptakeActualRoot +PlantNUptakeActualRoot ;

%Two approaches : nitrogen concentration in the freshly formed organs
%(ActPlantNperc) of nitrogen concentration averaged over the shoot
%(PercNShoot)
%ActPlantNPerc = (PlantNUptakeActual/TotDM)*100; % Nitrogen concentration in newly formed organs
PercNShoot = CumPlantNUptakeActualLeafStem/(TotDMLeaf+TotDMStem)*100;
%Calculate NFac. N correction factor on Pmax in PS model 
%NFac = 1 - (interp1(DVS_meas,Shoot_Nc,DVS) - ActPlantNPerc) / (interp1(DVS_meas,Shoot_Nc,DVS) - interp1(DVS_meas,Shoot_Nm,DVS));
%NFac = 1;
%NFac = 1 - (interp1(DVS_meas,Shoot_Nc,DVS) - PercNShoot) / (interp1(DVS_meas,Shoot_Nc,DVS) - interp1(DVS_meas,Shoot_Nm,DVS));
NFac = PercNShoot/ interp1(DVS_meas,Shoot_Nc,DVS);
NFac = min(NFac,1);
NFac = max(NFac, 0.05);
%NFac = 1-sqrt(1-NFac); %as in WAVE manuel page 6-12

%Create N-pool for retranslocation 
%if PlantNDemand <= NSupply
%NPool = NSupply-PlantNDemand
%else
%NPool = 0
%end

%Actual N under luxus consumption of 20 %
%ActPlantNPercL = ((PlantNDemandActual + 0.2*NPool)  / TotDM) * 100
%ActPlantNUptake = ActPlantNPerc*TotDM

%Creation of retranslocation pool
%ActNpool = NPool*0.2
%ActPlantNUptake = ActPlantNUptake + ActNPool

crop_result = [DATE,DAP,GP,RMAINT,MeanDayTemp,Te,TeSum,DVS,TotDM,DMRoot,DMLeaf,LeafArea,DMStem,DMCurd,TotDMLeaf,TotDMStem,TotDMCurd,TotDMRoot,TotDMShoot,TotLeafArea,LAICov,LAI,TotDMGreenLeaf, ...
                        TotDMDeadLeaf,TotRootDepth1,RCil,Transpiration,Evaporation,CumTranspiration,CumEvaporation,LeafStemNDemand,CurdNDemand,RootNDemand,PlantNDemand,TotNDemandLeafStem,TotNDemandCurd, ...
                        TotNDemandRoot,TotNDemandPlant];
crop_result2 = [CumPlantNUptakeActualCurd,CumPlantNUptakeActualLeaf,CumPlantNUptakeActualStem, CumPlantNUptakeActualLeafStem,...
    CumPlantNUptakeActualRoot,PlantNUptakeActualCurd,PlantNUptakeActualLeaf,PlantNUptakeActualStem, PlantNUptakeActualLeafStem,...
    PlantNUptakeActualRoot,NFac,PercNShoot];