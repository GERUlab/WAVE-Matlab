function  [rdens,drza]= calc_rdens(plant_uptake_param,ncs,dx,drz,t)

global CC
%Based on van Noordwijk(1991) and Feddes (1973)

rorad=plant_uptake_param(1);       %%rorad=nitro_param(4);
rdo=plant_uptake_param(2);         %%rdo=nitro_param(5);
g= plant_uptake_param(3);           %%%g= nitro_param(13);  
rnmaxp=plant_uptake_param(4);  %%rnmaxp=nitro_param(1);

% % %Average
% rdens0=2.8;   %%rdens0=nitro_param(2);
% alfa_rdens=0.031;  %%alfa_rdens=nitro_param(3);
% % 
% % %High
rdens0=4;   %%rdens0=nitro_param(2);
alfa_rdens=0.032;  %%alfa_rdens=nitro_param(3);
% % 
% % %  %Low
%  rdens0=1.6;   %%rdens0=nitro_param(2);
%  alfa_rdens=0.025;  %%alfa_rdens=nitro_param(3);

%if crop_type = 

%Determine maximum root length (root_length_max in the script)
root_length_max = 75;

if root_length_max ==0
    error('the maximum rooting depth is zero, check inputfile of crops')
end

if root_length_max > ncs*dx
    error('the maximum rooting depth exceeds length of profile, check inputfile of crops')
end

drza = drz;
drza_comp = ceil(drza/dx);
% i = 1:drza_comp;
% rdens(i) = rdens0*exp(-alfa_rdens.*i.*dx)*drza/root_length_max;
% i = drza_comp+1:ncs;
% rdens(i) = 0;

[CK, CumEvaporation,CumTranspiration,CurdFrac,Curd_Nc,DELT,DMCurd,DMLeaf,DMStem,...
    DTFAST,DVS,DVSAr,DVSArray,DVSArray2,DVS_meas,GREF,...
    Holding,Kd,LAICov,LeafCoverage,LeafFrac,LeafResFrac,NDAYS,NFAST,PDVS,...
    PGRED,PLM2,PropFact,Q10,QE,RDPF,RMRF,RMRL,RootFrac,Root_Nc,SLA,Shoot_Nc,StemFrac,...
    TAU1,TAU2,TMPG,Tb,TeSum,TeSumCurdIni,TeSumCurdMat,Tef,TotDMCurd,TotDMDeadLeaf,...
    TotDMGreenLeaf,TotDMLeaf,TotDMRoot,TotDMShoot,TotDMStem,TotLeafArea,...
    TotNDemandCurd,TotNDemandLeafStem,TotNDemandRoot,TotRootDepth2,VPDL,VirRD,XKK,XM,NFac,Shoot_Nm,initN] = in_cauli_parames();

NumLayer = drza_comp;
    RDDistProp(1:NumLayer,1) = 0;
    RDDistSum = 0;
    for j = 1:NumLayer
        RDDist = VirRD * NumLayer;
        RDDistProp(j) = interp1(RDDist,RDPF,j);
        RDDistSum = RDDistSum + RDDistProp(j);
    end

RDDistProp (NumLayer+1:ncs) = 0;
rdens = RDDistProp;
%rdens = rdens./CC;




    