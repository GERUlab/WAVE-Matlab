function [soil_om,ptup,cupt,pthyd,chyd,ptnit,cnit,ptvol,cvol,ptden,cden,ptsnorg,ptscorg,...
    tnorgs,tcorgs,ptmin,cmin]=In_initial_om(ncs)
%NITROGEN MODULE
%IN
%soil_parameters, ncs(number of compartiments)
%OUT
%C_MANP: carbone in manure (mg m-2)
%N_MANP: nitrogen in manure (mg m-2)
%C_LITP: carbone in litter (mg m-2)
%N_LITP: nitrogen in litter (mg m-2)
%C_HUMP: carbone in manure (mg m-2)
%N_HUMP: nitrogen in humus (mg m-2)
% ptsnorg : cumulative organic nitrogen sink
%ptscorg :  cumulative organic carbon sink
%ptupt : profile total solute uptake
%cupt : cumulative plant uptake in each compartiment
%ptmin :  profile total mineralisation
%cmin : cumulative mineralisation of nitrogen in each compartiment
%pthyd : profile total hydrolysis
%chyd : cumulative hydrolysis in each compartiment
%ptnit : profile total nitrification
%cnit : cumulative nitrification in each compartiment
%ptvol : profile total volatilisation
%cvol : cumulative volatilisation in each compartiment
%ptden : profile total denitrification
%cden : cumulative denitrification in each compartiment
%tnorgs : initial organic nitrogen profile total
%tcorgs : initial organic carbon profile total
%soil_om : actual organic carbon and notrogen in the profile

%CALL:none
%CALLED BY: WAVE_MAT
%--------------------------------------------------------------------------
%M.SALL 26/11/08

 %c_manp=soil_om(:,1);
 %n_manp=soil_om(:,2);
 %c_litp=soil_om(:,3);
 %n_litp=soil_om(:,4);
 %c_hump=soil_om(:,5);
 %n_hump=soil_om(:,6);
 
%%Résultats d'analyse de sol
%0.14% de carbone total sur le sol entre 0-10cm
 %0.12% d'azote total entre 0-10 cm
 %0.05% de COT entre 10 et 30 cm 
%0.02% d'azote total entre 10 et 30 cm  
 
  %0.14% de carbone total sur le sol entre 0-10cm soit pour une densité de
 %%1.52g/cm³ la concentration de C est 2128 kg par ha ou encore 21.28 mg/cm²
 %% si C/N est égal à 10% l'azote organique est 2.128 mg/cm² entre 0 et 10  %% cm
 % Avec le même raisonnement entre 10 et 30 cm
%% On trouve concentration de C = 7.6 mg/cm²
%%% Concentration de N= 0.76 mg/cm²

 
         % C_MAN    N_MAN    C_LIT    N_LIT        C_HUM    N_HUM
%C ----   ------    ------    ------    ------    ------ ------
 
  layer1= ones(10,1)*[2.128,0.2128,2.128,0.2128,17.024,1.7024]*0.1;     
  layer2= ones(20,1)*[0.76,0.076,0.76,0.076,6.08,0.608]*0.05;
  layer3= ones(70,1)*[0.1,0.01,0.1,0.01,0.1,0.01]*1E-8;
  soil_om=[layer1;layer2;layer3];
  
  
  
  
ptup = zeros(3,1); cupt = zeros(ncs,3);
ptmin = zeros(3,1); cmin = zeros(ncs,3);
pthyd = 0; chyd = zeros(ncs-1,1); 
ptnit = 0; cnit = zeros(ncs-1,1); 
ptvol = 0; cvol = zeros(ncs-1,1); 
ptden = 0; cden = zeros(ncs-1,1); 
ptsnorg = 0; %%initial sink is zero
ptscorg = 0; %%initial sink is zero
tcorgs = sum((soil_om(:,1)) + soil_om(:,3) + soil_om(:,5));                
tnorgs = sum((soil_om(:,2)) + soil_om(:,4) + soil_om(:,6));    
          
 
       
 %c_manp=soil_om(:,1);
 %n_manp=soil_om(:,2);
 %c_litp=soil_om(:,3);
 %n_litp=soil_om(:,4);
 %c_hump=soil_om(:,5);
 %n_hump=soil_om(:,6);
 
 
 
 
 