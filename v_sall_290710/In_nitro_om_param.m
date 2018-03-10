
function[om_param, miner_param, plant_uptake_param]=In_nitro_om_param(ncs);

%IN-NITROGEN_PARAM: parameters for the nitrogen module 
%IN:
%soil_parameters
%OUT:
%nitrogen module parameters
%CALL: none
%CALLED BY:WAVE_MAT main programme
%-------------------------------------------------------------------------
% Mamadou SALL 17/03/09

%%rorad, rd0, g, rnmaxp, rdens0, alfa_rdens
%% RNMAXP obtenu de l'étude expérimentale = 1 mg/cm²
plant_uptake_param=[0.01;0.01;1;1.9;0.8;0.5];

%%% ro, fe, fh
miner_param=[10;0.3;0.4];

           %% knitrif, kdenit, khydro, kvolat, klito, kmano, khumo,

               
lay1=ones(10,1)* [ 5.82     0.01     0.00    0.0    0.035  0.035  0.00007];
lay2= ones(20,1)* [5.82      0.01     0.00   0.0    0.035  0.035  0.00007];
lay3= ones(69,1)* [0.       0.01      0.    0.0    0.0    0.0    0.0000];
lay4=ones(1,1)* [0.       0.0      0.    0.0    0.0    0.0    0.0000];
om_param=[lay1;lay2;lay3;lay4];





%%Plant uptake param:
%rorad=plant_uptake_param(1) (dimension L);       
%rdo=plant_uptake_param(2)/ (dimension L);        
%g= plant_uptake_param(3)(exprimé en fraction sans dimension);            
%rnmaxp=plant_uptake_param(4) (dimension M/L^2/T);  
%rdens0=plant_uptake_param(5) (dimension L/L^3);   
%alfa_rdens=plant_uptake_param(6)(sans dimension);  
           
           
           

%nitro_param:

%rnmaxp                      %(MG M-2 DAY-3)maximum Nuptake
%rdens0                      %(CM L-1)root density at the soil surface
%alfa_rdens                  %(-) coefficient for calculation %%rdensi with rdens0 in script rooting_depth
%rorad                      %%%(MM)mean root radius
%rdo                          %(MM)average distance between soil solution and root surface
%ro                          %(-) C/N ratio
%fe                          %(-)efficiency %constant of the carbon assimilation
%fh       
%kdenit                     denitrification constant
%kvolat                   volatilisation constant
%knitrif                   nitrification constant
%khydro
%g                      %fraction of growing saison of potential uptake (-) a 








