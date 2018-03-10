
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
plant_uptake_param=[0.1;0.0;0.5;0.2224;0.1;0.5];

%%% ro, fe, fh
miner_param=[10;0.3;0.4];

%% knitrif, kdenit, khydro, kvolat, klito, kmano, khumo,
om_param=[     5.82    0.01     0.    0.0    0.035  0.035  0.00007       %%%5.82     0.01
               0.0       0.01     0.    0.0    0.035  0.035  0.00007            %%% 0.0      0.01
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
                0.       0.0      0.    0.0    0.0    0.0    0.0000
               0.       0.0      0.    0.0    0.0    0.0    0.0000];
               



%rorad=plant_uptake_param(1);       
%rdo=plant_uptake_param(2);        
%g= plant_uptake_param(3);           
%rnmaxp=plant_uptake_param(4);  
%rdens0=plant_uptake_param(5);   
%alfa_rdens=plant_uptake_param(6);  
           
           
           

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








