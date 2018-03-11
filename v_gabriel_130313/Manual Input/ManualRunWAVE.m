%Run WAVEmat without GUI
close all;clear all;clc
%for scenario_number = 1:10; %For loop that initiates scenario-analysis
%% Write input-parameterfiles (optional)    
write_inputfile_1punt
write_profile_data
scenario_number = 1;

%% Define input-files (each in line can represent a new scenario)
inputfile_soil = ('outputmatrix.xls');
inputfile_profile_soil = ('outputmatrix_profile.xls');
inputfile_crop = ('CropParams_Default_under40.txt');
inputfile_leek = ('LeekParams_Default_under40.txt');
inputfile_climate = 'y2009_2010_2011.txt';

% clearvars -except scenario_number inputfile_soil inputfile_profile_soil
% inputfile_crop inputfile_leek inputfile_climate 

%% Read in inputparameters, depending on scenario_number and name of inputfiles.   
[immobile,simsol,sim_nitro,simtemp,simplant,isucr,idvs,column_depth,compartiments_number,t,tmax,dt_start,dt_min,...
    dt_max,maxiter,err_tol,plant_date,crop_type,harvest_date,...
    plant_uptake_param,miner_param,nsol,PLM2] = read_inputfile(scenario_number, inputfile_soil);
    
[soil_parameters, solute_param1, solute_param2, om_param] =  read_input_profile(compartiments_number, inputfile_profile_soil, scenario_number); 
    om_param(:,end) = om_param(:,end);

%% Read in climatic data as defined in inputfile_climate to run crop
%simulations
climate = importdata(inputfile_climate); 

%% Extra definitions based on general inputdata
units={'cm' 'day' 'g' '°C'};  %current version only supports this units
dx = column_depth/compartiments_number;	        %space increment (L)
dx_inter(1) = dx/2;								%array containing the distance in between the nodes (L)
dx_inter(2:compartiments_number) = dx;
dx_inter(compartiments_number+1) = dx/2;
print_time=[t t+0.1 t+0.2 t+0.3 t+0.4 t+0.5 t+0.6 t+0.7 t+0.8 t+0.9 t+1];  %0 100 150 200 ];
print_node = [1:1:compartiments_number];  
ncbot=compartiments_number; %MAMS
ncomp = compartiments_number;
ncs = compartiments_number;


%% Initialization
tnode=[];wc_node=[];ph_node=[];SOLNODE=[];top_inf=[];trans=[];
wat_flxs=[];bot_flux=[];top_flux=[];front=[];bc=[];snode=[];evap=[];sink=[];
seep=0;case_breaking=0;runoff=0;
bctop_changed=0;boco_top_type=0;boco_top=0;
bcbot_changed=0;boco_bot_type=0;boco_bot=0;
cum_top_flxs = 0;cum_bot_flxs = 0;cum_sink_wat = 0;
ponded=0;pond=0;pond_from=9999;
phsurf=0;phbot=0;flxar=0;flxsbot=0;
runoff=0;first_time_bc=1;dt_changed_bc=0;next_dt_new_bc=0; iter=0;
flxsa=[]; flxsah=[]; pvela=[]; pvelah=[]; pvelo=[]; pveloh=[];wcio=[];
wciob=[]; wcma=[]; wcmah=[];wcmob=[]; wcmo=[]; wco=[]; wcob=[];
cm=[];cim=[];csol=[];wat_flxsa=[];diffus=[];tot_upt=0;
acsolmo=[];acsolio =[];cum_evap=[];cum_infiltr=[];bot_inf=[];potential_surface_flux=[];
cum_potential_surface_flux=0;potential_transp=[];cum_evap=0;
cum_infiltr=0;cum_pot_transp=0;cum_trans=0;applic_boolean =0;reservoir=0;dt_sol_count=0;first_time=1;
time_uptake=[];root_length_time = [];root_density_time2=[];root_density_time3=[];water_storage=[];
flxsa=[]; flxsah=[]; pvela=[]; pvelah=[]; pvelo=[]; pveloh=[];wcio=[];
wciob=[]; wcma=[]; wcmah=[];wcmob=[]; wcmo=[]; wco=[]; wcob=[];
cm=[];cim=[];csol=[];wat_flxsa=[];diffus=[];tot_upt=0;
acsolmo=[];uptakem=[];uptakei=[];minerm=[];mineri=[];nitrifm=[];nitrifi=[];denitm=[];
deniti=[];hydro_uream=[];hydro_ureai=[];volatm=[];volati=[];decnorg=0;deccorg=0;carbman=0;carblit=0;rnitman=0;rnitlit=0;
tflnorg=0;tflcorg=0;cco2o=0;cum_nit_sink=[];
tic;om_balance=[]; N_reaction_balance=zeros(1,10);;uptake_matrix=zeros(1,6);drza=0;uptakem = zeros(1,3);root_dens_time = [];
initsol=0;Nuptake_time=[];N_reaction_balance=[0 0 0 0 0 0 0 0 0 0];
cropdata = [];joachimtable=[0 0 0];joachimtable2=[0 0];

%disp('  time [min]  t. step[min]  nber of it.  real time [sec]');

%% Read Boundary Conditions
%Water
%[BOUNDARY_CONDITIONS_MATRIX,hs,phsa,arel,brel,initial_gwl]=In_Boundary_conditions(wat_bottom_BC_type,theta_table);
wat_bottom_BC_type = [4 0];
irri=[733774 0];

%Solute
if simsol
    %Defined in in_solute_boundary, called by solve_solute
    [conc,solute_applic,solboco_top_type,solboco_input_type]=In_solute_boundary(t,nsol);
    [ddepsol, wdepsol, fsol, conirsol, ppdepth]= in_solute_upper(nsol,dx,scenario_number);
In_solute_boundary.solute_applic = solute_applic;
In_solute_boundary.solboco_top_type = solboco_top_type;
In_solute_boundary.solboco_input_type = solboco_input_type ;
end

%Organic matter application
[om_appl]= in_om_applic() ;

%Boundary conditions temperature
%defined in in_boundary_cond_temp, called by solve_temperature_JV
temp_BC=[2 NaN 2 NaN];

%% Read initial conditions
%read initial pressure head
ph = In_initial_pressure_head(compartiments_number,soil_parameters,dx_inter);
if simsol
    
% Read initial solute
[cs,tcsolo_ini] = In_initial_solute(nsol,ncs);
tcsolo = tcsolo_ini;
    
%Read initial organic mater
[soil_om]=In_initial_om(soil_parameters,ncs);
end

% Read initial temperature
[temp,temptopp]=In_initial_temperature(ncomp);

cropinit.DWLeaf = [0.28 1.96 0.28 1.96 0.28 1.96];
cropinit.DWStem = [0.15 0.35 0.15 0.35 0.15 0.35];
cropinit.DWRoot = [0 0 0 0 0 0];
cropinit.SLA = [100 100 100 100 100 100];
cropinit.NLeaf =[5.4 4.5 5.4 4.5 5.4 4.5];
cropinit.NStem =[5.4 3 5.4 3 5.4 3];

cd ../

execute_WAVE