%Create inputfile with one point
immobile = 0;
simsol = 1; %boolean for solute simulation
sim_nitro = 1; %boolean for nitrogen simulation
sim_temp = 1; %boolean for temperature simulation
simplant = 1; ;% boolean for plant simulation
isucr = 1; %boolean for use of sucros (not used as for now)
idvs = 0; %boolean for use of development stages (not used as for now)

column_depth = 85; %depth of the column (in cm)
compartiments_number = 85%number of compartments

t = 733774; %initial time (day)skype
tmax =733774+365*3-1 ;%maximum time (day)
dt_start = 10^-4; %start time step (day)
dt_min = 10^-6; %minimum time step (day)
dt_max = 0.5; %maximum time step (day)
maxiter = 1000; %maximum number of iterations
err_tol = 10^-4; %error tolerance

rorad = 0.0224; %mean root radius
rd0 = 10^-5;
g = 1;
rnmaxp = 0.0024; %maximum nitrogen uptake
rdens0 =0.01; %root density near surface
alfardens =  0; %exponential factor for root density

r0 = 10; % biomass C-N ratio 
fe =0.3;%turn over efficiency
fh =0.4; %humification constant

plant_date = [61 175 441 540 799 898]+733774; %plant dates of each crop
crop_type = [1 2 1 2 1 2]; %crop type for each crop (1= cauliflower, 2=leek)
harvest_date = [160 308 530 686 873 1027]+733774; %harvest date for each crop
PLM2 = [2.85 15.58 2.85 15.58 2.85 15.58]; %plant density (plants/m²)

nsol = 3; %number of solutes

%Combine certain parameters from above in one matrix
plant_uptake_param = [rorad, rd0, g, rnmaxp, rdens0, alfardens];
miner_param = [r0, fe, fh];

outputmatrix = [immobile,simsol,sim_nitro,sim_temp,simplant,isucr,idvs,column_depth,compartiments_number,t,tmax,dt_start,dt_min,...
    dt_max,maxiter,err_tol,...
    plant_uptake_param,miner_param,nsol,plant_date, crop_type,harvest_date,PLM2]

%writ
xlswrite('outputmatrix.xls', outputmatrix)