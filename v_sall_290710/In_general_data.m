function [t,tmax,compartiments_number,dx_inter,dx, dt_start,dt_min,print_time,print_node,...
      wc_profile,maxiter,dt_max,column_depth,err_tol,simplant,simsol,isucr,idvs,sim_nitro,...
      immobile,units,plant_date,harvest_date,wc_change_tol] = In_general_data();

%General Input for the simulation
%
% IN:
%	none
% OUT:
%	t,tmax,compartiments_number,dx_inter,dx,
%  	dt_start,dt_min,maxiter,dt_max,column_depth
% plant_date, haverst_date,
% CALLS
%	none
%%     Specify the Units information 
%Column 1= geometry unit (m, cm, mm)
%Column 2 = time unit (day, h, min, s)
%Column 3 = mass unit (kg, g, mg)
%Column 4 = Temperature unit (°C)
units={'cm' 'day' 'g' '°C'};

%	Soil column informations
column_depth = 100;%origine								%depth of soil profile 
compartiments_number =100;						%number of (unsaturated) compartments
dx = column_depth/compartiments_number;	        %space increment 
dx_inter(1) = dx/2;								%array containing the distance in between the nodes (cm)
dx_inter(2:compartiments_number) = dx;
dx_inter(compartiments_number+1) = dx/2;


%	Time information
t = 0;					%Initial time 
tmax = 50; %%%316;			%Final time (min)
dt_start = 0.001;			%Initial time step 
dt_min = 1E-10;			%Minimum time step 
dt_max = 0.1;			%Maximum time step)%%%dt_max=0.01 convient pour reel

%	Iterations information
maxiter = 1000;			%maxinum number of Newton-Raphson iterations
err_tol=0.01;
wc_change_tol=0.002;


% Initialization for model output

print_time=[5];  
print_node =[20 40 60 80];                         %origine      
wc_profile=ones(length(print_time),compartiments_number);

% Simulations with (simplant=1) or without plants (simplant=0)
simplant=0;
%If simplant=1, do you want to use Sucros for crop growth modelling? If "yes" iscur=1, if "no" iscur=0
isucr=0;
%Do you want to use DVS units for crop stages? If "yes", idvs=1, if "no" idvs=0
idvs=0;
%If sipmplant specify the plant and harverst time, if not let them equal to 0  
plant_date=0;
harvest_date=76;

%Simulations of solute?
simsol=1;
% simulation of nitrogen?
sim_nitro=0;
% mmobile/immobile water concept
immobile=0; 