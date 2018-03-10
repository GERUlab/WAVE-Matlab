%Compare the different ways of calculating thermal conductivity

%Read in data as done in the main program
[t,tmax,compartiments_number,dx_inter,dx, dt_start,dt_min,print_time,print_node,...
        maxiter,dt_max,column_depth,err_tol,simplant,simsol,isucr,idvs,sim_nitro,...
      immobile,units,plant_date,harvest_date,simtemp,crop_type] = In_general_data();
  
ncomp =  compartiments_number;
soil_parameters=In_soil_parameter_data(ncomp);
[nsol,solute_param1,solute_param2]=in_solute_parameter(ncomp,0);

%Define a range of water contents
WC = linspace(0,0.45);
wcp = WC';

%Calculate thermal conductivities in three ways
%(input parameters come from In_temperature_param,
%In_temperature_param2,In_temperature_param3)
[allam, beta] = calc_thermal_conduc(units,ncomp,soil_parameters, solute_param1,WC,wcp);
[allam2,beta2] = calc_thermal_conduc2(units,ncomp,soil_parameters, solute_param1,WC,wcp);
[allam3, beta3] = calc_thermal_conduc3(units,ncomp,soil_parameters,...
solute_param1,WC,wcp);
[allam4, beta4] = calc_thermal_conduc4(units,ncomp,soil_parameters, solute_param1,WC,wcp);

WC(ncomp+1) = WC(ncomp);

%Make the plot
plot(WC, allam,WC, allam2,WC, allam3, WC, allam4*10^10)
xlabel('water content')
ylabel('thermal conductivity (J cm-1 day-1 C-1)')
legend('De Vries', 'Campbell', 'Chung', 'De Vries Daisy')
