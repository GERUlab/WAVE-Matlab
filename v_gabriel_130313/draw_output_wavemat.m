%File to plot the most import output of WAVEMAT

%% Time step information
draw_time(bc)

%% Soil properties
draw_soilprop(soil_parameters, phsa,units)

%% Output of water
%Fluxes
draw_fluxes(tnode, potential_surface_flux,potential_transp, top_inf,trans,bot_inf,...
    cum_potential_surface_flux,cum_pot_transp,top_flux, cum_trans ,bot_flux, cum_infiltr, cum_evap,units, simplant,water_storage)

%water content every timestep at speficied nodes(defined via print_node)
%Water profile every node at specified times (defined via print_time)
draw_output_water(wc_profile,ph_profile,tnode,ph_node,wc_node,print_node,...
    print_time,column_depth,dx,units,simplant)

%% Output of temperature
draw_temp(temp_profile,tnode,temp_node,compartiments_number,dx,...
    print_time,print_node,simtemp,column_depth,units);

%% Ouput of solute
%Solute balance
draw_sol_balance(solute_balance,units,nsol,N_reaction_balance,sim_nitro)

%Solute concentrations
draw_cm(cm_profile,cm_node,tnode,compartiments_number,dx,print_time,...
    print_node,nsol,simsol,column_depth,units);