
crop_climate = convert_climate_to_matlab(climate);
[ET0_time ,ET0_cm_per_day] =calc_ET02(climate);
[theta_table, temp_table] =convert_climate_to_daily(crop_climate,irri,ET0_cm_per_day);

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
%disp('  time [min]  t. step[min]  nber of it.  real time [sec]');

%% Read Boundary Conditions
%Water
[BOUNDARY_CONDITIONS_MATRIX,hs,phsa,arel,brel,initial_gwl]=In_Boundary_conditions(wat_bottom_BC_type,theta_table);

%Solute
if simsol
    %Defined in in_solute_boundary, called by solve_solute
%    [conc,solute_applic,solboco_top_type,solboco_input_type]=In_solute_boundary(t,nsol);
%    [ddepsol, wdepsol, fsol, conirsol, ppdepth]= in_solute_upper(nsol,dx,scenario_number);
In_solute_boundary.solute_applic = solute_applic;
In_solute_boundary.solboco_top_type = solboco_top_type;
In_solute_boundary.solboco_input_type = solboco_input_type ;
end

%Organic matter application
%[om_appl]= in_om_applic() ;

%Boundary conditions temperature
%defined in in_boundary_cond_temp, called by solve_temperature_JV

%% Read initial conditions
%read initial pressure head
%ph = In_initial_pressure_head(compartiments_number,soil_parameters,dx_inter);
if simsol
    
% Read initial solute
%[cs,tcsolo_ini] = In_initial_solute(nsol,ncs);
tcsolo = tcsolo_ini;
    
%Read initial organic mater
%[soil_om]=In_initial_om(soil_parameters,ncs);
end

% Read initial temperature
%[temp,temptopp]=In_initial_temperature(ncomp);


%% Create initial concentration, waterprofile, temperatureprofile
if simsol ==0
    nsol = 1;
end
[wc_profile, wc_node, cm_profile,temp_profile,cm_node,temp_node,om_profile, om_node,tc_profile,tc_node,solute_balance,ph_node,ph_profile]=...
    create_profile(nsol,print_time,compartiments_number,t);
tflsol=zeros(nsol,1);rleasa=zeros(nsol,1);tcsink=zeros(nsol,1);tsoli=zeros(nsol,1);
wcp = zeros(ncomp,1); cum_sink_wat2=zeros(1,ncomp);

tflsol=zeros(nsol,1);rleasa=zeros(nsol,1);tcsink=zeros(nsol,1);tsoli=zeros(nsol,1);
wcp = zeros(ncomp,1);
ptup = zeros(3,1); cupt = zeros(ncs,3);
ptmin = zeros(3,1); cmin = zeros(ncs,3);
pthyd = 0; chyd = zeros(ncs-1,1);
ptnit = 0; cnit = zeros(ncs-1,1);
ptvol = 0; cvol = zeros(ncs-1,1);
ptden = 0; cden = zeros(ncs-1,1);
ptsnorg = 0; %%initial sink is zero
ptscorg = 0; %%initial sink is zero
if sim_nitro
    tcorgs = sum((soil_om(:,1)) + soil_om(:,3) + soil_om(:,5));
    tnorgs = sum((soil_om(:,2)) + soil_om(:,4) + soil_om(:,6));
end
cropresults =[];cropdata =0;root_depth_time=[];
if simplant ==1
for i=1:size(plant_date,2)
output_crop{i}=[];
output_cropNup{i}=[];
end
else ncrop = NaN; 
end

%% Initial calculations
check_ini(soil_parameters,dt_start,dt_min,dt_max);
dt=dt_start;
stock_initial = calc_stock(ph,dx,soil_parameters,dt,0);
stock_max=sum(moist_ret(zeros(1,compartiments_number),soil_parameters,dt,0)*dx);
%u=figure;

%% Start simulation
while t < tmax

if round(t) ==t
if simplant
    ncrop = find_cropnumber (t,simplant,plant_date, harvest_date);
        if isnan(ncrop)==0 ;%Plant simulation and on the field
            flag_double_sim = 1;
        else  %Plant simulation, but currently not on the field
            flag_double_sim =0;
        end
    else 
        flag_double_sim = 0 ; %No plant simulation, never plant simulation
end

%% Use Booleans to run crop model if needed   
 if flag_double_sim ==1
%Initialize crop model if crop is planted     
if t== plant_date(ncrop)      %Check if DAP = 0 
    %==> Initialize
    TeSum = 0; TotDMLeaf = 0; TotDMStem = 0;TotDMCurd = 0;TotDMRoot = 0;TotDMShoot = 0;DVS = 0.004;TotDMGreenLeaf = 0;
    TotDMDeadLeaf = 0;TotNperLeaf = 0;TotNperStem = 0;CumEvaporation=0;CumTranspiration=0;RLD = zeros(100,1) ;   
    if crop_type(ncrop)==1
        CropParamsfile = importdata(inputfile_crop,',',1);
        CropParams = CropParamsfile.data(1,:);
        %Initialisation
        Holding1(1,:) = [0,0,0,0,0,0,0,0,0];   
   elseif crop_type(ncrop)==2
        CropParamsfile = importdata(inputfile_leek,',',1);
        CropParams = CropParamsfile.data(1,:); 
        %Initialisation
         Holding2(1,:) = [0,0,0,0,0,0,0,0,0];
   end
end
%Run crop simulation for one day
 if crop_type(ncrop)==1    
    [crop_result, RootResDay,crop_result2] =... 
            run_crop_modeles3(CropParams,output_crop,simplant, plant_date,...
            harvest_date,t, Nuptake_time,output_cropNup,PLM2(ncrop),cropinit,crop_climate);
     %Save output
     output_crop{ncrop}(end+1,:) =crop_result;
     Holding1 = [Holding1;RootResDay];
     output_cropNup{ncrop}(end+1,:) =crop_result2;
     drz =output_crop{ncrop}(end,25);% TotRootDepth1;
    fraction_plant = calc_fraction_plant(crop_type,ncrop,CropParams, crop_result);
    index = find(ET0_time== t);
    LAI = output_crop{ncrop}(end,22);
     epa = ET0_cm_per_day(index)*(1-exp(-0.6*LAI));
     esa = ET0_cm_per_day(index)-epa;
     NDemand = output_crop{ncrop}(end,34);
              
 elseif crop_type(ncrop)==2
    [crop_result, RootResDay,crop_result2] =... 
                 run_crop_model_leek2(CropParams,output_crop,simplant,...
                 plant_date, harvest_date,t, Nuptake_time,output_cropNup,PLM2(ncrop),cropinit,crop_climate);
        % Save output
    output_crop{ncrop}(end+1,:) =crop_result;
    Holding2 = [Holding2;RootResDay];
    output_cropNup{ncrop}(end+1,:) =crop_result2;
    drz =output_crop{ncrop}(end,26);% TotRootDepth1; 
    fraction_plant = calc_fraction_plant(crop_type,ncrop,CropParams, crop_result);
    index = find(ET0_time== t);
    LAI = output_crop{ncrop}(end,19);
     epa = ET0_cm_per_day(index)*(1-exp(-0.6*LAI));
     esa = ET0_cm_per_day(index)-epa;
    NDemand = output_cropNup{ncrop}(end,13);
 end  

tot_upt = 0; %reset total supply by soil to zero.
        
else %Means no crop simulation ==> all crop related parameter to zero
       index = find(ET0_time== t);
        esa = ET0_cm_per_day(index);
if isnan(esa)
    esa = 0;
end

        drz = 0; epa = 0;NDemand = 0;RLengthLa = 0;DATE = 0;tot_upt = 0;
end
 
if isnan(ncrop)==0
          crop_type_act = crop_type(ncrop);
          DAP = t-plant_date(ncrop);
          [rdens, drza] = calc_rdens2(DAP,drz,ncrop,crop_type_act);
else
    if simplant ==0
          rdens = NaN;
          drza = NaN;
          plant_date = NaN;
          harvest_date = NaN;
    else
          rdens = NaN;
          drza = NaN;        
    end
end
end


    % Evaluate bc at time t and t+dt, chooses right bc, checks if bc 
    % changed during previous iteration
    var_changed =[bctop_changed,bcbot_changed,dt_changed_bc,next_dt_new_bc,first_time_bc];
    
    [hs,phsa,flxar,flxa1,boco_top_type,boco_top,boco_bot_type,boco_bot,pond,...
        bctop_changed,bcbot_changed,first_time_bc,dt_changed_bc,next_dt_new_bc,phsurf,esa,...
        BOUNDARY_CONDITIONS_MATRIX,ncbot]=...
        soil_boundary_conditions_new(var_changed,t,ph,dt,phsurf,flxar,boco_top_type,...
        boco_top,boco_bot_type,boco_bot,pond,runoff,simplant,dt_min,dt_start,...
        BOUNDARY_CONDITIONS_MATRIX,esa,hs,phsa,ncbot,dx_inter,ncs,arel,brel,plant_date,harvest_date,wat_bottom_BC_type,theta_table);
    
    % Calculate pressure head(ph)
	[ph,wat_flxs,dt,iter,runoff,no_conv,WC,bctop_changed,bcbot_changed,pond,...
        phsurf,flxar,boco_top_type,boco_bot_type,rtex,EPRA,boco_top] =...
        solve_flow(ph,t,dt,compartiments_number,iter,boco_top_type,boco_top,boco_bot_type,...
        boco_bot,pond,flxa1,dt_start,dx,dx_inter,stock_max,seep,soil_parameters,ponded,pond_from,...
        hs,phsa,maxiter,dt_min,phsurf,phbot,flxsbot,flxar,bctop_changed,bcbot_changed,err_tol,simplant,...
        units,plant_date,harvest_date,ncbot,drz,epa);      
    
      %Calculate temperature
      if sim_temp
        [temp, temptopp] = solve_temperature_JV(dt, dx, ncomp, t, WC, wcp, temp, soil_parameters,temptopp,tmax,units,solute_param1,temp_BC,temp_table);

      end
      
      %Calculate solute
     if simsol
         if sim_temp ==0
          temp = ones(ncs,1).*7;
         end
         
         if sim_nitro ==0
              soil_om=ones(ncs,1)*[0,0,0,0,0,0];
              miner_param= [0 0 0];
              om_param= [0 0 0]; om_appl = 0; PLM2 = 0;
         end
        [wat_flxsa, wat_flxsah, pvela, pvelah, pveloh,pvelo, wcio, wciob, wcma, wcmah,wcmob, wcmo, wco, wcob,diffus,diffush]=...
         wat_sol(t,wat_flxs, wat_flxsa, immobile, ncs,solute_param1,solute_param2,pvela,WC, wcio, wcma,wcmo, wco,pvelo,diffus,nsol,dx,initsol,solboco_top_type);
         
[cs,csp,cm,cim,cmp,cimp,csol,soil_om,sol_sinkm,sol_sinki,tot_upt,acsolo,acsolio,acsolmo,sflxsa,...
            ddepsol, wdepsol, fsol, conirsol, ppdepth,reservoir,first_time,uptake_matrix,uptakem,...
            solute_applic,uptakei,minerm,mineri,nitrifm,nitrifi,denitm,deniti,hydro_uream,hydro_ureai,...
            volatm,volati,decnorg,deccorg,carbman,carblit,rnitman,rnitlit,om_appl,und,unc,rdens]=...
        solve_solute(t,dt,tmax,temp,dx,nsol,ncs,WC,rtex,cm,cim,csol,cs,acsolmo,acsolio,wat_flxs,ph,...
                    wat_flxsa, wat_flxsah,wcio, wciob, wcma, wcmah, wcmo, wcmob,wco,wcob,pvela,pvelah,...
                    pveloh,soil_om,immobile,soil_parameters,solute_param1,solute_param2,...
                    om_param,miner_param, plant_uptake_param, tot_upt,sim_nitro,tcsolo_ini,diffush,...
                    plant_date,harvest_date,isucr,simplant,ddepsol, wdepsol, fsol, conirsol,...
                    ppdepth,reservoir,tcsolo,first_time,om_appl,drz,NDemand,PLM2,1,initsol,...
                    solboco_input_type,RLengthLa,0,rdens,ncrop,In_solute_boundary,theta_table);

        [cberr_sol,tflsol,rleasa,tcsink,tsoli,tcsolo_ini,dsol,tcsolo]=sol_intgr(t,dt,dx,nsol,ncs,cm,cim,cmp,cimp,csol,wcmo,wcio,wcob,wco,solute_param1,...
            solute_param2,sol_sinki,sol_sinkm,acsolmo,acsolio,acsolo,sflxsa,tflsol,rleasa,tsoli,tcsink,tcsolo_ini,immobile,initsol);
        
                  if sim_nitro
              [ptsnorg, ptscorg ptup, cupt, pthyd, chyd, ptnit, cnit, ptvol, cvol, ptden, cden,...
                dnorg, dcorg, cberr_norg, cberr_corg,tflnorg,tflcorg,ptmin,cmin,om_balance] = ...
                nit_intg(ptsnorg,ptscorg, ptup, cupt, ptmin, cmin, pthyd, chyd, ptnit, cnit, ptvol,...
                 cvol, ptden, cden, soil_om,tcorgs,tnorgs,carbman,carblit,rnitman,rnitlit,decnorg,deccorg,...
                 uptakem,uptakei,minerm,mineri,nitrifm,nitrifi,denitm,deniti,hydro_uream,hydro_ureai,...
                volatm,volati,tflnorg,tflcorg,dx,dt,om_balance);
                  end

     end
   
     % Calculate output data
[delta_stock,mass_balance_error, cum_top_flxs, cum_bot_flxs, cum_sink_wat, cum_sink_wat2,sink_prof] = ...
   wat_intgr(ph,wat_flxs,cum_top_flxs ,cum_bot_flxs,cum_sink_wat, cum_sink_wat2, stock_initial,dt,dx,...
   compartiments_number,soil_parameters,rtex);
	% Display information
   	disp([datestr(t)]);%,dt,iter,boco_top_type,boco_top]);
    % Store data for model output 
[wc_profile,tnode,wc_node,ph_node,bot_flux,top_flux,top_inf,snode,evap,trans,sink,ph_profile,bot_inf, potential_surface_flux, cum_potential_surface_flux, potential_transp, cum_evap, cum_infiltr,cum_pot_transp,cum_trans,water_storage] = store_data(t,ph,...
   wc_profile,tnode,ph_node,wc_node,print_time,print_node,cum_bot_flxs,cum_top_flxs,...
   bot_flux,top_flux,no_conv,WC,wat_flxs,top_inf,esa,cum_sink_wat,snode,evap,trans,sink_prof,sink,ph_profile,bot_inf, potential_surface_flux, cum_potential_surface_flux, potential_transp, cum_evap, cum_infiltr,flxar,dt,epa,cum_pot_transp,cum_trans,cum_sink_wat2,simplant,water_storage,dx,theta_table);

if simsol
    %store Concentration
[cm_profile,cm_node]=store_CM(print_time,...
    print_node,cm,no_conv,t,cm_profile,cm_node,wc_profile,nsol,simsol);
    %Store organic matter pools
[om_profile,om_node]=store_om(print_time,...
    print_node,soil_om,no_conv,t,om_profile, om_node,sim_nitro);
%Store absolute values of solute
[tc_profile,tc_node]=store_tc(print_time,...
   print_node,tcsolo,no_conv,t,tc_profile,tc_node,nsol,simsol);
%Store Mass balance information for the solutes.
solute_balance = store_solute_balance(cberr_sol,tflsol,rleasa,tcsink,tsoli,tcsolo_ini,dsol,tcsolo,solute_balance,t,nsol,dt);
    
    if sim_nitro
        [cum_nit_sink] = store_nit(ptsnorg, ptscorg,ptup, cupt, pthyd, chyd, ptnit, cnit, ptvol, cvol, ptden, cden,...
            dnorg, dcorg, cberr_norg, cberr_corg,tflnorg,tflcorg,ptmin,cmin, cum_nit_sink);
    end
    
end

if simplant
    %Store uptake terms
    time_uptake(end+1,1:6) = [uptake_matrix];
    %Store root information
        root_length_time (end+1) = [drza];
        root_dens_time(end+1,:) = rdens;
        root_density_time2 (end+1,:) = uptakem(:,2)';
        root_density_time3 (end+1,:) = uptakem(:,3)';
end


if sim_temp
    [temp_node, temp_profile]=store_T(print_time,...
    print_node,no_conv,t, temp, temp_profile, temp_node);
end

% Calculation of delta t
   	t=t+dt; 
    [dt,dt_changed_bc]=calc_dt_new(t,dt,iter,dt_min,dt_max,BOUNDARY_CONDITIONS_MATRIX,rtex, wat_flxs,dx,ncs);
    if simsol
    [dt_sol,applic_boolean,dt_sol_count] = calc_dt_sol(t,nsol,dt, applic_boolean,dt_min, units,cs,csp,dt_sol_count,dx,reservoir,fsol, conirsol,wdepsol,ddepsol,In_solute_boundary, theta_table);
    else
    dt_sol = dt;
    end
    dt = min(dt, dt_sol);
   
% Evaluate bc
    bc(end+1,:)=[t,boco_top_type,boco_top,boco_bot_type,boco_bot,pond,runoff,flxa1,dt,iter];
   
    % Break when no convergence is reached
    if no_conv
        disp('No convergence');
        break
    end
    initsol=1;
%% Output of soil model to crop model
Nuptake_time(end+1,1) = DATE;
if tot_upt ==0 
Nuptake_time(end,2)=0 ;   
else
Nuptake_time(end,2) = tot_upt/PLM2(ncrop)*10^4; %Actual uptake (gram/plant)
end
Nuptake_time(end,3) = NDemand;  

N_reaction_balance(end+1,:) = [ptup(1) ptup(2) ptup(3) ptmin(1) ptmin(2) ptmin(3) pthyd ptnit ptvol ptden];
end

try
save([handles.dirpath handles.projectname '\output'])
catch
    save output
end
disp('model finished')



