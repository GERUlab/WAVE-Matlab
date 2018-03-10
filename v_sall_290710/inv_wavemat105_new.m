function [tnode wc_node] = inv_wavemat105_new(soil_parameters)

% INV_WAVEMAT105_new Main programme 
% Programme solving the inversion of  transient water flow equation 
%
%-------------------------------
% M. Vanclooster 17/2/2000
% modified by Vanclooster 5/4/08

% Initialization
tnode=[];wc_node=[];ph_node=[];SOLNODE=[];top_inf=[];trans=[];
wat_flxs=[];bot_flux=[];top_flux=[];front=[];bc=[];snode=[];evap=[];sink=[];
seep=0;case_breaking=0;runoff=0;
bctop_changed=0;boco_top_type=0;boco_top=0;
bcbot_changed=0;boco_bot_type=0;boco_bot=0;
cum_top_flxs = 0;cum_bot_flxs = 0;cum_sink_wat = 0;
ponded=0;pond=0;pond_from=9999;
phsurf=0;phbot=0;flxar=0;flxsbot=0;
runoff=0;first_time_bc=1;dt_changed_bc=0;next_dt_new_bc=0; iter=0;
tic;

% Read general data 
[t,tmax,compartiments_number,dx_inter,dx, dt_start,dt_min,print_time,...
    print_node,wc_profile,maxiter,dt_max,column_depth,err_tol,simplant,...
    simsol,isucr,idvs] = In_general_data;

% Read boundary conditions 
[BOUNDARY_CONDITIONS_MATRIX,hs,phsa]=In_Boundary_conditions;
[AKc,EPA,esa,DRZ,ARER,BRER]=In_ETsplit(t,simplant);

% Read initial conditions
ph = In_initial_pressure_head(compartiments_number,soil_parameters);

%Initial calculations 
check_ini(soil_parameters,dt_start,dt_min,dt_max);
dt=dt_start;
stock_initial = calc_stock(ph,dx,soil_parameters,dt,0);
stock_max=sum(moist_ret(zeros(1,compartiments_number),soil_parameters,dt,0)*dx);

%Start simulation 
while t < tmax
    % Evaluate bc at time t and t+dt, chooses right bc, checks if bc 
    % changed during previous iteration
    var_changed =[bctop_changed,bcbot_changed,dt_changed_bc,next_dt_new_bc,first_time_bc];
    
    [hs,phsa,flxar,flxa1,boco_top_type,boco_top,boco_bot_type,boco_bot,pond,...
        bctop_changed,bcbot_changed,first_time_bc,dt_changed_bc,next_dt_new_bc,phsurf,esa,...
        BOUNDARY_CONDITIONS_MATRIX]=...
        soil_boundary_conditions_new(var_changed,t,ph,dt,phsurf,flxar,boco_top_type,...
        boco_top,boco_bot_type,boco_bot,pond,runoff,simplant,dt_min,dt_start,...
        BOUNDARY_CONDITIONS_MATRIX,esa,hs,phsa);
    
    % Calculate ph
	[ph,wat_flxs,dt,iter,runoff,no_conv,WC,bctop_changed,bcbot_changed,pond,...
        phsurf,flxar,boco_top_type,boco_bot_type,rtex,EPRA] =...
        solve_flow(ph,t,dt,compartiments_number,iter,boco_top_type,boco_top,boco_bot_type,...
        boco_bot,pond,flxa1,dt_start,dx,dx_inter,stock_max,seep,soil_parameters,ponded,pond_from,...
        hs,phsa,maxiter,dt_min,phsurf,phbot,flxsbot,flxar,bctop_changed,bcbot_changed,err_tol,simplant);
    
    % Calculate output data
    [delta_stock,mass_balance_error, cum_top_flxs, cum_bot_flxs, cum_sink_wat,sink_prof] = ...
   	wat_intgr(ph,wat_flxs,cum_top_flxs ,cum_bot_flxs,cum_sink_wat, stock_initial,dt,dx,...
   	compartiments_number,soil_parameters,rtex);
    % Store data for model output 
    [wc_profile,tnode,wc_node,ph_node,bot_flux,top_flux,top_inf,snode,evap,trans,sink] = store_data(t,ph,wc_profile,tnode,ph_node,...
        wc_node,print_time,print_node,cum_bot_flxs,cum_top_flxs,bot_flux,top_flux,no_conv,WC,wat_flxs,top_inf,...
        esa,cum_sink_wat,snode,evap,trans,sink_prof,sink);  
    % Calculation of delta t
   	t=t+dt;
    [dt,dt_changed_bc] = calc_dt_new(t,dt,iter,dt_min,dt_max,BOUNDARY_CONDITIONS_MATRIX);
    % Evaluate bc
  	bc(end+1,:)=[t,boco_top_type,boco_top,boco_bot_type,boco_bot,pond,runoff,flxa1];
    % Break when no convergence is reached
    if no_conv
        disp('No convergence');
        break
    end
end