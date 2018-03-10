close all;clear all
% WAVEMAT105 Main programme 
% Programme solving the transient water flow equation 
%
% IN:
%	general_data (t,tmax), soil_parameter_data,initial_pressure_head
%	dx, compartiments_number, dt_start
%
% OUT:
%	wc_profile:water content profiles at different times choosen by the user in print_time (in this file).
%	wc_node: WC evolution in time at different depths choosen by the user in print_node.
%	ph_node: ph evol. in time at different depths choosen by the user in print_node [L][.
%	bot_flux: cumulative bottom flux leaving the soil profile [Lï¿½].
%
% CALL:
%	boundary_conditions_data,solve_flow, wat_intgr, store_data, draw_output
%	moist_ret, calc_stock,general_data,initial_pressure_head
% for the solutes:
%	SOL_PARAMETERS,initial_solute_concentration,solboco,calc_wat_sol,
%	sol_sink
%
%-------------------------------
% M. Vanclooster 17/2/2000
% modified by Javaux, 14/05/00
% modified by Vanclooster 4/10/07
%modified by Sall 12/04/09

% Initialization
tnode=[];wc_node=[];ph_node=[];SOLNODE=[];top_inf=[];trans=[];
wat_flxs=[];bot_flux=[];top_flux=[];front=[];bc=[];snode=[];evap=[];sink=[];cumul_evap=[];sink_wat=0;
seep=0;case_breaking=0;runoff=0;
bctop_changed=0;boco_top_type=0;boco_top=0;
bcbot_changed=0;boco_bot_type=0;boco_bot=0;
cum_top_flxs = 0;cum_bot_flxs = 0;cum_sink_wat = 0;cum_evap=0;cum_inf=0;
ponded=0;pond=0;pond_from=9999;
phsurf=0;phbot=0;flxar=0;flxsbot=0;
runoff=0;first_time_bc=1;dt_changed_bc=0;next_dt_new_bc=0; iter=0;
decnorg=0;deccorg=0;carbman=0;carblit=0;rnitman=0;rnitlit=0;
tflnorg=0;tflcorg=0;cco2o=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%MAMADOU
flxsa=[]; flxsah=[]; pvela=[]; pvelah=[]; pvelo=[]; pveloh=[];wcio=[]; 
wciob=[]; wcma=[]; wcmah=[];wcmob=[]; wcmo=[]; wco=[]; wcob=[];
cm=[];cim=[];csol=[];wat_flxsa=[];diffus=[];tot_upt=0;
acsolmo=[];uptakem=[];uptakei=[];minerm=[];mineri=[];nitrifm=[];nitrifi=[];denitm=[];
deniti=[];hydro_uream=[];hydro_ureai=[];volatm=[];volati=[];
%dxn=0;arel=2;brel=3;

tic;

%disp('  time [min]  t. step[min]  nber of it.  real time [sec]');	

% Read the data input data files
% Read general data 
[t,tmax,compartiments_number,dx_inter,dx, dt_start,dt_min,print_time,print_node,...
      wc_profile,maxiter,dt_max,column_depth,err_tol,simplant,simsol,isucr,idvs,sim_nitro,...
      immobile,units,plant_date,harvest_date,wc_change_tol] = In_general_data();
ncbot=compartiments_number;%%%%%%%%%%%%%%MAMS

% Read soil parameters
soil_parameters=In_soil_parameter_data(compartiments_number);
% Read solute parameters% 
[nsol,solute_param1,solute_param2,ncs,sol_upper_bc_type]=in_solute_parameter(compartiments_number,immobile);
%Read nitrogen parameters
[om_param, miner_param, plant_uptake_param]=In_nitro_om_param(ncs);

% Read boundary conditions 
[BOUNDARY_CONDITIONS_MATRIX,hs,phsa,arel,brel,initial_gwl,tprec]=In_Boundary_conditions;
[AKc,epa,esa,drz,Smax_param, Smax_given] = In_ETsplit(t,simplant,units,plant_date,harvest_date);
%Read solute application

% Read initial conditions
ph = In_initial_pressure_head(compartiments_number,soil_parameters,dx_inter,BOUNDARY_CONDITIONS_MATRIX,initial_gwl);
% Read solute and temperature initial conditions
[cs,tcsolo_ini] = In_initial_solute(nsol,ncs);
tflsol=zeros(nsol,1);rleasa=zeros(nsol,1);tcsink=zeros(nsol,1);tsoli=zeros(nsol,1);
[rocpp,tempp,lamdap,ttopp,tbotp,wcs]=In_initial_temperature(ncs,soil_parameters);
%Read initial organic mater
[soil_om,ptup,cupt,pthyd,chyd,ptnit,cnit,ptvol,cvol,ptden,cden,ptsnorg,ptscorg,...
    tnorgs,tcorgs,ptmin,cmin]=In_initial_om(ncs);
%Read initial concentration profile
[cm_profile,temp_profile,cm_node,temp_node,bot_sol_flux]=in_cmt_profile(nsol,print_time,compartiments_number);

%Initial calculations 
check_ini(soil_parameters,dt_start,dt_min,dt_max);
dt=dt_start;
stock_initial = calc_stock(ph,dx,soil_parameters,dt,0);
stock_max=sum(moist_ret(zeros(1,compartiments_number),soil_parameters,dt,0)*dx);
u=figure;

%Start simulation 
while t < tmax
    % Evaluate bc at time t and t+dt, chooses right bc, checks if bc 
    % changed during previous iteration
    var_changed =[bctop_changed,bcbot_changed,dt_changed_bc,next_dt_new_bc,first_time_bc];
    
    [hs,phsa,flxar,flxa1,boco_top_type,boco_top,boco_bot_type,boco_bot,pond,...
        bctop_changed,bcbot_changed,first_time_bc,dt_changed_bc,next_dt_new_bc,phsurf,esa,...
        BOUNDARY_CONDITIONS_MATRIX,ncbot]=...
        soil_boundary_conditions_new(var_changed,t,ph,dt,phsurf,flxar,boco_top_type,...
        boco_top,boco_bot_type,boco_bot,pond,runoff,simplant,dt_min,dt_start,...
        BOUNDARY_CONDITIONS_MATRIX,esa,hs,phsa,ncbot,dx_inter,ncs,arel,brel,tprec);
    
    % Calculate ph
	[ph,wat_flxs,dt,iter,runoff,no_conv,WC,bctop_changed,bcbot_changed,pond,...
  phsurf,flxar,boco_top_type,boco_bot_type,rtex,EPRA,esa] =...
   solve_flow(ph,t,dt,compartiments_number,iter,boco_top_type,boco_top,boco_bot_type,...
   boco_bot,pond,flxa1,dt_start,dx,dx_inter,stock_max,seep,soil_parameters,ponded,pond_from,...
   hs,phsa,maxiter,dt_min,phsurf,phbot,flxsbot,flxar,bctop_changed,bcbot_changed,err_tol,...
   simplant,units,plant_date,harvest_date,ncbot,wat_flxs);

       
      
      %CALCULATE TEMPERATURE
      [temp,rocp,lamda,ttop,tbot]=solve_temperature(t,dt,dx,WC,wcs,tempp,ttopp,tbotp,rocpp,lamdap,ncs);
    
      %CALCULATE SOLUTE
      if simsol
          [wat_flxsa, wat_flxsah, pvela, pvelah, pveloh,pvelo, wcio, wciob, wcma, wcmah,wcmob, wcmo,...
              wco, wcob,diffus,diffush]=...
              wat_sol(t,wat_flxs, wat_flxsa, immobile, ncs,solute_param1,solute_param2,pvela,WC, wcio,...
              wcma,wcmo, wco,pvelo,diffus,nsol,dx,sol_upper_bc_type);
          
          [cs,cm,cim,cmp,cimp,csol,soil_om,sol_sinkm,sol_sinki,tot_upt,acsolo,acsolio,acsolmo,sflxsa,...
    solute_applic,uptakem,uptakei,minerm,mineri,nitrifm,nitrifi,denitm,deniti,hydro_uream,hydro_ureai,...
    volatm,volati,decnorg,deccorg,carbman,carblit,rnitman,rnitlit,cco2o]=...
    solve_solute(t,dt,tmax,temp,dx,nsol,ncs,WC,rtex,cm,cim,csol,cs,acsolmo,wat_flxs,ph,wat_flxsa,...
    wat_flxsah,wcio, wciob, wcma, wcmah, wcmo, wcmob,wco,wcob,pvela,pvelah,pveloh,soil_om,immobile,...
    soil_parameters,solute_param1,solute_param2,om_param, miner_param, plant_uptake_param,tot_upt,...
    sim_nitro,tcsolo_ini,diffush,plant_date,harvest_date,isucr,simplant,tprec,uptakem,uptakei,minerm,...
    mineri,nitrifm,nitrifi,denitm,deniti,hydro_uream,hydro_ureai,volatm,volati,decnorg,deccorg,carbman,...
    carblit,rnitman,rnitlit,cco2o,sol_upper_bc_type);
          
          [cberr_sol,tflsol,rleasa,tcsink,tsoli,tcsolo_ini,dsol,bot_sol_flux]=sol_intgr(t,dt,dx,nsol,ncs,cm,cim,cmp,cimp,csol,...
              wcmo,wcio,wcob,wco,solute_param1,solute_param2,sol_sinki,sol_sinkm,acsolmo,acsolio,acsolo,sflxsa,...
              tflsol,rleasa,tsoli,tcsink,tcsolo_ini,immobile,bot_sol_flux);
          
          if sim_nitro
              [ptsnorg, ptscorg ptup, cupt, pthyd, chyd, ptnit, cnit, ptvol, cvol, ptden, cden,...
                dnorg, dcorg, cberr_norg, cberr_corg,tflnorg,tflcorg,ptmin,cmin] = ...
                nit_intg(ptsnorg,ptscorg, ptup, cupt, ptmin, cmin, pthyd, chyd, ptnit, cnit, ptvol,...
                 cvol, ptden, cden, soil_om,tcorgs,tnorgs,carbman,carblit,rnitman,rnitlit,decnorg,deccorg,...
                 uptakem,uptakei,minerm,mineri,nitrifm,nitrifi,denitm,deniti,hydro_uream,hydro_ureai,...
                volatm,volati,tflnorg,tflcorg,dx,dt);
          end
          
      end
      %%Graph
    subplot(1,4,1);plot(ph,length(ph):-1:1);axis([min(-600,min(ph)),100,0,length(ph)]);
        title('Pressure head vs. z');xlabel('pressure head [cm]');ylabel('z [cm]');grid on;
    subplot(1,4,2);plot(WC,length(WC):-1:1);axis([min(soil_parameters(:,1)),max(soil_parameters(:,2)),0,length(ph)]);
        title('Water content vs. z');grid on;xlabel('water content [cm^3/cm^3]');ylabel('z [cm]');
    %%%% plot solute
   
     %%%% plot temperature
    % subplot(1,4,3);plot(temp,length(temp):-1:1);axis([min(temp),max(temp),0,length(temp)]);
     %title('temp vs. z');xlabel('T[°C]');ylabel('z [cm]');grid on;
     %%%%%plot nitrogen
     if simsol
     subplot(1,4,4);plot(cm,length(cm):-1:1);axis([0,1000,0,length(cm)]);
     title('cm vs. z');xlabel('cm[mg/cm3]');ylabel('z [cm]');grid on; drawnow;
     end
     % Calculate output data
    [delta_stock,mass_balance_error, cum_top_flxs, cum_bot_flxs, cum_sink_wat,sink_prof,cum_evap,cum_inf,sink_wat] = ...
   	wat_intgr(ph,wat_flxs,cum_top_flxs ,cum_bot_flxs,cum_sink_wat, stock_initial,dt,dx,...
   	compartiments_number,soil_parameters,rtex,cum_evap,cum_inf,esa,sink_wat);
	% Display information
   	disp([t,dt,iter,boco_top_type,boco_top]);
    % Store data for model output 
    [wc_profile,tnode,wc_node,ph_node,bot_flux,top_flux,top_inf,snode,evap,trans,sink,cumul_evap] = store_data(t,ph,wc_profile,tnode,ph_node,...
        wc_node,print_time,print_node,cum_bot_flxs,cum_top_flxs,bot_flux,top_flux,no_conv,WC,wat_flxs,top_inf,...
        esa,cum_sink_wat,snode,evap,trans,sink_prof,sink,ncs,cum_evap,cum_inf, cumul_evap,sink_wat);  
    %store CM and T
  [cm_profile,temp_profile,cm_node,temp_node]=store_cmt(print_time,...
   print_node,cm,temp,no_conv,t,cm_profile,temp_profile,cm_node,temp_node,wc_profile,nsol,simsol);
    % Calculation of delta t
   	t=t+dt; 
    tempp=temp; rocpp=rocp; lamdap=lamda; ttopp=ttop;
    [dt,dt_changed_bc]=calc_dt_new(t,dt,iter,dt_min,dt_max,BOUNDARY_CONDITIONS_MATRIX,nsol,...
    tprec,wat_flxs,simsol,compartiments_number,rtex,dx,wc_change_tol,sol_upper_bc_type,dt_start);

    % Evaluate bc
  	bc(end+1,:)=[t,boco_top_type,boco_top,boco_bot_type,boco_bot,pond,runoff,flxa1];
   
    % Break when no convergence is reached
    if no_conv
        disp('No convergence');
        break
    end
end
close(u);
toc;

% Model output (draw result)
draw_output(wc_profile,tnode,ph_node,wc_node,bot_flux,top_flux,...
    compartiments_number,dx,dx_inter,print_time,print_node,top_inf,bc,...
    simplant,snode,evap,trans,sink,column_depth,units,cumul_evap);
%draw CM and T
draw_cmt(cm_profile,temp_profile,cm_node,tnode,temp_node,compartiments_number,dx,dx_inter,...
    print_time,print_node,nsol,simsol,column_depth,units,bot_sol_flux);







