function [tnode, top_inf,trans, bot_flux, top_flux,evap,cum_evap,...
    cum_infiltr,bot_inf,potential_surface_flux,...
    potential_transp,cum_potential_surface_flux,...
    cum_pot_transp,cum_trans,water_storage,root_length_time,om_appl,...
    wat_flxs,pvela, pvelah,wcma, wcmah,wat_flxsa,...
    pvelo, pveloh,wcio, wciob, wcmob, wcmo, wco, wcob,...
    decnorg,deccorg,...
    nitrifm,nitrifi,denitm,deniti,hydro_uream,hydro_ureai,...
    volatm,volati,bc,time_uptake,cum_nit_sink,om_balance,...
    N_reaction_balance,...
    snode,root_density_time2,root_density_time3,root_dens_time,...
    sink,bctop_changed,boco_top_type,boco_top,bcbot_changed,...
    boco_bot_type,boco_bot,cum_top_flxs,cum_bot_flxs,cum_sink_wat,...
    ponded,pond,pond_from,phsurf,phbot,flxar,flxsbot,runoff,...
    dt_changed_bc,next_dt_new_bc, iter,...
    first_time_bc,applic_boolean,tot_upt,carbman,carblit,rnitman,...
    rnitlit,tflnorg,tflcorg,cco2o,diffus,csol,acsolmo,acsolio,...
    reservoir,cm,cim,uptakem,uptakei,minerm,mineri,...
    seep,case_breaking,flxsa,flxsah,dt_sol_count,...
    first_time,drza,uptake_matrix,front,initsol] = initialize();



%% Initialization
%Empty and 1 row, time(+-1) columns
tnode=[];top_inf=[];trans=[];bot_flux=[];top_flux=[];evap=[];
cum_evap=[];cum_infiltr=[];bot_inf=[];potential_surface_flux=[];
potential_transp=[];cum_potential_surface_flux=0;cum_evap= 0;
cum_infiltr=0;cum_pot_transp=0;cum_trans=0;water_storage=[];

root_length_time = [];
om_appl=[];

%empty and 1 row, ncomp+1 columns
wat_flxs=[];pvela=[]; pvelah=[];wcma=[]; wcmah=[];wat_flxsa=[];

%empty and 1 row, ncomp columns
pvelo=[]; pveloh=[];wcio=[]; wciob=[]; wcmob=[]; wcmo=[]; wco=[]; wcob=[];
decnorg=0;
deccorg=0;

%empty and ncomp-1 rows, 1 column
nitrifm=[];nitrifi=[];
denitm=[];deniti=[];
hydro_uream=[];hydro_ureai=[];
volatm=[];volati=[];

%Empty and time rows, x columns
bc=[];time_uptake=[];cum_nit_sink=[];om_balance=[];
N_reaction_balance=zeros(1,10);

%Empty and ncomp rows, time columns,
snode=[];root_density_time2=[];root_density_time3=[];root_dens_time = [];

%1 number empty
sink=[];

% 1 number zero
bctop_changed=[0];boco_top_type=[0];boco_top=[0];
bcbot_changed=[0];boco_bot_type=[0];boco_bot=[0];cum_top_flxs = 0;
cum_bot_flxs = [0]; cum_sink_wat = 0;ponded=0;pond=0;pond_from=9999;
phsurf=0;phbot=0;flxar=0;flxsbot=0;runoff=0;
dt_changed_bc=0;next_dt_new_bc=0; iter=0;
first_time_bc=1;applic_boolean =0;tot_upt=0;
carbman=[0];carblit=[0];rnitman=0;rnitlit=0;
tflnorg=0;tflcorg=0;cco2o=0;

%ncomp rows,nsol columns
diffus=[];csol=[];
acsolmo=[];
acsolio=[];
reservoir=0;
cim=[];
cm=[];
uptakem=[];uptakei=[];
minerm=[];mineri=[];

% ??
seep=0;case_breaking=0;runoff=0;
flxsa=[]; flxsah=[];tot_upt=0;dt_sol_count=0;
first_time=1;drza=0;
uptake_matrix=zeros(1,6);applic_boolean =0; front = []; initsol=0;