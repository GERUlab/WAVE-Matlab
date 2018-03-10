function [delta_stock,mass_balance_error, cum_top_flxs, cum_bot_flxs, cum_sink_wat,sink_prof,cum_evap,cum_inf,sink_wat] = ...
   wat_intgr(ph,wat_flxs,cum_top_flxs ,cum_bot_flxs,cum_sink_wat, stock_initial,dt,dx,...
   compartiments_number,soil_parameters,rtex,cum_evap,cum_inf,esa,sink_wat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the profile mass balance
% M. Vanclooster 17/2/2000
%
% IN: 
%	ph = the soil pressure head
% 	wat_flxs = the soil moisture flux
%  cum_top_flxs = cumulative water flux at the top 
%  cum_bot_flxs = cumulative water flux at the bottom 
%  cum_sink_wat= cumulative sink of water 
% 	stock_initial = the initial soil moisture storage
%	dx,dt,compartiments_number
% OUT: 
% 	Mass balance error
% CALL:
% 	calc_stock,root_extraction_rates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modified by M Sall


% Actual storage (cm)
stock_actual= calc_stock(ph,dx,soil_parameters,dt,0);

% Cumulative flux at the top of the soil profile
cum_top_flxs = cum_top_flxs+wat_flxs(1)*dt;

% Cumulative flux at the bottom of the soil profile
cum_bot_flxs = cum_bot_flxs+wat_flxs(compartiments_number+1)*dt;

% Cumulative transpiration 
sink_prof=cumsum(rtex*dt)*dx;      %% dx modifiי par sall le 09/03/03
cum_sink_wat = cum_sink_wat+ sink_prof; 

%%%%%%%%%%%%%%%%%%%%%%%%% Ajout cumulative transpiration over time only and
%%%%%%%%%%%%%%%%%%%%%%%%% not over depth

sink_wat=sink_wat+sum(rtex*dt)*dx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% actual infiltration and evaporation
if wat_flxs(1) > -1.0E-10
cum_evap=cum_evap+wat_flxs(1)*dt;
    %if rinf > 1.0e-10   %%%%%%rinf = potential infilltration
    %    cum_evap=cum_evap+rinf*dt
    %    cum_inf=cum_inf+rinf*dt
    %end
else
    cum_inf=cum_inf-wat_flxs(1)*dt+esa*dt;
    cum_evap=cum_evap+esa*dt;
end
%%%%ששששששששששששששששששששששששששששששששששששש
% Change of storage (cm)
delta_stock = stock_actual-stock_initial;

% Mass balance error (cm)
%%mass_balance_error = delta_stock + cum_top_flxs - cum_bot_flxs-cum_sink_wat; %Modified by Mamadou
mass_balance_error = delta_stock + cum_top_flxs - cum_bot_flxs + sink_wat;

