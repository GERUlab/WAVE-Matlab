function [delta_stock,mass_balance_error, cum_top_flxs, cum_bot_flxs, cum_sink_wat, cum_sink_wat2,sink_prof] = ...
   wat_intgr(ph,wat_flxs,cum_top_flxs ,cum_bot_flxs,cum_sink_wat, cum_sink_wat2, stock_initial,dt,dx,...
   compartiments_number,soil_parameters,rtex);
% Calculates the profile mass balance
%%%%%%%%%%%%%%%%%%%%
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

% Actual storage (cm)
stock_actual= calc_stock(ph,dx,soil_parameters,dt,0);

% Cumulative flux at the top of the soil profile
% = Cumulative actual surface flux
cum_top_flxs = cum_top_flxs+wat_flxs(1)*dt;

% Cumulative flux at the bottom of the soil profile
cum_bot_flxs = cum_bot_flxs+wat_flxs(compartiments_number+1)*dt;

% Cumulative actual transpiration via actual transpiration (sink_prof)
sink_prof=sum(rtex);
cum_sink_wat = cum_sink_wat+ sum(rtex)*dt;
cum_sink_wat2 = cum_sink_wat2 +dt*rtex;

% Change of storage (cm)
delta_stock = stock_actual-stock_initial;

% Mass balance error (cm)
mass_balance_error = delta_stock + cum_top_flxs - cum_bot_flxs - cum_sink_wat;