function [wc_profile,tnode,wc_node,ph_node,bot_flux,top_flux,top_inf,snode,evap,trans,sink,ph_profile,bot_inf, potential_surface_flux, cum_potential_surface_flux, potential_transp, cum_evap, cum_infiltr,cum_pot_transp,cum_trans,water_storage] = store_data(t,ph,...
   wc_profile,tnode,ph_node,wc_node,print_time,print_node,cum_bot_flxs,cum_top_flxs,...
   bot_flux,top_flux,no_conv,WC,wat_flxs,top_inf,esa,cum_sink_wat,snode,evap,trans,sink_prof,sink,ph_profile,bot_inf, potential_surface_flux, cum_potential_surface_flux, potential_transp, cum_evap, cum_infiltr,flxar,dt,epa,cum_pot_transp, cum_trans,cum_sink_wat2,simplant,water_storage,dx,theta_table)

%Store the output of solve_flow for each time step
%IN:
% 	t,ph,wc_profile,tnode,ph_node,wc_node,print_time,
%	print_node,cum_bot_flxs,bot_flux,no_conv,soil_parameters
%OUT:
%	wc_profile,tnode,wc_node,ph_node,bot_flux
%CALLS:
% moist_ret
% if no convergence, then  output are not stored
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if no_conv~=1		
	%wc_profile information
	for i=1:length(print_time)
        if (t>print_time(i))&(wc_profile(i,:)==1)
            wc_profile(i,:) = WC;
            ph_profile(i,:)= ph;
            sink(i,:)=sink_prof;
        end
   	end
   
%Observation nodes
tnode(end+1) = t;
wc_node(:,end+1)=(WC(print_node)).';
ph_node(:,end+1)=(ph(print_node)).';
if isempty(snode)
    snode(:,end+1)=(cum_sink_wat2(print_node)).';
else
    snode(:,end+1)=(cum_sink_wat2(print_node)).'-snode(:,end);
end
%Cumulative bottom flux
bot_flux(end+1)=cum_bot_flxs;
%Instanate bottom flux
bot_inf(end+1) = wat_flxs(end);

%Cumulative top flux = cumulative actual surface flux
top_flux(end+1)=-cum_top_flxs;
%INstantanate top flux = Actual surface flux
top_inf(end+1)=-wat_flxs(1);
%Potential surface flux
potential_surface_flux(end+1) = -flxar;
%Cumulative potential surface flux
cum_potential_surface_flux(end+1) = cum_potential_surface_flux(end)+flxar*dt;

%Potential transpiration
potential_transp(end+1) = epa;
%Cumulative potential transpiration
cum_pot_transp(end+1) = cum_pot_transp(end)+epa*dt;
%Actual transpiration;
trans(end+1) = sink_prof(end)*dx;
%cumulative actual transpiration
cum_trans(end+1)=sum(cum_sink_wat)*dx;


%Potential evaporation
evap(end+1)=esa;

%Potential infiltration (=rinf)(not for output, but calculating cum
%infiltra and evap)
climate = theta_table;
i=max(find(climate(:,1)<=t));
rinf = (climate(i,3)+climate(i,4));

%Cumulative actual infiltration and actual evaporation
if wat_flxs(1) > 0 && esa ~=0
   cum_evap(end+1) = cum_evap(end) + abs(wat_flxs(1))*dt;
   cum_infiltr(end+1) = cum_infiltr(end);
    if rinf > 1*10-10
       cum_evap(end) = cum_evap(end) + rinf*dt;
       cum_infiltr(end) = cum_infiltr(end) + rinf*dt;
    end
else
   cum_infiltr(end+1) = cum_infiltr(end) +abs(wat_flxs(1))*dt+esa*dt;
   cum_evap(end+1) = cum_evap(end)+esa*dt;
end

else 
   	no_conv=0;
  % Write in check.out times where problems of convergence occured
  	 %	fid=fopen('check.out','w');
%		fprintf(fid,'%s\n',t);
end

%Water storage in the profile
water_storage(end+1) = sum(WC.*dx);