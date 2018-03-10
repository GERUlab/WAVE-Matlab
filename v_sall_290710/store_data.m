function [wc_profile,tnode,wc_node,ph_node,bot_flux,top_flux,top_inf,snode,evap,trans,sink,cumul_evap] = store_data(t,ph,...
   wc_profile,tnode,ph_node,wc_node,print_time,print_node,cum_bot_flxs,cum_top_flxs,...
   bot_flux,top_flux,no_conv,WC,wat_flxs,top_inf,esa,cum_sink_wat,snode,evap,trans,sink_prof,sink,ncs,cum_evap,...
   cum_inf,cumul_evap,sink_wat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            sink(i,:)=sink_prof;
        end
   	end

   
%Observation nodes
tnode(end+1) = t;
wc_node(:,end+1)=(WC(print_node)).';
ph_node(:,end+1)=(ph(print_node)).';
if isempty(snode)
    snode(:,end+1)=(cum_sink_wat(print_node)).';
else
    snode(:,end+1)=(cum_sink_wat(print_node)).'-snode(:,end);
end

%Bottom flux
bot_flux(end+1)=cum_bot_flxs;
%transpiration
%trans(end+1)=sum(cum_sink_wat);Modified Mamadou
trans(end+1)=(cum_sink_wat(ncs));

%evaporation
evap(end+1)=esa;

%cumulative evaporation
cumul_evap(end+1)=cum_evap;

%Cumulative top flux
top_flux(end+1)=-cum_top_flxs;
%INstantanate top flux
top_inf(end+1)=-wat_flxs(1);
else 
   	no_conv=0;
  % Write in check.out times where problems of convergence occured
  	 %	fid=fopen('check.out','w');
%		fprintf(fid,'%s\n',t);
end
