function [cm_profile,temp_profile,cm_node,temp_node,bot_sol_flux]=in_cmt_profile(nsol,print_time,compartiments_number);
    
%% initialisation of the and selected node and profile information for solute and temperature
%IN: 
%--------------------------------------
% M. Sall 17/03/09

for j=1:nsol
cm_profile{j,1}=ones(length(print_time),compartiments_number);%%MAMADOU
cm_node{j,1}=[];
temp_node=[];
bot_sol_flux{j,1}=[];
end
temp_profile=ones(length(print_time),compartiments_number);%%MAMADOU


