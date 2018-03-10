 function [wc_profile, wc_node, cm_profile,temp_profile,cm_node,temp_node,om_profile, om_node,tc_profile,tc_node,solute_balance,ph_node,ph_profile]=...
    create_profile(nsol,print_time,compartiments_number,t)

%Create profile variables which will be updated every timestep and printing
%time

%Water
wc_node =[];
wc_profile=ones(length(print_time),compartiments_number);

ph_node = [];
ph_profile = ones(length(print_time),compartiments_number);

%Solute
for j=1:nsol
cm_profile{j,1}=ones(length(print_time),compartiments_number);
cm_node{j,1}=[];
solute_balance{j,1}=zeros(1,6);
solute_balance{j,1}(1,1)=t;
end

%Temperature
temp_node=[];
temp_profile=ones(length(print_time),compartiments_number);

%Organic matter 
%(7 pools: C_manure N_manure C_litter N_litter C_humus N_humus CO2)
for j=1:7
om_profile{j,1}=ones(length(print_time),compartiments_number);
om_node{j,1}=[];
end

for j=1:nsol
tc_profile{j,1}=ones(length(print_time),compartiments_number);%%MAMADOU
tc_node{j,1}=[];
end