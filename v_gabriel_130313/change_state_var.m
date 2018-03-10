function [ph1, ph2, ph_node1_new, ph_node2_new,wc_node2_new,wc_node1_new] = change_state_var(wcmo1,wcmo2,ph_node1_new, ph_node2_new,...
    fract1_new, fract2_new,soil_parameters, wc_node2_new,wc_node1_new)

if fract1_new==1
ph1 = calc_ph_from_wc(wcmo1, soil_parameters);
ph_node1_new(:,end+1) = ph1';
wc_node1_new(:,end+1) = wcmo1';
wc_node2_new(:,end+1) = wcmo1'.*0;
ph2 =[];
ph_node2_new = ph_node2_new;

elseif fract2_new==1
ph2 = calc_ph_from_wc(wcmo2, soil_parameters);
ph_node2_new(:,end+1) = ph2';
wc_node2_new(:,end+1) = wcmo2';
ph1=[];
ph_node1_new = ph_node1_new;

else
ph1 = calc_ph_from_wc(wcmo1, soil_parameters);  
ph2 = calc_ph_from_wc(wcmo2, soil_parameters);
ph_node2_new(:,end+1) = ph2';
ph_node1_new(:,end+1) = ph1';
wc_node2_new(:,end+1) = wcmo2';
wc_node1_new(:,end+1) = wcmo1';
end


