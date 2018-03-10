function draw_temp(temp_profile,tnode,temp_node,compartiments_number,dx,...
    print_time,print_node,simtemp,column_depth,units);

%plot the results of the simulation for temperature
if simtemp
Z=0:-dx:-dx*(compartiments_number-1);
z_node=ceil([print_node*column_depth/compartiments_number]);

figure;
plot(temp_profile,Z);
grid on;
title('Profile information : Temperature','FontSize',14,'FontWeight','bold');
ylabel(['Depth (',units{1,1},')']); xlabel(['Temperature (',units{1,4},')']);
legend([repmat('t=',length(print_time),1),num2str(print_time'),repmat(units{1,2},length(print_time),1)]);

figure;
plot(tnode,temp_node);
grid on;
title('Observation nodes : temperature','FontSize',14,'FontWeight','bold');
ylabel(['Temperature (',units{1,4},')']); xlabel(['Time (',units{1,2},')']);
legend([repmat('z=',length(print_node),1),num2str(z_node'),repmat(units{1,1},length(print_node),1)]);
end