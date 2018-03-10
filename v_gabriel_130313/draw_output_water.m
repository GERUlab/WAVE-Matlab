function draw_output_water(wc_profile,ph_profile,tnode,ph_node,wc_node,print_node, print_time,column_depth,dx,units,simplant)

%Plot the results of the water simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compartiments_number = column_depth;

z_node=ceil([print_node*column_depth/compartiments_number]);  
%Preliminary calculation
Z=0:-dx:-dx*(compartiments_number-1);

%Profile information
figure;
plot(wc_profile,Z);
grid on;
title('Profile information : water content','FontSize',14,'FontWeight','bold');
ylabel(['Depth (',units{1,1},')']); xlabel(['water content (', [units{1,1},'^3 ',units{1,1}],'^-^3',')'])
legend([repmat('t=',length(print_time),1),num2str(print_time'),repmat(units{1,2},length(print_time),1)]);

% if simplant==1
%     figure;plot(sink,Z);grid on;
%     title('Profile information : sink term','FontSize',14,'FontWeight','bold');
%     ylabel(['Depth (',units{1,1},')']); xlabel(['sink term (', [units{1,1},' ',units{1,2}],')'])
%     legend([repmat('t=',length(print_time),1),num2str(print_time'),repmat(units{1,2},length(print_time),1)]);
% end

figure;
plot(ph_profile,Z);
grid on;
title('Profile information : pressure head','FontSize',14,'FontWeight','bold');
ylabel(['Depth (',units{1,1},')']); xlabel(['pressure head (', [units{1,1}, ')']]);
legend([repmat('t=',length(print_time),1),num2str(print_time'),repmat(units{1,2},length(print_time),1)]);


%Observation points
figure;
plot(tnode,ph_node);
grid on;
title('Observation nodes : pressure heads','FontSize',14,'FontWeight','bold');
ylabel(['Pressure heads (',units{1,1},')']); xlabel(['Time (',units{1,2},')']);
legend([repmat('z=',length(print_node),1),num2str(z_node'),repmat(units{1,1},length(print_node),1)]);  % Modified by M.Sall z_node (prevously print_node) 

figure;
plot(tnode,wc_node);
grid on;
title('Observation nodes : water content','FontSize',14,'FontWeight','bold');
ylabel(['water content (', [units{1,1},'^3 ',units{1,1}],'^-^3',')']); xlabel(['Time (',units{1,2},')']);
legend([repmat('z=',length(print_node),1),num2str(z_node'),repmat(units{1,1},length(print_node),1)]);  % Modified by M.Sall

% if simplant==1
%     figure
%     plot(tnode,snode);
%     grid on;
%     title('Observation nodes : sink','FontSize',14,'FontWeight','bold');
% ylabel(['sink term (', [units{1,1},' ',units{1,2}],')']);xlabel(['Time (',units{1,2},')']);
% legend([repmat('z=',length(print_node),1),num2str(z_node'),repmat(units{1,1},length(print_node),1)]);  % Modified by M.Sall
% 
% end