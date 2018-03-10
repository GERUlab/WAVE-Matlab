function draw_output(wc_profile,tnode,ph_node,wc_node,bot_flux,top_flux,compartiments_number,dx,dx_inter,...
   print_time,print_node,top_inf,bc,simplant,snode,evap,trans,sink,column_depth,units,cumul_evap)

%Plot the results of the simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z_node=cumsum(dx_inter); 
z_node=Z_node(1,print_node);

%Preliminary calculation
Z=0:-dx:-dx*(compartiments_number-1);

%Profile information
figure;
plot(wc_profile,Z);
grid on;
title('Profile information : water content','FontSize',14,'FontWeight','bold');
ylabel(['Depth (',units{1,1},')']); xlabel(['water content (', [units{1,1},'^3 ',units{1,1}],'^-^3',')'])
legend([repmat('t=',length(print_time),1),num2str(print_time'),repmat(' min.',length(print_time),1)]);

if simplant==1
    figure;plot(sink,Z);grid on;
    title('Profile information : sink term','FontSize',14,'FontWeight','bold');
    ylabel(['Depth (',units{1,1},')']); xlabel(['sink term (', [units{1,1},' ',units{1,2}],')']);
    legend([repmat('t=',length(print_time),1),num2str(print_time'),repmat(' min.',length(print_time),1)]);
end

%Observation points
figure;
plot(tnode,ph_node);
grid on;
title('Observation nodes : pressure heads','FontSize',14,'FontWeight','bold');
ylabel(['Pressure heads (',units{1,1},')']); xlabel(['Time (',units{1,2},')']);
legend([repmat('z=',length(print_node),1),num2str(z_node'),repmat(' cm',length(print_node),1)]);  % Modified by M.Sall z_node (prevously print_node) 

figure;
plot(tnode,wc_node);
grid on;
title('Observation nodes : water content','FontSize',14,'FontWeight','bold');
ylabel(['water content (', [units{1,1},'^3 ',units{1,1}],'^-^3',')']); xlabel(['Time (',units{1,2},')']);
legend([repmat('z=',length(print_node),1),num2str(z_node'),repmat(' cm',length(print_node),1)]);  % Modified by M.Sall

if simplant==1
    figure
    plot(tnode,snode);
    grid on;
    title('Observation nodes : sink','FontSize',14,'FontWeight','bold');
ylabel(['sink term (', [units{1,1},' ',units{1,2}],')']);xlabel(['Time (',units{1,2},')']);
legend([repmat('z=',length(print_node),1),num2str(z_node'),repmat(' cm',length(print_node),1)]);  % Modified by M.Sall

end
%Boundary fluxes
figure;
if simplant==1
    plot(tnode,bot_flux,tnode,evap,tnode,trans,tnode,cumul_evap);
    ylabel(['Output fluxes (', [units{1,1},' ',units{1,2}],'^-^1',')']);xlabel(['Time (',units{1,2},')']);
    legend('Drainage', 'evaporation flux', 'transpiration flux','cumul_evap');
else
    plot(tnode,bot_flux);
    ylabel(['Bottom flux (', [units{1,1},' ',units{1,2}],'^-^1',')']);xlabel(['Time (',units{1,2},')']);
end
grid on;
title('Output flux','FontSize',14,'FontWeight','bold');



figure;
plot(tnode,top_inf,'o',bc(:,1),-bc(:,8),'-+');
ylim([-0.0050 max([top_inf(:);-bc(:,8)])*1.2]);
grid on;
title('Upper boundary conditions','FontSize',14,'FontWeight','bold');
ylabel(['Instantaneous flux (', [units{1,1},' ',units{1,2}],'^-^1',')']);xlabel(['Time (',units{1,2},')']);
legend('actual upper flux','imposed upper flux');
