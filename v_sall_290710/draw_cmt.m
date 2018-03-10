function draw_cmt(cm_profile,temp_profile,cm_node,tnode,temp_node,compartiments_number,dx,dx_inter,...
    print_time,print_node,nsol,simsol,column_depth,units,bot_sol_flux);

%module solute
%plot the results of the simulation for solute concentration and
%temperature

%%%Call: none
%Called by: wavemat105_sol


%-------------------------------------
%M. Sall 17/03/09

Z=0:-dx:-dx*(compartiments_number-1);
Z_node=cumsum(dx_inter); 
z_node=Z_node(1,print_node);
if simsol
    for j=1:nsol;
%Profile information
    figure;
    plot(cm_profile{j,1},Z);
    grid on;
    title(['Profile information concentration, solute:',num2str(j)],'FontSize',14,'FontWeight','bold');
    ylabel(['Depth (',units{1,1},')']);xlabel(['concentration (', [units{1,3},' ',units{1,1}],'^-^3',')']);
    legend([repmat('t=',length(print_time),1),num2str(print_time'),repmat(' min.',length(print_time),1)]);
    end
end
figure;
plot(temp_profile,Z);
grid on;
title('Profile information : Temperature','FontSize',14,'FontWeight','bold');
ylabel(['Depth (',units{1,1},')']); xlabel(['Temperature (',units{1,4},')']);
legend([repmat('t=',length(print_time),1),num2str(print_time'),repmat(' min.',length(print_time),1)]);


if simsol
    %selected nodes information 
    for j=1:nsol;
    figure;
    plot(tnode,cm_node{j,1});
    grid on;
    title(['Obs. nodes concentration, Solute:',num2str(j)],'FontSize',14,'FontWeight','bold');
    ylabel(['concentration (', [units{1,3},' ',units{1,1}],'^-^3',')']); xlabel(['Time (',units{1,2},')']);
    legend([repmat('z=',length(print_node),1),num2str(z_node'),repmat(' cm',length(print_node),1)]);
  
    %ששששששששששששplot solute bottom flux for al solutes
    figure;
    plot(tnode,bot_sol_flux{j,1});
    title(['cumulative solute bottom flux:',num2str(j)])
    legend(['solute1:',num2str(j)]);
    
    end
end
figure;
plot(tnode,temp_node);
grid on;
title('Observation nodes : temperature','FontSize',14,'FontWeight','bold');
ylabel(['Temperature (',units{1,4},')']); xlabel(['Time (',units{1,2},')']);
legend([repmat('z=',length(print_node),1),num2str(z_node'),repmat(' cm',length(print_node),1)]);




