function draw_cm(cm_profile,cm_node,tnode,compartiments_number,dx,print_time,...
    print_node,nsol,simsol,column_depth,units);

%Plot the results of the simulation for solute concentration
%-------------------------------------
%M. Sall 17/03/09

Z=0:-dx:-dx*(compartiments_number-1);
z_node=ceil([print_node*column_depth/compartiments_number]);     
if simsol
    for j=1:nsol;
%Profile information
    figure;
    plot(cm_profile{j,1},Z);
    grid on;
    title(['Profile information concentration, solute:',num2str(j)],'FontSize',14,'FontWeight','bold');
    ylabel(['Depth (',units{1,1},')']);xlabel(['concentration (', [units{1,3},' ',units{1,1}],'^-^3',')']);
    legend([repmat('t=',length(print_time),1),num2str(print_time'), repmat(units{1,2},length(print_time),1)]);
    end
end


%selected nodes information 
if simsol
    for j=1:nsol;
    figure;
    plot(tnode,cm_node{j,1});
    grid on;
    title(['Obs. nodes concentration, Solute:',num2str(j)],'FontSize',14,'FontWeight','bold');
    ylabel(['concentration (', [units{1,3},' ',units{1,1}],'^-^3',')']); xlabel(['Time (',units{1,2},')']);
    legend([repmat('z=',length(print_node),1),num2str(z_node'),repmat(units{1,1},length(print_node),1)]);
    end
end

