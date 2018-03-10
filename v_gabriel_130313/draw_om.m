function draw_om(om_profile,om_node,tnode,compartiments_number,dx,...
    print_time,print_node,nsol,sim_nitro,column_depth,units);

%Plot the results of the simulation for organic matter decay
%(manure,litter, humus)

%Nodal informations are summarized for all pools and then plotted
%-------------------------------------
%c_manp n_manp c_litp n_litp c_hump n_hump
Z=0:-dx:-dx*(compartiments_number-1);
z_node=ceil([print_node*column_depth/compartiments_number]);     

%Profile information
if sim_nitro
    figure
    for j=1:6;
    %Profile information
    subplot(3,2,j)
    plot(om_profile{j,1},Z);
    grid on;
    if j==1
        title ('Carbon manure')
    elseif j==2
        title ('Nitrogen manure')
    elseif j==3
        title ('Carbon litter')        
    elseif j==4
        title ('Nitrogen litter')
    elseif j==5
        title ('Carbon humus')        
    elseif j==6
        title ('Nitrogen humus')
    end         
    ylabel(['Depth (',units{1,1},')']);
    xlabel(['concentration (', units{1,3}, '/', units{1,1}, '^2)']);
    legend([repmat('t=',length(print_time),1),num2str(print_time'),repmat(' min.',length(print_time),1)]);
    hold on
    end
    hold off
end

%selected nodes information 
if length(print_node)>20
    error('output can not be put in graph for more than 20 nodes, please reduce the number of print_nodes ,,in in_general_data')
end


if sim_nitro
    for i=1:length(print_node)
    figure
    var =[];
    for j=1:6;
        var(:,j) = om_node{j,1}(i,:);
    end
    plot(tnode,var);
    grid on;
    title(['Obs. nodes concentration, OM:node = ',num2str(print_node(i))],'FontSize',14,'FontWeight','bold');
    ylabel(['concentration (', units{1,3}, '/', units{1,1}, '^2)']); 
    xlabel(['Time (',units{1,2},')']);
    legend('C_{man}', 'N_{man}', 'C_{lit}', 'N_{lit}', 'C_{hum}' ,'N_{hum}')
    end
end

%selected nodes information 
if sim_nitro
    for i=1:length(print_node)
    figure
    var =[];
    for j=1:1:7;
        var(:,j) = om_node{j,1}(i,:);
    end
            sumpools = var(:,1)+var(:,3)+var(:,5)+var(:,7);
    plot(tnode,var(:,1:2:end), tnode, sumpools);
    grid on;
    title(['Obs. nodes concentration, OM:node = ',num2str(print_node(i))],'FontSize',14,'FontWeight','bold');
    ylabel(['concentration (', units{1,3}, '/', units{1,1}, '^2)']); 
    xlabel(['Time (',units{1,2},')']);
    legend('C_{man}', 'C_{lit}', 'C_{hum}' ,'CO2', 'sum of all carbon pools')
    end
end


