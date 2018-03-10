function [om_profile,om_node]=store_om(print_time,...
    print_node,soil_om,no_conv,t,om_profile, om_node,sim_nitro);

%Store organic carbon and nitrogen content for output
if no_conv~=1;		
    for j=1:7 %there are 6 organic carbon and nitrogen pools
        for i=1:length(print_time);
            if (t>print_time(i))&(om_profile{j}(i,:)==1);
                if sim_nitro
                    om_profile{j}(i,:) = soil_om(:,j);
                end
            end
        end
    if sim_nitro
        om_node{j,1}(:,end+1)=(soil_om(print_node,j));
    end
    end
end
