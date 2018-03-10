function [tc_profile,tc_node]=store_tc(print_time,...
    print_node,tcsolo,no_conv,t,tc_profile,tc_node,nsol,simsol);

%Store concentration for output (in M L-2)
if no_conv~=1;		
	%wc_profile information;
    for j=1:nsol
        for i=1:length(print_time);
            if (t>print_time(i))&(tc_profile{j}(i,:)==1);
                if simsol
                    tc_profile{j}(i,:) = tcsolo(:,j);
                end
            end
        end
        if simsol
            tc_node{j,1}(:,end+1)=(tcsolo(print_node,j));
        end
    end
end
