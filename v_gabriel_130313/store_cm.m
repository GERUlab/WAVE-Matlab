function [cm_profile,cm_node]=store_CM(print_time,...
    print_node,cm,no_conv,t,cm_profile,cm_node,wc_profile,nsol,simsol)

%Store concentration and temperature for output
if no_conv~=1;		
	%wc_profile information;
    for j=1:nsol
        for i=1:length(print_time);
            if (t>print_time(i))&(cm_profile{j}(i,:)==1);
                if simsol
                    cm_profile{j}(i,:) = cm(:,j);
                end
            end
        end
        if simsol
            cm_node{j,1}(:,end+1)=(cm(print_node,j));
        end
    end
end
