function [cm_profile,temp_profile,cm_node,temp_node]=store_CMT(print_time,...
    print_node,cm,temp,no_conv,t,cm_profile,temp_profile,cm_node,temp_node,wc_profile,nsol,simsol);

%MODULE NITROGEN



if no_conv~=1;		
	%wc_profile information;
    for j=1:nsol
        for i=1:length(print_time);
            if (t>print_time(i))&(cm_profile{j}(i,:)==1);
                if simsol
                    cm_profile{j}(i,:) = cm(:,j);
                end
            temp_profile(i,:)=temp;
            end
        end
        if simsol
cm_node{j,1}(:,end+1)=(cm(print_node,j));
        end
    end
    temp_node(:,end+1)=(temp(print_node));
end
