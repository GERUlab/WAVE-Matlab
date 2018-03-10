function [temp_node, temp_profile]=store_T(print_time,...
    print_node,no_conv,t, temp, temp_profile, temp_node)

%Store concentration and temperature for output
if no_conv~=1;		
	%wc_profile information;
        for i=1:length(print_time);
            if (t>print_time(i))&(temp_profile(i,:)==1);
            temp_profile(i,:)=temp;
            end
        end
    temp_node(:,end+1)=(temp(print_node));
end