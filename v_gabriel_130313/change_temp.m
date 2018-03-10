function [temp_node1_new, temp_node2_new] = change_temp(temp_node1_new, temp_node2_new,...
    fract1_new, fract2_new,temp1,temp2)

if fract1_new==1
temp_node1_new(:,end+1) = temp1';
temp_node2_new(:,end+1) = temp1'*0;

elseif fract2_new==1
temp_node2_new(:,end+1) = temp2';
temp_node1_new(:,end+1) = temp2'*0;
else
temp_node2_new(:,end+1) = temp2';
temp_node1_new(:,end+1) = temp1';
end
