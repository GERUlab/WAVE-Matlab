function [tc_node1_new, tc_node2_new] = change_tc(tc_node1_new, tc_node2_new,...
    fract1_new, fract2_new,tcsolo1,tcsolo2)

for i=1:3
if fract1_new==1
tc_node1_new{i,1}(:,end+1) = tcsolo1(:,i);
tc_node2_new{i,1}(:,end+1) = tcsolo1(:,i)*0;

elseif fract2_new==1
tc_node2_new{i,1}(:,end+1) = tcsolo2(:,i);
tc_node1_new{i,1}(:,end+1) = tcsolo2(:,i)*0;

else
tc_node2_new{i,1}(:,end+1) = tcsolo2(:,i);
tc_node1_new{i,1}(:,end+1) = tcsolo1(:,i);
end
end


