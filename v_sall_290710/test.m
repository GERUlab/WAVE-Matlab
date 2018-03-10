%script auxiliaire
close all; 
save wccalc wc_node;
load wccalc;
load nl30;
tcalc=[1:1:256];
plot(nl30(:,1),nl30(:,2),tcalc',wc_node')