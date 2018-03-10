function [rdens, drza] = calc_rdens2(DAP,drz,ncrop,crop_type)

% [tunused,tmax,compartiments_number,dx_inter,dx, dt_start,dt_min,print_time,print_node,...
%         maxiter,dt_max,column_depth,err_tol,simplant,simsol,isucr,idvs,sim_nitro,...
%       immobile,units,plant_date,harvest_date,simtemp,crop_type] = In_general_data();


compartiments_number = 85;

if crop_type ==1
%for DAP = 1:100
rootDAP = [0 20 40 68 120];
z1_param = [2.98 2.98 1.77 1.68 1.68];
z2_param = [1.36 1.36 5.24 7.57 7.57];

rootdensparam1 =interp1(rootDAP,z1_param,DAP);
rootdensparam2 =interp1(rootDAP,z2_param,DAP);
%depth = 1:1:compartiments_number;
depth = 1:1:85;


rdens = gampdf(depth,rootdensparam1,rootdensparam2);
% plot(depth,rdens)
% hold on
% drawnow
% end
rdens(drz+1:end) = 0;
drza= drz;
maxrdens = max(rdens);
if DAP < 35
rdens = 25*rdens/maxrdens;
elseif DAP > 55 
    rdens = 30*rdens/maxrdens;
else
    rdens = 40*rdens/maxrdens;
end


elseif crop_type==2
%for DAP = 0:10:140
rootDAP = [0 32 68 102 180];
 z1_param = [7.45 7.45 3.05 2.63 2.63];
 z2_param = [0.5 1.96 4.99 5.34 5.34];


rootdensparam1 =interp1(rootDAP,z1_param,DAP);
rootdensparam2 =interp1(rootDAP,z2_param,DAP);
depth = 1:1:compartiments_number;
%depth =1:1:85;
rdens = gampdf(depth,rootdensparam1,rootdensparam2);
drz = round(drz);
rdens(drz+1:end) = 0;
%  plot(depth,rdens)
%  hold on
%  drawnow
%drza= drz;
drza= drz;
maxrdens = max(rdens);
if DAP < 35
rdens = 25*rdens/maxrdens;
elseif DAP > 55 
    rdens = 30*rdens/maxrdens;
else
    rdens = 40*rdens/maxrdens;
end

else
    disp('oei')
end


end
