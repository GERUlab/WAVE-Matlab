%Check the functionality of the calculation of the reduction function (eq
%2-64 in manual)
clear
clc
close all

%Read input as is done in the main program
climate = in_climatic_data();
[t,tmax,compartiments_number,dx_inter,dx, dt_start,dt_min,print_time,print_node,...
        maxiter,dt_max,column_depth,err_tol,simplant,simsol,isucr,idvs,sim_nitro,...
      immobile,units,plant_date,harvest_date,simtemp,crop_type] = In_general_data();
[drz, rlai, kc, Smax_given,Smax_param,root_length_max] = in_crop_number1(t,plant_date,harvest_date,1);

%Define a range of transpiration rates and pressure heads
%first three will be used for linear calculations, rest hyperbolic. Values chosen so that below
%lower limit, in between lower and upper limit and above upper limit.
different_epa = [0.05 0.25 0.75 0.05 0.25 0.75]; 

ph=0:-01:Smax_param(5)*2;

%Do the calucations
for i = 1: 6
     if i<=3 
     ilinr =1;
 else
     ilinr = 0;
 end
epa = different_epa(i);

p0 = Smax_param(1);p1 = Smax_param(2);p2l = Smax_param(3);   
p2h = Smax_param(4); p3 = Smax_param(5);epa_high = Smax_param(6);
epa_low = Smax_param(7);  ilinr = Smax_param(8); 

p2 = p2h;
 if epa < epa_low
     p2 = p2l;
 elseif epa < epa_high
     p2 = p2h +((epa_high - epa)/(epa_high - epa_low))*(p2l - p2h); %Based on fortran
 end


%Determination of alpha_h for each compartiment
dummy1 = ph >=p0;
dummy2 = ph>=p1 & ph < p0;
dummy3 =  ph >= p2 & ph < p1;
dummy4 = ph >= p3 & ph < p2 ;
dummy5 = ph < p3;

if ilinr
alpha_h =...
    dummy1.*0 +...
    dummy2.*(p0-ph)./(p0-p1) +...
    dummy3*1+...
    dummy4.*ilinr.*((p3-ph)./(p3-p2))+...
    dummy5.*0;
else
    dummy1.*0 +...
    dummy2.*(p0-ph)./(p0-p1) +...
    dummy3*1+...
    dummy4.*((1-ilinr)*(10.^((p2-ph)./p3)))+...
    dummy5.*((1-ilinr)*(10.^((p2-ph)./p3)))
    %Theory: dummy5.*0; but no smooth curve as in figure 2-9
end

alpha_h = alpha_h';
plot_var{i} = alpha_h';
end

subplot(2,2,1)
plot(abs(ph), plot_var{1}(:,:),abs(ph), plot_var{2}(:,:),abs(ph), plot_var{3}(:,:) )
title('linear')
xlabel('pressure head (cm)')
ylabel('reduction function alpha (-)')
legend('below lower', 'in between', 'above upper')
subplot(2,2,3)
plot(abs(ph), plot_var{4}(:,:),abs(ph), plot_var{5}(:,:),abs(ph), plot_var{6}(:,:) )
title('hyperbolic')
xlabel('pressure head (cm)')
ylabel('reduction function alpha (-)')
legend('below lower', 'in between', 'above upper')
subplot(2,2,2)
semilogx(abs(ph), plot_var{1}(:,:),abs(ph), plot_var{2}(:,:),abs(ph), plot_var{3}(:,:) )
title('linear')
xlabel('pF(log 10 cm)')
ylabel('reduction function alpha (-)')
legend('below lower', 'in between', 'above upper')
subplot(2,2,4)
semilogx(abs(ph), plot_var{4}(:,:),abs(ph), plot_var{5}(:,:),abs(ph), plot_var{6}(:,:) )
title('hyperbolic')
xlabel('pF(log 10 cm)')
ylabel('reduction function alpha (-)')
legend('below lower', 'in between', 'above upper')

% %Compare calclation methods of p2
% epa = epa_low:0.01:epa_high;
%      p2_H = p2h+(p2l-p2h)/(epa_high - epa_low)./(epa_high-epa); %based on Hydrus
%      p2_F = p2h +((epa_high - epa)./(epa_high - epa_low))*(p2l - p2h); %Based on fortran
%      figure
% plot(p2_H, p2_F); xlabel('hydrus'); ylabel('Fortran');
% 
% figure
% plot(epa, p2_H); xlabel('epa'); ylabel('hydrus')
% plot(epa, p2_F); xlabel('epa'); title('fortran code'); ylabel('p2')
