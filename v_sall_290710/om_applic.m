function[om_appl]=om_applic(t,dt);

%NITROGEN MODULE: dates of application of C and N 
%IN:
%t=time, dt=time increment
%OUT:
%om_appl: amount of organic carbon (C) and organic nitrogen (N) application at time t
%CALL:none
%CALLED BY: miner_immob
%-------------------------------------------------------------------------
% Mamadou SALL 17/03/09

% time(column1) vs. application of carbon (column2(mg m^-2)) and nitrogen(column3(mg m-2))   

%% om application converted from mg/m² to mg/cm²
                       %%t    C        N
om_application= [  00   225   22.5
                   02   0  0
                   03   0   0];
           
               
    % JJ=find(floor(t-dt)<=om_application(:,1) & om_application(:,1)<=ceil(t));
 %if isempty(JJ);
 %    N_application=0;C_application=0;
 %else
 %    C_application=om_application(JJ,2);
  %   N_application=om_application(JJ,3);   
 %end

 i=max(find(om_application(:,1)<=t));
if numel(i)~=0
     C_application=om_application(i,2)*dt;
     N_application=om_application(i,3)*dt; 
else
     C_application=0;
     N_application=0; 
end
  
 om_appl=[C_application, N_application];
 
 %% Remarks: In Wave.FOR, C_application = carborg and N_application = rnitorg
