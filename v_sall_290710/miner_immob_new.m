function [minerm,mineri,rcarbmin,soil_om]=...
     miner_immob(t,dt,temp,ncs,WC,csol,ph,soil_om,wco,dx,cm,cim,soil_parameters,om_param, miner_param, plant_uptake_param);
      
           
% MODULE NITROGEN
%miner_immob: calculate the nitrogen mineralisation or immobilisation
%IN:
%t= time(min)
%dt = the time increment (min)
%temp=temperature (°C)
%WC:water content at time t, wcp: water content at time t-dt
%cmp:solute mobile concentration at time t-dt
%soil parameters for water flow
%nitro_param =nitrogen parameters (from script in_nitro_om_param)
%om_param = organic matter parameters (from script in_nitro_om_param)
%c_manp,n_manp,c_litp,n_litp,c_hump,n_hump: carbone and nitrogen in manure,
%litter and humus (from in initial om)
%OUT:
%minerm= amount of mineralisation or immobilisation in the mobile zone
%minerm= amount of mineralisation or immobilisation in the immobile zone
%rcarbmin= organic carbon

%CALL: OM_APPLIC
%CALLED BY: SOLVE_SOLUTE

%-------------------------------------------------------------------------
% Mamadou SALL 17/03/09

%% SOIL PARAMETER
wcr=soil_parameters(:,1); 
wcs=soil_parameters(:,2);
alfa=soil_parameters(:,3);
 n=soil_parameters(:,4);
 m=1-1./n;
 %%%%% NITROGEN AND ORGANIC MATER PARAMETER
 ro=miner_param(1);  
 fe=miner_param(2);  
 fh=miner_param(3);  
 k_lito=om_param(:,5); 
 k_mano=om_param(:,6);  
 k_humo=om_param(:,7);  
%%%%%%%%%%%%%%%%INITIAL OM
 c_manp=soil_om(:,1);
 n_manp=soil_om(:,2);
 c_litp=soil_om(:,3);
 n_litp=soil_om(:,4);
 c_hump=soil_om(:,5);
 n_hump=soil_om(:,6);
 
 
volume=dx;
% Find the addition of organic carbon and organic nitrogen due to organic fertilisation
%--------------------------------------------------------------------------
[om_appl]=om_applic(t,dt);
%--------------------------------------------------------------------------
% Consider that organic matter is filled in the first compartiment 
%organic matter application(om_applic) is added in the manure pool and litter pool (see Nitboco.FOR) 
carbman = 0.3*om_appl(1,1); % organic carbon added in the manure pool of the first compariment
rnitman = carbman*(1/ro);   % organic nitrogen added in the manure pool of the first compariment
carblit = om_appl(1,1)-carbman; % organic carbon added in the litter pool of the first compariment
rnitlit = om_appl(1,2)-rnitman; % organic nitrogen added in the litter pool of the first compariment

% new values of organic carbone and organic nitrogen in the first compartiment
c_manp(1)=c_manp(1)+carbman;
c_litp(1)= c_litp(1)+carblit; 
n_manp(1)= n_manp(1)+rnitman;
n_litp(1)= n_litp(1)+rnitlit;

cnlito=n_litp;
cnmano=n_manp;
cchumo=c_hump;
cnhumo=n_hump;
ccmano=c_manp;
cclito=c_litp;


% CALCULATE reduction factor FOR ACTUAL WATER CONTENT and temperature
% temp is the soil temperature (°C)

se=(wco'-wcr)./(wcs-wcr);
thetah=(wcr+(wcs-wcr)./(1+(alfa*100).^n).^m)';       %%WC for ph=100
theta1=(wcr+(wcs-wcr)./(1+(alfa.*1000).^n).^m)';    %%WC for ph=1000
thetaw= (wcr+(wcs-wcr)./(1+(alfa.*150000).^n).^m)';  %%% WC f
%% reduction factor for temperature
redtemp=3.^((temp-16)/10);
for i=1:ncs;
 %%   reduction factor for moisture
    if wcs(i) > wco(i) & wco(i) >= thetah(i)
    redmoist(i)=0.5+ 0.5*((wcs(i)-wco(i))/(wcs(i)-thetah(i)));
    elseif thetah(i) > wco(i) & wco(i) >=theta1(i)
    redmoist(i)=1;
    elseif theta1(i) >wco(i) & wco(i) >= thetaw(i)
    redmoist(i)=((wco(i)-thetaw(i))/(thetal(i)-thetaw(i)));
    elseif wco(i)>= wcs(i) 
    redmoist(i)=0.5;
    elseif wco(i)< thetaw(i)
    redmoist(i)=0;
    end
end
%%%reduction factor
red_fact=redtemp'.*redmoist';
%%%effective coefficient for the actual time
pxm =n_manp/c_manp-fe/ro;     %%rmano=pxm in Wave.FOR ratio defining mineralisation or immobilisation for liiter
pxl=n_litp/c_litp -fe/ro;    %%rlito=pxl in Wave.FOR ratio defining mineralisation or immobilisation for litter
%eff_rate_lit=red_fact.*k_lito;         %%%k_lit =eff_rate_lit
%eff_rate_man=red_fact.*k_mano;         %%%k_man=eff_rate_man   
%eff_rate_hum=red_fact.*k_humo;          %%%%k_hum=eff_rate_hum

rklit=k_lito;
rkman=k_mano;
rkhum=k_humo;
redmoi=redmoist;


for i=1:ncs
    red_fact = redmoi(i) * redtemp(i);
    pxl=(cnlito(i)/cclito(i)-(fe/ro));
    pxm=(cnmano(i)/ccmano(i)-(fe/ro));
    eff_rate_lit = rklit(i)*red_fact;
    eff_rate_man = rkman(i)*red_fact;
    eff_rate_hum = rkhum(i)*red_fact;
    check_mim = ( eff_rate_hum*cnhumo(i)+ pxm*eff_rate_man*ccmano(i)+...
    pxl*eff_rate_lit*cclito(i))*dt;
    if check_mim >= 0.
        change_nh4 = check_mim;
        change_no3 = 0.;
        eff_rate_lit = rklit(i)*red_fact;
        eff_rate_man = rkman(i)*red_fact;
    else
        if((csolo(i,2)-1.8)*wco(i)*dx >= (-check_mim)) 
            change_nh4 = check_mim;
            change_no3 = 0;
            eff_rate_lit = rklit(i)*red_fact;
            eff_rate_man = rkman(i)*red_fact;
        else
            change_nh4 = -(csolo(i,2)-1.8)*wco(i)*dx;
            chno3 = (csolo(i,2)-1.8)*wco(i)*dx+check_mim;
            if((csolo(i,3)-6.2)*wco(i)*dx >= (-chno3)) ;
                change_no3 = chno3;
                eff_rate_lit = rklit(i)*red_fact;
                eff_rate_man = rkman(i)*red_fact;
            else
                change_no3 = -(csolo(i,3)-6.2d0)*wco(i)*dx;
                shortage_min_n = (csolo(i,3)-6.2d0)*wco(i)*dx+chno3;
                if pxm < 0 & pxl >= 0 
                    vimmb = pxm*eff_rate_man*ccmano(i)*dt;
                    vimma = vimmb-shortage_min_n;
                    eff_rate_lit = rklit(i)*red_fact;
                    eff_rate_man = vimma/(pxm*ccmano(i)*dt)
                elseif pxl < 0 & pxm >= 0
                    vimmb = pxl*eff_rate_lit*cclito(i)*dt;
                    vimma = vimmb-shortage_min_n
                    eff_rate_lit = vimma/(pxl*cclito(i)*dt);
                    eff_rate_man = rkman(i)*red_fact;
                elseif pxl< 0 & pxm < 0 
                    vimmb = pxl*eff_rate_lit*cclito(i)*dt;
                    vimma = vimmb-shortage_min_n;
                    if vimma > 0 
                        vimmb = pxm*eff_rate_man*ccmano(i)*dt;
                        vimma = vimmb+vimma;
                        eff_rate_lit = 0;
                        eff_rate_man = vimma/(pxm*ccmano(i)*dt);
                    else
                        error ('troubles in rminimm')
                        %call stop_simulation ('programme stopped : check err_file')
                    end
                else
                    error('troubles in rminimm')
                    % stop_simulation ('programme stopped : check err_file')
                end
            end
        end
    end
    % the gain/losses of the nitrogen litter pool
    % through mineralisation/immobilisation and nitrogen humification
    cnlito(i)=cnlito(i)+(-pxl*eff_rate_lit-fe*fh/ro*eff_rate_lit)*cclito(i)*dt;
    %c the gain/losses of the nitrogen litter pool
    %c through mineralisation/immobilisation and nitrogen humification
    cnmano(i)=cnmano(i)+(-pxm*eff_rate_man-fe*fh/ro*eff_rate_man)*ccmano(i)*dt;
    % the production of co2
   % cco2o(i)=cco2o(i)+(1-fe)*(eff_rate_man*ccmano(i)+eff_rate_lit*cclito(i))*dt+eff_rate_hum*cchumo(i)*dt;
    
   % the production of carbon humus
    cchumo(i)=cchumo(i)+fe*fh*(eff_rate_man*ccmano(i)+eff_rate_lit* cclito(i))*dt-eff_rate_hum*cchumo(i)*dt;
    % the changes in the n-humus pool through mineralisation and n-litter humification
    cnhumo(i)=cnhumo(i)+(fe*fh/ro*(eff_rate_lit*cclito(i)+eff_rate_man*ccmano(i))-eff_rate_hum*cnhumo(i))*dt;
    % the change in the manure pool
    ccmano(i)=ccmano(i)+(-fe*fh-(1-fe))*eff_rate_man*ccmano(i)*dt;
    % the change in the litter pool
    cclito(i)=cclito(i)+(-fe*fh-(1-fe))*eff_rate_lit*cclito(i)*dt;
    % retain the rate of change of the ammonia and nitrate species
    % to be incorporated in the sequential transformation process (mg/(day*liter))
    minerm(i,2)=change_nh4/(dt*volume);%%%rminm(i,2)=change_nh4/(dt*volume)
    minerm(i,3)=change_no3/(dt*volume);
    mineri(i,2)=0.0;
    mineri(i,3)=0.0;
    rcarbmin(i)=(1-fe)*(eff_rate_man*ccmano(i)+...
        eff_rate_lit*cclito(i))/volume + eff_rate_hum*cchumo(i)/volume;
end

n_litp=cnlito;
n_manp=cnmano;
c_hump=cchumo;
n_hump=cnhumo;
c_manp=ccmano;
c_litp=cclito;

soil_om=[c_manp n_manp c_litp n_litp c_hump n_hump];

