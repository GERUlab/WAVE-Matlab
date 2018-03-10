function [minerm,mineri,rcarbmin,soil_om]=...
     miner_immob(t,dt,temp,ncs,WC,csol,ph,soil_om,wco,dx,cm,cim,soil_parameters,nitro_param,om_param);
      
%Calculate the nitrogen mineralisation or immobilisation

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
 wcs=soil_parameters(:,2);
 wcr=soil_parameters(:,1);
 alfa=soil_parameters(:,3);
 n=soil_parameters(:,4);
 m=1-1./n;
 %%%%% NITROGEN AND ORGANIC MATER PARAMETER
 ro=nitro_param(6);
 fe=nitro_param(7);
 fh=nitro_param(8);
 k_lito=om_param(:,1);
 k_mano=om_param(:,2);
 k_humo=om_param(:,3);
%%%%%%%%%%%%%%%%INITIAL OM
 c_manp=soil_om(:,1);
 n_manp=soil_om(:,2);
 c_litp=soil_om(:,3);
 n_litp=soil_om(:,4);
 c_hump=soil_om(:,5);
 n_hump=soil_om(:,6);
 cco2p = soil_om(:,7);
 
 
volume=dx;
% Find the addition of organic carbon and organic nitrogen due to organic fertilisation
%--------------------------------------------------------------------------

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


% CALCULATE reduction factor FOR ACTUAL WATER CONTENT and temperature
% temp is the soil temperature (°C)
rmano=n_manp/c_manp-fe/ro;     %%pxm in Wave.FOR
rlito=n_litp/c_litp -fe/ro;    %%pxl in Wave.FOR
se=(wco'-wcr)./(wcs-wcr);
wc1=(wcr+(wcs-wcr)./(1+(alfa*100).^n).^m)';       %%WC for ph=100
wc2=(wcr+(wcs-wcr)./(1+(alfa.*1000).^n).^m)';    %%WC for ph=1000
wc3= (wcr+(wcs-wcr)./(1+(alfa.*15000).^n).^m)';  %%% WC f
%% reduction factor for temperature
redtemp=3.^((temp-16)/10);
for i=1:ncs;
 %%   reduction factor for moisture
   if wco(i)>wcs(i)                     %%%0-wcs
       redmoist(i)=0.5;
   elseif wco(i)<wcs(i) & wco(i)>wc1(i);         %%%wcs-100
        redmoist(i)=0.5+0.5*(wcs(i)-wco(i))/(wcs(i)-wc2(i));
   elseif wc3(i)<wco(i) & wco(i)<wc2(i);        %%%% 1000-15000
        redmoist(i)=(WC(i)-wc3(i))/(wc1(i)-wc3(i)); 
   elseif wco(i)<wc3(i)                       %%%%%%%%%%%%>15000
        redmoist(i)=0;
   else                                     %%%%%%%%%%%%%%%100-1000
       redmoist(i)=1;                   
   end
end
%%%reduction factor
red_fact=redtemp'.*redmoist';
%%%effective coefficient for the actual time

k_lit=red_fact.*k_lito; 
k_man=red_fact.*k_mano; 
k_hum=red_fact.*k_humo;

 %CALCULATE AMOUNT OF MINERALISATION

for p=[1:ncs];
    miner(p)= dt*((n_litp(p)./c_litp(p) - fe/ro).*k_lit(p).*c_litp(p)+ ...
    (n_manp(p)./c_manp(p) - fe/ro).*k_man(p).*c_manp(p)+... 
    k_hum(p).*n_hump(p));                 %+NH4P 

    if miner(p)>=0;                                       % MINERALISATION
        ch_nh4(p)=miner(p);          %%change on nh4= increase of nh4
        ch_no3(i)=0;                         %%%%%no change on no3
    else                          %%%%in miner <0 then immobilisation will be done
      if csol(p,2)*wco(p)*dx >= abs(miner(p));      %immobilisation and USE NH4
          ch_nh4(p)=-abs(miner(p));    %%% decrease of nh4
          ch_no3(p)=0;                 %%% no change on no3
      else                                   %%%if there is not enough nh4
          ch_nh4(p)=-csol(p,2)*wco(p)*dx;              %%% the remaining nh4 is immobilised
             if csol(p,3)*wco(p)*dx >abs(miner)-csol(p,2)*wco(p)*dx;            %% no3 is greater than the immobilisation gap
                ch_no3(p)=-(abs(miner)-csol(p,2)*wco(p)*dx);      %%no3 supply for the gap of immobilisation
             else                                    %% there is not enough no3 to fill the gap
               ch_no3(p)=-csol(p,3)*wco(p)*dx;                      %the remaining no3 is used
               %%%Shortage for mineralisation voir fortan programme
               shortage_min(p)=csol(p,3)*wco(p)*dx+ch_no3(p);
               if rmano(p)<0 & rlito(p)>0
                   vimmb(p) = rmano(p)*k_man(p)*c_manp(p)*dt;
                   vimma(p) = vimmb(p)-shortage_min(p); 
                   k_lit(p,1) = k_lito(p)*red_fact(p);
                   k_man(p,1) = vimma(p)/(rmano(p)*c_manp(p)*dt);
               elseif rmano(p)>0 & rlito(p)<0  
                                 
                   
                   vimmb(p) = rlito(p)*k_lit(p)*c_litp(p)*dt;
                   vimma(p) = vimmb(p)-shortage_min(p);  
				   k_lit(p,1)= vimma(p)/(rlito(p)*c_litp(p)*dt);
				   k_man(p,1)= k_mano(p)*red_fact(p);
               elseif rmano(p) < 0 & rlito(p) < 0 
                   vimmb(p) = rlito(p)*k_lit(p)*c_litp(p)*dt;
				   vimma(p) = vimmb(p)-shortage_min(p);  
                   if vimma(p) > 0 
                       vimmb(p)= rmano(p)*k_man_man*c_manp(p)*dt;
					   vimma(p)= vimmb(p)+vimma(p);
					   k_lit(p,1)= 0;
					   k_man(p,1)= vimma/(rmano*c_manp(p)*dt);
                   else
                       error('troubles in miner_immob shortage');
                   end
               end
             end
      end
    end
end

 
%CALCULATE organic mater(OM) diminution or increase due to MINERALISATION
%and IMMOBILISATION   
  % Equation 6-1 wave manual
  % Equation 6-2
%%the gain/losses of the nitrogen litter pool through mineralisation/immobilisation and nitrogen humification
n_litp=n_litp+ dt*((1-fh)*fe*1/ro-n_litp./c_litp).*k_lit.*c_litp;

%the gain/losses of the nitrogen litter pool through mineralisation/immobilisation and nitrogen humification% Equation 6-4
n_manp=n_manp+ dt*((1-fh)*fe*1/ro-n_manp./c_manp).*k_man.*c_manp;

%%the production of co2

cco2o=cco2p + (1-fe)*(k_man.*c_manp(i)+k_lit.*c_litp(i))*dt+k_hum.*c_hump(i)*dt;

%%the production of carbon humus% Equation 6-5
c_hump=c_hump+dt*(fe*fh*(k_lit.*c_litp+k_man.*c_manp)-k_hum.*c_hump);

%%the changes in the n-humus pool through mineralisation and n-litter humification
n_hump=n_hump+dt*(((fe*fh)/ro)*(k_lit.*c_litp+k_man.*c_manp)-k_hum.*n_hump);

%% the change in the manure pool   Equation 6-3
c_manp=c_manp+ dt*((1-fh)*fe-1)*k_man.*c_manp;

%the change in the litter pool
c_litp=c_litp+ dt*((1-fh)*fe-1)*k_lit.*c_litp;

%retain the rate of change of the ammonia and nitrate species to be incorporated in the sequential transformation process  (mg/(day*liter))
%there is not change in urea due to mineralisation 
ch_urea=zeros(ncs,1);
%matrix containing mineralisation amount in the mobile zone for this time step
minerm=[ch_urea,ch_nh4',ch_no3']/(dt*volume);
%%not mineralisation for immobile zone
mineri=zeros(ncs,3);
rcarbmin=(1-fe)*(k_man.*c_manp+k_lit.*c_litp)/volume + k_hum.*c_hump/volume;

soil_om=[c_manp n_manp c_litp n_litp c_hump n_hump cco2o];   



