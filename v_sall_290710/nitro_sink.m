function[sol_sink,sol_sinkm,sol_sinki,nitrifm,nitrifi,hydro_uream,hydro_ureai,uptakem,uptakei,denitm,...
    deniti,volatm,volati,decnorg,deccorg]=...
    nitro_sink(t,dt,tmax,temp,nsol,ncs,cm,cim,WC,wcmo,wcio,decsolm,decsoli,minerm,mineri,rcarbmin,...
    uptakem,uptakei,soil_parameters,solute_param1,solute_param2,om_param, miner_param, plant_uptake_param,csol,immobile)
                        

%MODULE NITROGEN
%nitro_sink:calculate the amount of sink term
%IN:
%t, dt and tmax= the time, time increment and maximum time(min)
%temp= the temperature
%nsol and ncs= number of solutes and number of compartiments
%cm and cim= solute mobile and immobile concentration   
%WC=water content(-)
%wcmo and wcio = water content in the mobile and immobile region
%decsolm= amount of decayed solute in the mobile region (from script solute_sink)
%sol_sinki= amount of decayed solute in the immobile region (from script solute_sink)
%minerm, mineri= amount of mineralisation or immobilisation (from script miner_immob)
%rcarbmin= organic carbon
%uptakem,uptakei= N uptakke in the mobile and immobile region (from scrpt nitro_uptake)
%soil_parameters= parameters of the soil for water flow
%solute_param1,solute_param2= parameters of the solute (from script in_nitro_param)
%nitro_param,om_param= parameters of the nitrogen and organic matter (from %script in_nitro_param)
%OUT:
%sol_sink= solute sink term
%sol_sinkm and sol_sinki = solute sink term for the mobile and immobile region

%CALL:none
%CALLED BY: SOLVE_SOLUTE
%-------------------------------------------------------------------------
% Mamadou SALL 17/03/09

%%Solute decaying

%%%REM sol_sink is negative if it induces reduction of the solute concentration
%% terms of sol_sink are negative if they induce reduction of the solute

%% concentration
wcr=soil_parameters(:,1);
wcs=soil_parameters(:,2);
alfa=soil_parameters(:,3);
n=soil_parameters(:,4);
m=1-1./n;
alf=solute_param1(:,4);
knitrif=om_param(:,1);  
kdenit=om_param(:,2);   
khydro=om_param(:,3);  
kvolat=om_param(:,4); 
kd=(solute_param2(:,nsol+1:2*nsol));
bd=(solute_param1(:,1));
f=(solute_param1(:,5));
%sol_sinkm = decsolm;

%%calculate reduction factors for temperature and moisture
%critical water content for reduction factor
se=(WC'-wcr)./(wcs-wcr);
thetah=(wcr+(wcs-wcr)./(1+(alfa*100).^n).^m)';       %%WC for ph=100
theta1=(wcr+(wcs-wcr)./(1+(alfa.*1000).^n).^m)';    %%WC for ph=1000
thetaw= (wcr+(wcs-wcr)./(1+(alfa.*150000).^n).^m)';  %%% WC f
%% reduction factor for temperature
redtemp(:,1)=3.^((temp-16)/10);
for i=1:ncs;
 %%   reduction factor for moisture
    if wcs(i) > WC(i) & WC(i) >= thetah(i)
    redmoist(i)=0.5+ 0.5*((wcs(i)-WC(i))/(wcs(i)-thetah(i)));
    elseif thetah(i) > WC(i) & WC(i) >=theta1(i)
    redmoist(i)=1;
    elseif theta1(i) >WC(i) & WC(i) >= thetaw(i)
    redmoist(i)=(WC(i)-thetaw(i))/(theta1(i)-thetaw(i));
    elseif WC(i)>= wcs(i) 
    redmoist(i)=0.5;
    elseif WC(i)< thetaw(i)
    redmoist(i)=0;
    end
% red_kdenit for denitrification (reduction of the dentrification constant)
    rsat(i)=WC(i)/wcs(i);
    if rsat(i)-0.80 < 0.0001
    redkdenit(i)=0;
    else
    redkdenit(i)=((rsat(i)-0.80)/0.20).^2;
    end
end
redkdenit=redkdenit';
redmoist=redmoist';

%denirification
kdenit_act=redtemp.*redkdenit.*kdenit;
denitm=kdenit_act(1:ncs-1).*cm(1:ncs-1,3).*(wcmo(1:ncs-1))';
deniti=kdenit_act(1:ncs-1).*cim(1:ncs-1,3).*(wcio(1:ncs-1))';

%%nitrification
knitr_act = redtemp.*redmoist.*knitrif;
nitrifm = knitr_act(1:ncs-1).*cm(1:ncs-1,2).*(wcmo(1:ncs-1))';
nitrifi = knitr_act(1:ncs-1).*cim(1:ncs-1,2).*(wcio(1:ncs-1))';

%%volatilisation
kvolat_act = redtemp.*redmoist.*kvolat;
volatm = kvolat_act(1:ncs-1).*cm(1:ncs-1,2).*(wcmo(1:ncs-1))';
volati = kvolat_act(1:ncs-1).*cim(1:ncs-1,2).*(wcio(1:ncs-1))';

%%urea hydrolysis
khydro_act=redtemp.*redmoist.*khydro;
hydro_uream=khydro_act(1:ncs-1).*cm(1:ncs-1,1).*(wcmo(1:ncs-1))';
hydro_ureai=khydro_act(1:ncs-1).*cim(1:ncs-1,1).*(wcio(1:ncs-1))';

%calculate the sink term
%%%urea
decsolm(1:ncs-1,1)= decsolm(1:ncs-1,1)- hydro_uream;  %%in the mobile water zone
decsoli(1:ncs-1,1)= decsoli(1:ncs-1,1)- hydro_ureai;
%%ammonia
decsolm(1:ncs-1,2)= decsolm(1:ncs-1,2)+ minerm(1:ncs-1,2)- uptakem(1:ncs-1,2) -nitrifm -volatm + hydro_uream;  %%in the mobile water zone
decsoli(1:ncs-1,2)= decsoli(1:ncs-1,2)+ mineri(1:ncs-1,2)- uptakei(1:ncs-1,2) -nitrifi - volati+ hydro_ureai; 
%%nitrate
decsolm(1:ncs-1,3)= decsolm(1:ncs-1,3)+ minerm(1:ncs-1,3)- uptakem(1:ncs-1,3)+ nitrifm - denitm;  %%in the mobile water zone
decsoli(1:ncs-1,3)= decsoli(1:ncs-1,3)+ mineri(1:ncs-1,3)- uptakei(1:ncs-1,3)+ nitrifi - deniti; 

%%total sol_sink
                                    %%%sol_sink = sol_sinkm+sol_sinki;

%%organic carbon sink
deccorg = -rcarbmin;
%% organic nitrogen sink
decnorg = -minerm(:,3)-mineri(:,3)-minerm(:,2)-mineri(:,2);

%%organic nitrogen sink

%% SINK TERM CORRECTION (from subroutine nit_correct_sink)

for sp = 1:nsol;
%		the mobile sinkterm
		for i=1:ncs-1;
			sol_sinkm(i,sp)=decsolm(i,sp)-alf(i)*(cm(i,sp)-cim(i,sp));
			supply=cm(i,sp)*(wcmo(i)+kd(i,sp)*bd(i)*f(i))/dt;
			if sol_sinkm(i,sp)<0;
				demand=-sol_sinkm(i,sp);
				if demand > supply;
					deficit=supply-demand;
					sol_sinkm(i,sp)= min(sol_sinkm(i,sp)-deficit,0);
					if sp==1; 
                        hydro_uream(i)= max(hydro_uream(i)+deficit,0);
                    end
					if sp==2 
                        nitrifm(i)= max(nitrifm(i)+deficit,0);
                    end
					if sp==3;
                        uptakem(i,3)= max(uptakem(i,3)+deficit,0);
                    end
                end
            end
%		the immobile sink term
               sol_sinki=zeros(ncs-1,3);
			if immobile;
				sol_sinki(i,sp)=decsoli(i,sp);
				supply=csol(i,sp)*(wcio(i)+kd(i,sp)*bd(i)* (1-f(i)))/dt;
				if sol_sinki(i,sp)<0;
					demand=-sol_sinki(i,sp);
					if demand >supply;
						deficit=supply-demand;
						sol_sinki(i,sp)= min(sol_sinki(i,sp)-deficit,0);
						if sp ==1;
                            hydro_ureai(i)=max(hydro_ureai(i)+deficit,0.0);
                        end
						if sp==2;
                            nitrifi(i)=max(rnitrifi(i)+deficit,0);
                        end
						if sp==3;
                            uptakei(i,3)=max(uptakei(i,3)+deficit,0);
                        end
                    end
                end
            end
        end
end
sol_sinkm(ncs,:)=zeros(1,nsol);
sol_sinki(ncs,:)=zeros(1,nsol);
sol_sink = sol_sinkm+sol_sinki;



