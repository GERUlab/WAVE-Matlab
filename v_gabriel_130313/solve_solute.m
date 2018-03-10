function        [cs,csp,cm,cim,cmp,cimp,csol,soil_om,sol_sinkm,sol_sinki,tot_upt,acsolo,acsolio,acsolmo,sflxsa,...
            ddepsol, wdepsol, fsol, conirsol, ppdepth,reservoir,first_time,uptake_matrix,uptakem,...
            solute_applic,uptakei,minerm,mineri,nitrifm,nitrifi,denitm,deniti,hydro_uream,hydro_ureai,...
            volatm,volati,decnorg,deccorg,carbman,carblit,rnitman,rnitlit,om_appl,und,unc,rdens]=...
        solve_solute(t,dt,tmax,temp,dx,nsol,ncs,WC,rtex,cm,cim,csol,cs,acsolmo,acsolio,wat_flxs,ph,...
                    wat_flxsa, wat_flxsah,wcio, wciob, wcma, wcmah, wcmo, wcmob,wco,wcob,pvela,pvelah,...
                    pveloh,soil_om,immobile,soil_parameters,solute_param1,solute_param2,...
                    om_param,miner_param, plant_uptake_param, tot_upt,sim_nitro,tcsolo_ini,diffush,...
                    plant_date,harvest_date,isucr,simplant,ddepsol, wdepsol, fsol, conirsol,...
                    ppdepth,reservoir,tcsolo,first_time,om_appl,drz,NDemand,PLM2,PropFact,initsol,...
                    solboco_input_type,RLengthLa,fraction_time,rdens,ncrop,In_solute_boundary,theta_table);
  
% Solve the 1-D convection dispersion equation 
%using finite differences (cf. WAVE MANUEL)
%%%IN 
%t= time(min)
%dt = the time increment (min)
%temp=temperature (°C)
%nsol=number of solute spieces
%WC:water content at time t, wco: water content=WC
%wat_flxs:(water fluxs ;ph (pressure head)
%solute_param1:some parameters of the solutes [bd,lm,ratio,alf,f]
%solute_param2:some parameters of the solutes [rates,kd,dif,ar,br]
%nitro_param:parameters for nitrogen simulation (from script in_nitro_om_param) 
%soil_parameters: parameters of the soil
%om_param: parameters of the organic matter (from script in_nitro_om_param)
%rtex : root water extraction calculated in the water flow module
%wat_flxsa : water_flxs between the nodes (from script watsol)
%wat_flxsah:water_flxs between the nodes at half time (from script watsol)
%wcio: water content in the immobile region  (from script watso)
%wciob:water content in the immobile region at time t-dt (from script watsol)
%wcma:water content in the mobile region between nodes (from script watsol)
%wcmah: water content in the mobile region between nodes at half time (from script watsol)
%wcmo: water content in the mobile region at the nodes (from script watsol)
%wcmob:water content in the mobile region at the nodes at time t-dt (from script watsol)
%wco: water content at the nodes (from script watsol)
%pvela: pore velocity between the nodes (from script watsol)
%pvelah: pore velocity between the nodes at half time (from script watsol)
%pveloh: pore velocity at the nodes at half time (from script watsol)

%%%OUT
%cm:solute mobile concentration at time t
%cim:solute immobile concentration%at time t
%cs:solute concentration at the soil surface
%c_manp,n_manp,c_litp,n_litp,c_hump,n_hump: carbone and nitrogen in manure,
%litter and humus at time t
%diffus: diffusivity
%sol_sinkm: solute sink for the mobile water
%tot_upt: cumulative nitrogen uptake

%CALL:solute_boundary, solute_sink, miner_immob, nitro_uptake, nitro_sink
%CALLED BY wave_mat main programme
%-------------------------------------------------
%% M. Sall 17/03/09
%Read inputparameters
bd=(solute_param1(:,1));
lm=solute_param1(:,2);
ratio=solute_param1(:,3);
alf=solute_param1(:,4);
lm(ncs+1)=lm(ncs);
f=(solute_param1(:,5));
rates=(solute_param2(:,1:nsol));
kd=(solute_param2(:,nsol+1:2*nsol));
dif=solute_param2(:,2*nsol+1:3*nsol);
ar=solute_param2(:,3*nsol+1:4*nsol);
br=solute_param2(:,4*nsol+1:5*nsol);


%calculate beta coeficients
i = 1:ncs;
betha1(i,1) = wat_flxsa(i)<0;
betha4(i,1) = 1- betha1(i,1);
betha2(i,1)= wat_flxsa(i+1) <= 0 ;
betha3(i,1) = 1- betha2(i,1);

%Determine initial values in the first time step
[cm,cim,csol,acsolmo,acsolio,initsol]=solute_data (cm,cim,csol,acsolmo,acsolio,t,dx,tcsolo_ini,wco,wcmo,wcma,wcio,immobile,solute_param1,solute_param2,nsol,ncs,pvela,initsol);

%Determine the sink term 
[decsolm,decsoli,sol_sinkm,sol_sinki]=solute_sink(dt,nsol,ncs,cm,cim,wcmo,...
    wcio,solute_param1,solute_param2,immobile);
%%
%%Nitrogen module

if sim_nitro;
%     [minerm,mineri,rcarbmin,soil_om]=...
%      miner_immob(t,dt,temp,ncs,WC,csol,ph,soil_om,wco,dx,cm,cim,soil_parameters,nitro_param,om_param);
 [minerm,mineri,rcarbmin,soil_om,om_appl,carbman, rnitman, carblit, rnitlit]=...
     miner_immob_new(t,dt,temp,ncs,WC,csol,ph,soil_om,wco,dx,cm,cim,soil_parameters,om_param, miner_param, plant_uptake_param,om_appl,initsol);

if simplant
[uptakem,uptakei,tot_upt, first_time, uptake_matrix,rdens,und,unc]=nitro_uptake(csol,dt,dx,harvest_date,...
        isucr, ncs,nsol,soil_parameters,plant_uptake_param,solute_param1,solute_param2,...
        plant_date,simplant,t,rtex,tot_upt,wco,first_time,drz,NDemand,PLM2,PropFact,RLengthLa,rdens,ncrop);
%     plot(rdens)
% hold on
% drawnow
% [uptakem,uptakei,tot_upt, first_time, uptake_matrix,rdens]=nitro_uptake2(csol,dt,dx,harvest_date,...
%         isucr, ncs,nsol,soil_parameters,plant_uptake_param,solute_param1,solute_param2,...
%         plant_date,simplant,t,rtex,tot_upt,wco,first_time,drz,NDemand,PLM2,PropFact,RLengthLa);

% [uptakem,uptakei,tot_upt]=nitro_uptake_wavefor(csol,dt,dx,harvest_date,...
%         isucr, ncs,nsol,soil_parameters,om_param, miner_param, plant_uptake_param, solute_param1,solute_param2,...
%         plant_date,simplant,t,rtex,tot_upt,wco);
else
uptakem=zeros(ncs,3);
uptakei=zeros(ncs,3);
tot_upt=[];
uptake_matrix = zeros(1,6);
rdens=[];
unc = zeros(ncs,3); und = zeros(ncs,3);

end

    [sol_sink,sol_sinkm,sol_sinki,nitrifm,nitrifi,hydro_uream,hydro_ureai,uptakem,uptakei,denitm,...
    deniti,volatm,volati,decnorg,deccorg]=...
    nitro_sink(t,dt,tmax,temp,nsol,ncs,cm,cim,WC,wcmo,wcio,decsolm,decsoli,minerm,mineri,rcarbmin,...
    uptakem,uptakei,soil_parameters,solute_param1,solute_param2,om_param, miner_param, plant_uptake_param,csol,immobile);

else
uptake_matrix = zeros(1,6);
minerm= [];mineri=[];nitrifm=[];nitrifi=[];denitm=[];deniti=[];hydro_uream=[];hydro_ureai=[];...
            volatm=[];volati=[];decnorg=[];deccorg=[];carbman=[];carblit=[];rnitman=[];rnitlit=[];cco2o=[];drza=[];
uptakem=zeros(ncs,3);
uptakei=zeros(ncs,3);
tot_upt=0;
uptake_matrix = zeros(1,6);
rdens=[];
unc = zeros(ncs,3); und = zeros(ncs,3);
end

%% actual value to passed value
solute_applic = In_solute_boundary.solute_applic;
solboco_top_type = In_solute_boundary.solboco_top_type;
solboco_input_type = In_solute_boundary.solboco_input_type;
i=max(find(solute_applic(:,1)<=t));
if numel(i)~=0
conc=(solute_applic(i,2:end));
else
    conc=zeros(1,nsol);
end

%[conc,solute_applic,solboco_top_type,solboco_input_type]=In_solute_boundary(t,nsol);
if solboco_input_type ==2
csp=cs./abs(wat_flxsah(1))/dt;
if sum(isnan(csp))== nsol
    csp = zeros(1,nsol);
end
else
csp =cs;
end
cmp=cm;cimp=cim; % cm and cimp = precedent mobile and immobile soute concentration
acsolmb=acsolmo;
%set cmp zero gradient
cmp(ncs+1,:)=cmp(ncs,:);

%% Boundary Condition
%[conc,solute_applic,solboco_top_type,solboco_input_type]=In_solute_boundary(t,nsol);
if solboco_input_type == 1 % As by M. Sall expressed in M L-3
 if wat_flxs(1)<=0
     cs=conc;  
 else
     cs=zeros(1,nsol);
 end
elseif solboco_input_type ==2 %As by J. Vansteenkiste expressed in M L-2
   [fsol, conirsol, solsur] = solute_upper(t, fsol, conirsol, wdepsol, ddepsol,dt,fraction_time,theta_table);

if wat_flxs(1) <0
    cs = (reservoir + solsur)./abs(wat_flxsah(1))/dt;
    cs_temp = (reservoir + solsur);
    reservoir = zeros(1,nsol);
    if sum(cs) ~=0
        %disp('')
    end
else
    reservoir = (reservoir + solsur);
        if sum(reservoir) ~=0
        %disp('')
    end
    cs =zeros(1,nsol);
    cs_temp =cs;
end 
end

%Determine solute redistribution due to plowing
if solboco_input_type == 2
[csol, cm, cim, acsolmo, acsolio,ppdepth,tcsolo] = plowing(t,solute_param1, solute_param2,csol, cm, cim, acsolmo,...
     acsolio, wcio, wco, nsol, immobile,wcmo,tcsolo,dx,ppdepth,solsur);
end


cm = inv_sol_matrix(immobile,ncs,cs,f,dt, pveloh, wat_flxsah, wcmo, wcmob, dx, diffush, wcmah, pvelah,betha1,betha2, betha3,betha4, csp, cmp, bd, kd, sol_sinkm,nsol);

for j=1:nsol
     %%%%    
    for I=1:ncs;
        if cm(I,j)<0
           cm(I,j)=0;
        end
        %Adsorbed mass
    acsolmo(I,j)=kd(I,j)*f(I)*cm(I,j);
 %SOLUTE IMMOBILE CONCENTRATION
        if immobile
            ag(I,j)=wciob(I)/dt+(1-f(I)).*bd(I).*kd(I,j)/dt-alf(I)/2;
            bg(I,j)=wcio(I)/dt+(1-f(I))*bd(I)*kd(I,j)/dt+alf(I)/2;        
            cim(I,j)=(ag(I,j)*cimp(I,j)+(alf((I))*cm(I,j))/2 +...
     				(alf((I))*cmp(I,j))/2)/bg(I,j)+sol_sinki(I,j)/bg(I,j);
			acsolio(I,j)=cim(I,j)*kd((I),j)*(1-f((I)));
			csol(I,j)=(cim(I,j)*wciob(I)+cm(I,j)*wcmob(I))/wcob(I);
			acsolo(I,j)=acsolio(I,j)+acsolmo(I,j);
        else
            cim(I,j)=0.0;
			acsolio(I,j)=0.0;
			csol(I,j)=cm(I,j);
			acsolo(I,j)=acsolmo(I,j);
        end
    end
    
    for i=2:ncs;   
		if wat_flxsah(i)<0  
				sflxsa(i,j)=(wat_flxsah(i)*(cm((i-1),j)+... 
     				cmp((i-1),j))/2)*dt;                
		else
				sflxsa(i,j)=(wat_flxsah(i)*(cm(i,j)+...  
     				cmp(i,j))/2)*dt;   
        end
    end
	if wat_flxsah(1)<0;
  			sflxsa(1,j) = (cs(j)+csp(j))/2*abs(wat_flxsah(1))*dt+(diffush(1,j))*...
                ((cs(j)+csp(j))-(cm(1,j)+cmp(1,j)))/2/dx*(wcmah(1))*dt;
	else
			sflxsa(1,j)=0.0;
    end
		volume= dx;
	if  abs(cm(ncs,j)-cmp(ncs,j)>0.001  & wat_flxsah(ncs+1)>0.00000000001)
			sflxsa((ncs+1),j)=((-cmp(ncs,j)+cm(ncs,j))*...
     		wcmo(ncs)*dx)-(acsolmb(ncs,j)-acsolmo(ncs,j))*...
     		bd(1)*volume;
    elseif wat_flxsah(ncs+1)<0.0
				sflxsa((ncs+1),j)=(wat_flxsah(ncs+1)*...
     			((cm(ncs,j)+ cmp(ncs,j))/2))*dt;
	else
				sflxsa(ncs+1,j)=(wat_flxsah(ncs+1)*cmp((ncs+1),j))*dt;
    end
end

if solboco_input_type ==2
cs= cs_temp;
end

