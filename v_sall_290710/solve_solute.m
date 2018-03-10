function[cs,cm,cim,cmp,cimp,csol,soil_om,sol_sinkm,sol_sinki,tot_upt,acsolo,acsolio,acsolmo,sflxsa,...
    solute_applic,uptakem,uptakei,minerm,mineri,nitrifm,nitrifi,denitm,deniti,hydro_uream,hydro_ureai,...
    volatm,volati,decnorg,deccorg,carbman,carblit,rnitman,rnitlit,cco2o]=...
    solve_solute(t,dt,tmax,temp,dx,nsol,ncs,WC,rtex,cm,cim,csol,cs,acsolmo,wat_flxs,ph,wat_flxsa,...
    wat_flxsah,wcio, wciob, wcma, wcmah, wcmo, wcmob,wco,wcob,pvela,pvelah,pveloh,soil_om,immobile,...
    soil_parameters,solute_param1,solute_param2,om_param, miner_param, plant_uptake_param,tot_upt,...
    sim_nitro,tcsolo_ini,diffush,plant_date,harvest_date,isucr,simplant,tprec,uptakem,uptakei,minerm,...
    mineri,nitrifm,nitrifi,denitm,deniti,hydro_uream,hydro_ureai,volatm,volati,decnorg,deccorg,carbman,...
    carblit,rnitman,rnitlit,cco2o,sol_upper_bc_type);
    
    
%%%%solve_solute:Solve the 1-D convection dispersion equation using finite
%%%%differences (cf. WAVE MANUEL)
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
% M. Sall 17/03/09


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
for I = 1:ncs;
    if wat_flxsa(I) < 0;
        betha1(I,1)=1;
		betha4(I,1)=0;
    end    
    if wat_flxsa(I) >= 0
			betha1(I,1)=0;
			betha4(I,1)=1;
    end
    
    if  wat_flxsa(I+1) <= 0 
			betha2(I,1)=1;
			betha3(I,1)=0;
    end      
    if wat_flxsa(I+1) > 0
			betha2(I,1)=0;
			betha3(I,1)=1;
    end
end

%-------------------------------------------------------------------------
[cm,cim,csol,acsolmo]= solute_data (cm,cim,csol,acsolmo,t,dx,tcsolo_ini,wco,wcmo,wcma,wcio,immobile,solute_param1,solute_param2,nsol,ncs,pvela);
[conc,solute_applic]=In_solute_boundary(t,nsol,tprec,wat_flxs,sol_upper_bc_type);
[decsolm,decsoli,sol_sinkm,sol_sinki]=solute_sink(dt,nsol,ncs,cm,cim,wcmo,...
    wcio,solute_param1,solute_param2,immobile);
if sim_nitro;
    [minerm,mineri,rcarbmin,soil_om,carbman,carblit,rnitman,rnitlit,cco2o]=...
     miner_immob(t,dt,temp,ncs,WC,csol,ph,soil_om,wco,dx,cm,cim,soil_parameters,om_param, miner_param,plant_uptake_param);
    [uptakem,uptakei,tot_upt]=...
    nitro_uptake_wavefor(csol,dt,dx,harvest_date,isucr, ncs,nsol,soil_parameters,om_param, miner_param,...
    plant_uptake_param,solute_param1,solute_param2,plant_date,simplant,t,rtex,tot_upt,wco);
  
    [sol_sink,sol_sinkm,sol_sinki,nitrifm,nitrifi,hydro_uream,hydro_ureai,uptakem,uptakei,denitm,...
    deniti,volatm,volati,decnorg,deccorg]=...
    nitro_sink(t,dt,tmax,temp,nsol,ncs,cm,cim,WC,wcmo,wcio,decsolm,decsoli,minerm,mineri,rcarbmin,...
    uptakem,uptakei,soil_parameters,solute_param1,solute_param2,om_param, miner_param, plant_uptake_param,csol,immobile);
end

%--------------------------------------------------------------------------
%% actual value tto passed value
csp=cs;
cmp=cm;cimp=cim; % cm and cimp = precedent mobile and immobile soute concentration
acsolmb=acsolmo;
%set cmp zero gradient
cmp(ncs+1,:)=cmp(ncs,:);

if sol_upper_bc_type==1;   %%%'concentration_flux_bc'
    if wat_flxs(1)<0
        cs=conc/abs(wat_flxs(1));
    else
        cs=zeros(1,nsol);
    end
elseif sol_upper_bc_type==2;    %%%'concentration_bc'
    if wat_flxs(1)<=0
        cs=conc;
    else
        cs=zeros(1,nsol);
    end
end

for j=[1:nsol];                     
   
%TOP VALUES
     I=1;
ab(1,1)=  dt*pveloh(1)*wat_flxsah(1)*(wcmo(1)-wcmob(1))/(16*dx^2)+...        
    (diffush(1,j))*((wcmah(1)))/(2*dx^2);
bb(1,1)= dt*pveloh(1)*pvelah(2)*(wcmo(I)-wcmob(I))/(16*dx^2)+...
    (diffush((2),j))*(wcmah(2))/(2*dx^2);

cb(1,1)=wat_flxsah(1)/(2*dx);                           

db(1,1)=wat_flxsah(2)/(2*dx);

bl(1,1)= wcmo(1)/dt +ab(1,1)+bb(1,1)+(bd(1)*f(1)*kd(1,j))/dt-betha2(1)*db(1,1)+betha4(1)*cb(1,1);

cl(1,1)= -bb(1,1)-betha3(I)*db(1,1);
dl(1,1)= csp(j)*(ab(1,1) - betha1(I)*cb(1,1))+ cs(j)*(ab(1,1)-betha1(1)*cb(1,1))+...
    cmp(1,j)*(wcmob(1)/dt -ab(1,1)-bb(1,1)+(bd(I)*f(1)*kd(1,j))/dt+betha2(1)*db(1,1)-betha4(1)*cb(1,1))+...
    cmp((2),j)*(bb(1,1)+betha3(1)*db(1,1))+ sol_sinkm(1,j); 


%INTERMEDIATE VALUES
   I= [2:ncs-1]; % for
       
ab(I,1)= (diffush(I,j))'.*(wcmah(I))./(2*dx^2) +...
dt*pveloh(I).*((pvelah(I))).*(wcmo(I)-wcmob(I))./(16*dx^2);

bb(I,1)= (diffush((I+1),j))'.*wcmah(I+1)./(2*dx^2) +... 
dt*pveloh(I).*((pvelah(I+1))).*(wcmo(I)-wcmob(I))./(16*dx^2);

cb(I,1)=wat_flxsah(I)./(2*dx);

db(I,1)=wat_flxsah(I+1)./(2*dx);

al(I-1,1)= -ab(I,1)+betha1(I).*cb(I,1);

bl(I,1)= (wcmo(I))'./dt +ab(I,1)+bb(I,1)+(bd(I).*f(I).*kd(I,j))/dt-betha2(I).*db(I,1)+betha4(I).*cb(I,1);

cl(I,1)= -bb(I,1)-betha3(I).*db(I,1);

dl(I,1)= cmp((I-1),j).*(ab(I,1) - betha1(I).*cb(I,1))+ ...
    cmp(I,j).*((wcmob(I))'/dt -ab(I,1)-bb(I,1)+(bd(I).*f(I).*kd(I,j))./dt+betha2(I).*db(I,1)-betha4(I).*cb(I,1))+...
    cmp((I+1),j).*(bb(I,1)+betha3(I).*db(I,1))+ sol_sinkm(I,j); 
     
    
 %BOTTOM VALUES
     I=ncs;
ab(ncs,1)= ((diffush(ncs,j)))*(wcmah(ncs))/(2*dx^2) +...
dt*pveloh(ncs)*((pvelah(ncs)))*(wcmo(ncs)-wcmob(ncs))/(16*dx^2);

bb(ncs,1)= 0;  %%%%%%%vmm(ncs+1/2)?

cb(ncs,1)=(wat_flxsah(ncs))/(2*dx);

db(ncs,1)=(wat_flxsah(ncs+1))/(2*dx);         

al(ncs-1,1)= -ab(ncs,1)+betha1(ncs)*cb(ncs,1); 

bl(ncs,1)= wcmo(ncs)/dt +ab(ncs,1)+bb(ncs,1)+(bd(ncs)*f(ncs)*kd(ncs,j))/dt -...
    betha2(ncs)*db(ncs,1)+betha4(ncs)*cb(ncs,1);

dl(ncs,1)=cmp((ncs-1),j)*(ab(ncs,1)-betha1(ncs)*cb(ncs,1)) +...
cmp(ncs,j)*(wcmob(ncs)/dt-ab(ncs,1)+(bd(ncs)*f(ncs)*kd(ncs,j)/dt)+...
    betha2(ncs)*db(ncs,1)-bb(ncs,1)-betha4(ncs)*cb(ncs,1))+...
    cmp(ncs+1,j)*(bb(ncs,1) +betha3(ncs)*db(ncs,1));

 %solution of the tridiagonal matrix

 cm(:,j)=(inv(diag(bl)+diag(cl,1)+diag(al,-1))*dl).';
   
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
    
    if sol_upper_bc_type==1;
        if wat_flxsah(1)<0;
            sflxsa(1,j) = (cs(j)+csp(j))/2*abs(wat_flxs(1))*dt+(diffush(1,j))*...
                ((cs(j)+csp(j))-(cm(1,j)+cmp(1,j)))/2/dx*(wcmah(1))*dt;
            
        else
            sflxsa(1,j)=0.0;
        end
    elseif sol_upper_bc_type==2;
        if wat_flxsah(1)<=0;
            sflxsa(1,j) = (cs(j)+csp(j))/2*abs(wat_flxs(1))*dt+(diffush(1,j))*...
                ((cs(j)+csp(j))-(cm(1,j)+cmp(1,j)))/2/dx*(wcmah(1))*dt;
        else
            sflxsa(1,j)=0.0;
        end
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
%%%
end 

