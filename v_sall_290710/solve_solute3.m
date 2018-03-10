function[cs,cm,cim,cmp,cimp,csol,c_manp,n_manp,c_litp,n_litp,c_hump,n_hump,sol_sinkm,sol_sinki,tot_upt,acsolo,acsolio,acsolmo,sflxsa,solute_applic]=...
    solve_solute(t,dt,tmax,temp,dx,nsol,ncomp,WC,rtex,cm,cim,csol,cs,wat_flxs,ph,...
    wat_flxsa, wat_flxsah,wcio, wciob, wcma, wcmah, wcmo, wcmob,wco,wcob,pvela,pvelah,...
    pveloh,c_manp,n_manp,c_litp,n_litp,c_hump,n_hump,immobile,soil_parameters,solute_param1,solute_param2,...
    nitro_param,om_param,tot_upt,sim_nitro,tcsolo_ini,diffush);
    
    
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
lm(ncomp+1)=lm(ncomp);
f=(solute_param1(:,5));
rates=(solute_param2(:,1:nsol));
kd=(solute_param2(:,nsol+1:2*nsol));
dif=solute_param2(:,2*nsol+1:3*nsol);
ar=solute_param2(:,3*nsol+1:4*nsol);
br=solute_param2(:,4*nsol+1:5*nsol);


%calculate beta coeficients 
for I = 1:ncomp;
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
[cm,cim,csol]=solute_data (cm,cim,csol,t,dx,tcsolo_ini,wco,wcmo,wcma,wcio,immobile,solute_param1,solute_param2,nsol,ncomp,pvela);
[conc,solute_applic]=solute_boundary(t,nsol);
[decsolm,decsoli,sol_sinkm,sol_sinki]=solute_sink(dt,nsol,ncomp,cm,cim,wcmo,...
    wcio,solute_param1,solute_param2,immobile);
if sim_nitro;
    [minerm,mineri,rcarbmin,c_manp,n_manp,c_litp,n_litp,c_hump,n_hump]=...
     miner_immob(t,dt,temp,ncomp,WC,csol,ph,c_manp,n_manp,...
     c_litp,n_litp,c_hump,n_hump,wco,dx,cm,cim,soil_parameters,nitro_param,om_param);

    [uptakem,uptakei,tot_upt] = nitro_uptake(t,dt,tmax,dx,ncomp,nsol,wco,cm,...
     rtex,tot_upt,soil_parameters,solute_param1,solute_param2,nitro_param);
  
    [sol_sink,sol_sinkm,sol_sinki]=nitro_sink(t,dt,tmax,temp,nsol,ncomp,cm,cim,...
    WC,wcmo,wcio,decsolm,sol_sinki,minerm,mineri,rcarbmin,uptakem,uptakei,...
    soil_parameters,solute_param1,solute_param2,nitro_param,om_param);
end

%--------------------------------------------------------------------------

csp=cs;
cmp=cm;cimp=cim; % cm and cimp = precedent mobile and immobile soute concentration
%set cmp zero gradient
cmp(ncomp+1,:)=cmp(ncomp,:);

if wat_flxs(1)<=0
    cs=conc;  
else
    cs=zeros(1,nsol);
end


for j=[1:nsol];                         
    
%TOP VALUES
     I=1;
ab(1,1)= (diffush(1,j))*((wcmah(1)))/(2*dx^2);%+...
          
    
bb(1,1)= (diffush((2),j))*(wcmah(2))/(2*dx^2);%+...
         

cb(1,1)=wat_flxsah(1)/(2*dx)+...
    dt*pveloh(1)*wat_flxsah(1)*(wcmo(1)-wcmob(1))/(16*dx^2); %cb(1,1)=4*cb(1,1);
    
db(1,1)=wat_flxsah(2)/(2*dx)+...
dt*pveloh(1)*pvelah(2)*(wcmo(I)-wcmob(I))/(16*dx^2);%db(1,1)=4*db(1,1);

bl(1,1)= wcmo(1)/dt +ab(1,1)+bb(1,1)+(bd(1)*f(1)*kd(1,j))/dt-betha2(1)*db(1,1)+betha4(1)*cb(1,1);

cl(1,1)= -bb(1,1)-betha3(I)*db(1,1);
dl(1,1)= csp(j)*(ab(1,1) - betha1(I)*cb(1,1))+ cs(j)*(ab(1,1)-betha1(1)*cb(1,1))+...
    cmp(1,j)*(wcmob(1)/dt -ab(1,1)-bb(1,1)+(bd(I)*f(1)*kd(1,j))/dt+betha2(1)*db(1,1)-betha4(1)*cb(1,1))+...
    cmp(2,j)*(bb(1,1)+betha3(1)*db(1,1))+ sol_sinkm(1,j); 


%INTERMEDIATE VALUES
   I= [2:ncomp-1]; % for
       
ab(I,1)= (diffush(I,j))'.*(wcmah(I))./(2*dx^2);%+...

bb(I,1)= (diffush((I+1),j))'.*wcmah(I+1)./(2*dx^2);%+... 
 

cb(I,1)=wat_flxsah(I)./(2*dx)+...
dt*pveloh(I).*((pvelah(I))).*(wcmo(I)-wcmob(I))./(16*dx^2);%cb(I,1)=4*cb(I,1); 

db(I,1)=wat_flxsah(I+1)./(2*dx)+...
   dt*pveloh(I).*((pvelah(I+1))).*(wcmo(I)-wcmob(I))./(16*dx^2);%db(I,1)=4*db(I,1);

al(I-1,1)= -ab(I,1)+betha1(I).*cb(I,1);

bl(I,1)= (wcmo(I))'./dt +ab(I,1)+bb(I,1)+(bd(I).*f(I).*kd(I,j))/dt-betha2(I).*db(I,1)+betha4(I).*cb(I,1);

cl(I,1)= -bb(I,1)-betha3(I).*db(I,1);

dl(I,1)= cmp((I-1),j).*(ab(I,1) - betha1(I).*cb(I,1))+ ...
    cmp(I,j).*((wcmob(I))'/dt -ab(I,1)-bb(I,1)+(bd(I).*f(I).*kd(I,j))./dt+betha2(I).*db(I,1)-betha4(I).*cb(I,1))+...
    cmp((I+1),j).*(bb(I,1)+betha3(I).*db(I,1))+ sol_sinkm(I,j); 
     
    
 %BOTTOM VALUES
     I=ncomp;
ab(ncomp,1)= ((diffush(ncomp,j)))*(wcmah(ncomp))/(2*dx^2)+...
 dt*pveloh(ncomp)*((pvelah(ncomp)))*(wcmo(ncomp)-wcmob(ncomp))/(16*dx^2);%ab(ncomp,1)=0;

bb(ncomp,1)= 0;  %%%%%%%vmm(ncomp+1/2)?

cb(ncomp,1)=(wat_flxsah(ncomp))/(2*dx);
    

db(ncomp,1)=(wat_flxsah(ncomp+1))/(2*dx);         

al(ncomp-1,1)= -ab(ncomp,1)+betha1(ncomp)*cb(ncomp,1); 

bl(ncomp,1)= wcmo(ncomp)/dt +ab(ncomp,1)+bb(ncomp,1)+(bd(ncomp)*f(ncomp)*kd(ncomp,j))/dt -...
    betha2(ncomp)*db(ncomp,1)+betha4(ncomp)*cb(ncomp,1);

dl(ncomp,1)=cmp((ncomp-1),j)*(ab(ncomp,1)-betha1(ncomp)*cb(ncomp,1)) +...
cmp(ncomp,j)*(wcmob(ncomp)/dt-ab(ncomp,1)+(bd(ncomp)*f(ncomp)*kd(ncomp,j)/dt)+...
    betha2(ncomp)*db(ncomp,1)-bb(ncomp,1)-betha4(ncomp)*cb(ncomp,1))+...
    cmp(ncomp+1,j)*(bb(ncomp,1) +betha3(ncomp)*db(ncomp,1));

 %solution of the tridiagonal matrix

 cm(:,j)=(inv(diag(bl)+diag(cl,1)+diag(al,-1))*dl).'; %%
   cm(:,j)=1.01*cm(:,j);
 %%%%    
    for I=1:ncomp;
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
    
    for i=2:ncomp;   
		if wat_flxsah(i)<0  
				sflxsa(i,j)=(wat_flxsah(i)*(cm((i-1),j)+... 
     				cmp((i-1),j))/2)*dt;                
		else
				sflxsa(i,j)=(wat_flxsah(i)*(cm(i,j)+...  
     				cmp(i,j))/2)*dt;   
        end
    end
	if wat_flxsah(1)<0;
			sflxsa(1,j) = (cs(j)+csp(j))/2*abs(wat_flxs(1))*dt+(diffush(1,j))*...
                ((cs(j)+csp(j))-(cm(1,j)+cmp(1,j)))/2/dx*(wcmah(1))*dt;
	else
			sflxsa(1,j)=0.0;
    end
		volume= dx;
	if  abs(cm(ncomp,j)-cmp(ncomp,j)>0.001  & wat_flxsah(ncomp+1)>0.00000000001)
			sflxsa((ncomp+1),j)=((-cmp(ncomp,j)+cm(ncomp,j))*...
     		wcmo(ncomp)*dx)-(acsolmb(ncomp,j)-acsolmo(ncomp,j))*...
     		bd(1)*volume;
    elseif wat_flxsah(ncomp+1)<0.0
				sflxsa((ncomp+1),j)=(wat_flxsah(ncomp+1)*...
     			((cm(ncomp,j)+ cmp(ncomp,j))/2))*dt;
	else
				sflxsa(ncomp+1,j)=(wat_flxsah(ncomp+1)*cmp((ncomp+1),j))*dt;
    end
%%%
end 

