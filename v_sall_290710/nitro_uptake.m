function[uptakem,uptakei,tot_upt]=nitro_uptake(csol,dt,dx,harvest_date,...
        isucr, ncs,nsol,soil_parameters,om_param, miner_param, plant_uptake_param, solute_param1,solute_param2,...
        plant_date,simplant,t,rtex,tot_upt,wco);

 
        

%%%%%%%%%%%%%%%%%%%%%%
%NITROGEN MODULE
%nitro_uptake calculate the plant N uptake
%IN:
%t,dt,tmax= time, time increment and maximum time(min)
%nsol, ncs=number of solute spieces and number of compartiments
%WC:water content at time t, wcp: water content at time t-dt
%cmp:solute mobile concentration 
%wco= water content
%soil_parameters= parameters of the soil
%solute_param1,solute_param2 = parameters of the solute (from scrip%in_soliue_parameters)
%nitro_param = nitrogen parameters (from scrip in_nitro_om_param)
%rtex= root water extraction calculated in the water flow module

%:
%tot_upt= cumulative N uptake
%uptakem and uptakei = N uptake in the mobile and immobile region
%CALLED BY solve_solute
%--------------------------------------------------------------------------
%M.SALL 17/03/09


ratio=solute_param1(:,3);
dif=solute_param2(:,2*nsol+1:3*nsol);
ar=solute_param2(:,3*nsol+1:4*nsol);
br=solute_param2(:,4*nsol+1:5*nsol);
rorad=plant_uptake_param(1);       %%rorad=nitro_param(4);
rdo=plant_uptake_param(2);         %%rdo=nitro_param(5);
g= plant_uptake_param(3);           %%%g= nitro_param(13);  
rnmaxp=plant_uptake_param(4);  %%rnmaxp=nitro_param(1);
rdens0=plant_uptake_param(5);   %%rdens0=nitro_param(2);
alfa_rdens=plant_uptake_param(6);  %%alfa_rdens=nitro_param(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%        


isucr=0; %%%%%PROVISOIRE MAMADOU

rlncl=0.005;
rcsolo=[0 1.8 6.2];

volume=dx;
if simplant ==0 | t< plant_date | t> harvest_date;
    for sp=2:3;
        for i=1:ncs;
        uptakem(i,sp)=0;
        uptakei(i,sp)=0;
        end
    end
end

if simplant & t>=plant_date & t<harvest_date;
    if isucr;
        anlv = sucr_anlv();
    else
        totdem = rnmaxp/100;
        g = g * ((harvest_date - plant_date)+1);
        fact = rnmaxp/(((g/2) - (g/3))*(g^2));
    end
    
    if isucr;
        rdens= rdens_sucros;
    else
    %    rdens=ones (ncs)*rdens0;     %%% A CULCULER ULTERIEUREMENT ();
        rdens= [ones(20,1)*0.05; zeros(80,1)];
    end

%if t == harvest_date
%	rfract_tot = 0;
%	for i= 1:ncs;
%	rfract_tot=rfract_tot+rdens(i);
%	rdens_left(i) = rdens(i);
 %   end;
%	rfract_tot = max(rfract_tot);
%end;
 
    if isucr;
	xncle= calc_xncle();
	xncst=calc_xncst();
	xncrt=calc_xncrt();
	rmncl=0.5*xncle;
            if rmncl == 0;
			error ('rmncl is zero');
            end
	 plant_weights(wso,wlv,wrt,wst,wlvg);
	ancl=anlv/wlvg;
    red = calc_nitreductgrow();
	fndef=1-sqrt(1-red*red);
	rndeml= dmax1(0.d0,(wlv*xncle-anlv));
	rndemst= dmax1(0.d0,(wst*xncst-anst));
	rndemrt= dmax1(0.d0,(wrt*xncrt-anrt));
	xncso=calc_xncso();
	rndemso= min(max(0,(wso*xncso-anso))*fndef,rndeml+rndemst);
	totdem= max(0.d0,rndeml+rndemst+rndemrt)*100;
        if totdem == 0;
		uptakem(1:ncs,1:nsol)=0;
		uptakei(1:ncs,1:nsol)=0;
        end
    else
        if(rnmaxp-tot_upt)< 0;
        uptakem(1:ncs,1:nsol)=0;
        uptakei(1:ncs,1:nsol)=0;
        end

	time_since_pl = t-(plant_date)+1;
	totdem= max(fact*time_since_pl*(g-time_since_pl),0)* dt;
        if rnmaxp < tot_upt + totdem ;
       totdem = rnmaxp - tot_upt;
        end
    end
tconn=0;
    for sp=2:3;
        for i=1:ncs;
    hcsolo(sp) = max(0,(csol(i,sp)-rcsolo(sp)));
    unc(i,sp) = rtex(i)*hcsolo(sp).*wco(i)*dt*dx;
    tconn=tconn+(unc(i,sp));
        end
    end
		
    if tconn > totdem;
   unc(1:ncs,1:nsol)=unc(1:ncs,1:nsol)*totdem/tconn;
   und(1:ncs,1:nsol) =0;
   tconn=totdem;
   tund=0;
    else
   tpdnup=totdem-tconn;
   tund=0;
    for sp=2:3;
        for i=1:ncs
       hcsolo(sp) = max(0,(csol(i,sp)-rcsolo(sp)));
	   diffus_rm(i,sp) = dif(i,sp).*ar(i,sp).*exp(br(i,sp).*wco(i))/wco(i);
	   und(i,sp)=rdens(i)*rorad*2*pi*diffus_rm(i,sp)*hcsolo(sp)*wco(i)*dx*dt/(rdo*100000);
	   tund=tund+und(i,sp);
        end
    end
    if tund > tpdnup;
        for sp = 2,3;
   	 und(1:ncs,sp)=und(1:ncs,sp)*tpdnup/tund;
        end
     tund=tpdnup;
    end
    end
%calculate the total uptake (mg/m**2)
totvr=(tconn+tund);
tot_upt=tot_upt+totvr;
    if isucr;
    rnupl=totvr*rndeml/(totdem+0.000000001);
    rnupst=totvr*rndemst/(totdem+0.0000000001);
	rnuprt=totvr*rndemrt/(totdem+0.0000000001);
%calculate new concentrations in the different parts
	anlv=anlv+(rnupl-rndemso*anlv/(anlv+anst));
	anst=anst+(rnupst-rndemso*anst/(anlv+anst));
	anrt=anrt+rnuprt;
	anso=anso+rndemso;
    end
%calculate the uptake rates
    for sp=2:3;
        for i=1:ncs;
        uptakem(i,sp)=(unc(i,sp)+und(i,sp)).*ratio(i)/(dx*dt);
        uptakei(i,sp)=(unc(i,sp)+und(i,sp)).*(1-ratio(i))/(dx*dt);
        end
    end

end
