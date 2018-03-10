function[uptakem,uptakei,tot_upt]=nitro_uptake(csol,dt,dx,harvest_date,...
        isucr, ncs,nsol,soil_parameters,om_param, miner_param, plant_uptake_param, solute_param1,solute_param2,...
        plant_date,simplant,t,rtex,tot_upt,wco);


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

totup=tot_upt;
volume=dx;
rcsolo=[0 1.8 6.2];

if simplant ==0 | t< plant_date | t> harvest_date;
    for i=1:ncs
        for sp = 2:3
            sunm(i,sp)=0.00;
            suni(i,sp)=0.00;
        end
    end
else
    totdem = rnmaxp/100;
    g = g * ((harvest_date - plant_date)); %%%+1.0);
    factor = rnmaxp/(((g/2.0) - (g/3.0))*(g^2));
    
    
    rdens= [ones(20,1)*0.05; zeros(80,1)];
    
    %if ( dint(t).eq.t.and.idint(t).eq.harvest_date-1)
    %    rfract_tot = 0.0;
    %    for i= 1,ncs
    %        rfract_tot=rfract_tot+rdens(i);
    %        rdens_left(i) = rdens(i);
    %    end
    %    rfract_tot = max(rfract_tot);
    %end
    
    if (rnmaxp-totup)< 0.1
        for i=1:ncs
            for sp = 2:3
                sunm(i,sp)=0.00;
                suni(i,sp)=0.00;
            end
        end
    end
    
    time_since_pl = t-(plant_date); %+1.0;
    totdem= max(factor*time_since_pl*(g-time_since_pl),0.0)* dt;
    if rnmaxp < (totup + totdem )
        totdem = rnmaxp - totup;
    end
    tconn=0.0;
    for sp=2:3
        for i=1:ncs
            hcsolo(sp) = max(0.0,(csol(i,sp)-0.001*rcsolo(sp)));
            unc(i,sp) = rtex(i)*hcsolo(sp)*wco(i)*dt*dx;
            tconn=tconn+unc(i,sp);
        end
    end
    
    if tconn > totdem
        for sp=2:3
            for j=1: ncs
                unc(j,sp)=unc(j,sp)*totdem/tconn;
                und(j,sp) =0.0;
            end
        end
        tconn=totdem;
        tund=0.0;
    else
        tpdnup=totdem-tconn;
        tund=0.0;
        for sp=2:3
            for i=1:ncs
                hcsolo(sp) = max(0.0,(csol(i,sp)-0.001*rcsolo(sp)));
                %        diffus_rm(i,sp) = parasol(4,nla(i),sp)* dexp(parasol(5,nla(i),sp)*wco(i))*parasol(3,nla(i),sp)/ wco(i);
                diffus_rm(i,sp) = dif(i,sp).*ar(i,sp).*exp(br(i,sp).*wco(i))/wco(i);
                und(i,sp)=rdens(i)*rorad*2*pi*diffus_rm(i,sp)*hcsolo(sp)*wco(i)*dx*dt/(rdo); %%%*100000.0)suprimé car les unités sonr harmonisé
                tund=tund+und(i,sp);
            end
        end
        if tund > tpdnup
            for sp = 2:3
                for j=1: ncs
                    und(j,sp)=und(j,sp)*tpdnup/tund;
                end
            end
            tund=tpdnup;
        end
    end
    % calculate the total uptake (mg/m**2)
    totvr=(tconn+tund);
    totup=totup+totvr;
    % calculate the uptake rates
    for sp=2:3
        for i=1:ncs
            sunm(i,sp)=(unc(i,sp)+und(i,sp))*ratio(i)/(dx*dt);
            suni(i,sp)=(unc(i,sp)+und(i,sp))*(1-ratio(i))/(dx*dt);
        end
    end
end
tot_upt=totup;
uptakem = sunm;
uptakei=suni;
