function [cberr_sol,tflsol,rleasa,tcsink,tsoli,tcsolo_ini,dsol,bot_sol_flux]=sol_intgr(t,dt,dx,nsol,ncs,cm,cim,cmp,cimp,csol,wcmo,wcio,wcob,wco,solute_param1,...
    solute_param2,sol_sinki,sol_sinkm,acsolmo,acsolio,acsolo,sflxsa,tflsol,rleasa,tsoli,tcsink,tcsolo_ini,immobile,bot_sol_flux);
%%
%tsoli= inital mass solute in the profile
%tflsol total solute inflow in the soil surface
%rleasa total solute inflow at the bottom
%tsink= total sink term

bd=solute_param1(:,1);
alf=solute_param1(:,4);
volume=dx;


if t==0
    tsoli=sum(tcsolo_ini);
end

for j = 1:nsol  %%     do sp = 1,nr_of_sol
%		calculate new total values of solute on the nodes
		for i=1:ncs    %%%		do i=1,ncs
			if immobile
				tcsolmo(i,j)=(cm(i,j)*wcmo(i)+(acsolmo(i,j)*bd((i))))*volume;
				tcsolio(i,j)=(cim(i,j)*wcio(i)+(acsolio(i,j)*bd((i))))*volume;
				tcsolo(i,j)=(tcsolmo(i,j)+tcsolio(i,j));
			else
				tcsolo(i,j)=(csol(i,j)*wco(i)+(acsolo(i,j)*	bd((i))))*volume;
				tcsolmo(i,j)=tcsolo(i,j);
				tcsolio(i,j)=0.0;
            end
        end
  
        ptsol=0.0;
		ptsoli=0.0;
		ptsolm=0.0;
		ptsinki=0.0;
		ptsinkm=0.0;
		for i=1:ncs;               
            ptsoli=ptsoli+tcsolio(i,j);
            ptsolm=ptsolm+tcsolmo(i,j);
            ptsol=ptsol+tcsolo(i,j);
            tsinki(i,j) = (sol_sinki(i,j))*volume*dt;
            tsinkm(i,j) = (sol_sinkm(i,j)+alf((i))*(cm(i,j)-cim(i,j)))*volume*dt;
            ptsinki = ptsinki +(sol_sinki(i,j))*volume*dt;
        end
       
        
        ptsinkm= sum(sol_sinkm(:,j)+(alf.*((cm(1:ncs,j)+cmp(1:ncs,j))/2-(cim(:,j)+...
                cimp(:,j))/2)))*volume*dt; 
        
        ptsink=ptsinki+ptsinkm ;
        
        
        
%           calculate cumulative values for the whole profile
%		netto flux
		rleasa(j)=rleasa(j)+sflxsa(ncs+1,j);
       bot_sol_flux{j}(:,end+1)=rleasa(j);
%		cumulative downward flux (negativ)
       %if sflxsa(ncs+1,j)<=0.0 
          %dleasa(sp)=dleasa(j)+sflxsa(ncs+1,j);
       %end
%		cumulativ upward flux (positiv)
       %if sflxsa(ncs+1,j)>=0.0 
          %pleasa(j)=pleasa(j)+sflxsa(ncs+1,j);
       %end
%		inflow
		tflsol(j) = tflsol(j)+sflxsa(1,j);
%		inflow during one day
	%	if(t- dint(t).eq.0.) solinfl(sp)=0.0
	%	solinfl(sp)=solinfl(sp)+sflxsa(1,sp)
%		sink first order
		tcsink(j) = tcsink(j) + ptsink;
%		change in the system from the start till the present time step
          dsol(j)=ptsol-tsoli(j);
%		mass balance error 
		cberr_sol(j)=dsol(j) - (rleasa(j)+tflsol(j)+tcsink(j));
end

                 
        
        
        