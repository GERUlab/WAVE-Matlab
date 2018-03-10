function [cberr_sol,tflsol,rleasa,tcsink,tsoli,tcsolo_ini,dsol,tcsolo,initsol]=sol_intgr(t,dt,dx,nsol,ncs,cm,cim,cmp,cimp,csol,wcmo,wcio,wcob,wco,solute_param1,...
    solute_param2,sol_sinki,sol_sinkm,acsolmo,acsolio,acsolo,sflxsa,tflsol,rleasa,tsoli,tcsink,tcsolo_ini,immobile,initsol);

%Calculates the mass balance for solute

%tsoli= inital mass solute in the profile
%tflsol total solute inflow in the soil surface
%rleasa total solute inflow at the bottom
%tsink= total sink term

bd=solute_param1(:,1);
alf=solute_param1(:,4);
volume=dx;
if initsol == 2
    joachim = 1;
    initsol = 0;
else
    joachim = 0;
end
if initsol==0
    tsoli=sum(tcsolo_ini);
    initsol = 1;
end
        ptsol=zeros(1,nsol);
		ptsoli=zeros(1,nsol);
		ptsolm=zeros(1,nsol);
		ptsinki=zeros(1,nsol);
		ptsinkm=zeros(1,nsol);
    for j = 1:nsol  %%     do sp = 1,nr_of_sol
%		calculate new total values of solute on the nodes
		for i=1:ncs    %%%		do i=1,ncs
			if immobile
				tcsolmo(i,j)=(cm(i,j)*wcmo(i)+(acsolmo(i,j)*...
     			bd((i))))*volume;
				tcsolio(i,j)=(cim(i,j)*wcio(i)+(acsolio(i,j)*...
     			bd((i))))*volume;
				tcsolo(i,j)=(tcsolmo(i,j)+tcsolio(i,j));
			else
				tcsolo(i,j)=(csol(i,j)*wco(i)+(acsolo(i,j)*...
     			bd((i))))*volume;
				tcsolmo(i,j)=tcsolo(i,j);
				tcsolio(i,j)=0.0;
            end
        end
  

		for i=1:ncs;               
            ptsoli(j)=ptsoli(j)+tcsolio(i,j);
            ptsolm(j)=ptsolm(j)+tcsolmo(i,j);
            ptsol(j)=ptsol(j)+tcsolo(i,j);
            tsinki(i,j) = (sol_sinki(i,j))*volume*dt;
            tsinkm(i,j) = (sol_sinkm(i,j)+alf((i))*...
            (cm(i,j)-cim(i,j)))*volume*dt;
            ptsinki(j) = ptsinki(j) +(sol_sinki(i,j))*volume*dt;

        end
       
        
        ptsinkm(j)= sum(sol_sinkm(:,j)+(alf.*((cm(1:ncs,j)+cmp(1:ncs,j))/2-(cim(:,j)+...
                cimp(:,j))/2)))*volume*dt; 
        
        ptsink(j)=ptsinki(j)+ptsinkm(j) ;
        
        
        
%      calculate cumulative values for the whole profile
%		netto flux
		rleasa(j)=rleasa(j)+sflxsa(ncs+1,j);
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
		tcsink(j) = tcsink(j) + ptsink(j);
%		change in the system from the start till the present time step
          dsol(j)=ptsol(j)-tsoli(j);
%		mass balance error 
		cberr_sol(j)=dsol(j) - (rleasa(j)+tflsol(j)+tcsink(j));
    end

 if  joachim == 1
    initsol =2;
end      
       
        
        
        