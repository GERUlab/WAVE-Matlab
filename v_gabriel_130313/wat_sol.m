function[wat_flxsa, wat_flxsah, pvela, pvelah, pveloh,pvelo, wcio, wciob, wcma, wcmah,wcmob, wcmo, wco, wcob,diffus,diffush,initsol]=...
         wat_sol(t,wat_flxs, wat_flxsa, immobile, ncs,solute_param1,solute_param2,pvela,WC, wcio, wcma,wcmo, wco,pvelo,diffus,nsol,dx,initsol,solboco_top_type);

%Prepare data from the water flow module  for the solute module     
%mainly inspired from WATSOL.FOR         
%IN
%t:time
%wat_flxs:water fluxs 
%immobile: choice for immobile water concept (see in_solute_parameter)
%solute_param1: parameters of the solutes
%solute_param2: parameters of the solutes
%ncs: number of compartiments
%WC: water content
%nsol:number of solutes
 
%OUT
%wat_flxsa : water_flxs between the nodes
%wat_flxsah:water_flxs between the nodes at half time
%wcio: water content in the immobile region 
%wciob:water content in the immobile region at time t-dt
%wcma:water content in the mobile region between nodes 
%wcmah: water content in the mobile region between nodes at half time
%wcmo: water content in the mobile region at the nodes 
%wcmob:water content in the mobile region at the nodes at time t-dt
%wco: water content at the nodes
%pvela: pore velocity between the nodes
%pvelah: pore velocity between the nodes at half time
%pveloh: pore velocity at the nodes at half time
%pvelo: pore velocity between the nodes 
%diffush: apparent diffusivity
%CALL: none
%CALLED BY: solve_solute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M. Sall 10/04/09

ar=solute_param2(:,3*nsol+1:4*nsol);
br=solute_param2(:,4*nsol+1:5*nsol);
dif=solute_param2(:,2*nsol+1:3*nsol);
lm=solute_param1(:,2);
lm(ncs+1)=lm(ncs);
ratio=(solute_param1(:,3))';
if initsol == 2
    joachim = 1;
    initsol = 0;
else
    joachim = 0;
end
if initsol == 0 %first time in this program
%if t==0 

%		calculation of the fluxes and water contents across the nodes
 
		wat_flxsa=wat_flxs;  
	if(immobile) 		
		wcmo=ratio.*WC;
		wcio=WC-wcmo;
    else          
		wcmo=WC; 
		wcio=zeros(1,ncs);
    end

    if wat_flxs(1) <=0      
		if immobile   
            wcma(1)=0.999;
		else
			wcma(1)=1.0;
        end
	else 
		if immobile   
			wcma(1)= ratio(1).*WC(1);  
		else
			wcma(1)=WC(1);
        end
    end
	if immobile   
        wcma(1,2:ncs)=(wcmo(1,2:ncs)+wcmo(1,1:ncs-1))/2;
        wcma(ncs+1)=wcma(ncs);
    else  
		wcma(1,2:ncs)=WC(1,2:ncs);
		wcma(ncs+1)=wcma(ncs);
    end
		wco(1,1:ncs)=WC(1,1:ncs);

%		calculation of the pore velocity across the nodes
 
	pvela(1,1:ncs+1)=wat_flxsa(1,1:ncs+1)./wcma(1,1:ncs+1);
	wat_flxso(1,1:ncs)=(wat_flxsa(1,1:ncs)+wat_flxsa(1,2:ncs+1))/2; 
	pvelo(1,1:ncs)=wat_flxso(1,1:ncs)./wcmo(1,1:ncs);       
    
    dm=lm.*abs(pvela');
   
    for j=1:nsol
        de(1:ncs,j)=(dif(1:ncs,j)).*ar(1:ncs,j).*exp(br(1:ncs,j).*(wcmo(1,1:ncs))'); 
        diffus(1:ncs,j)=de(1:ncs,j)./(wcma(1,1:ncs))'+dm(1:ncs);
    end
    
end

%%    calculation of the fluxes and water contents across the nodes

wat_flxsab(1,1:ncs+1) = wat_flxsa(1,1:ncs+1);  
wat_flxsa(1,1:ncs+1)=wat_flxs(1,1:ncs+1);  
wcmob(1,1:ncs) = wcmo(1,1:ncs);
wciob(1,1:ncs) = wcio(1,1:ncs);

if immobile    
	wcmo(1,1:ncs)=ratio(1,1:ncs).*WC(1,1:ncs);
	wcio(1,1:ncs)=WC(1,1:ncs)-wcmo(1,1:ncs);
else
	wcmo(1,1:ncs)=WC(1,1:ncs);
	wcio=zeros(1,ncs);
end
	wcmab(1,1:ncs+1) = wcma(1,1:ncs+1);

if wat_flxsa(1)<=0   
	if immobile 
		wcma(1)=0.999;
	else
		wcma(1)=1.0;
    end
else 
	if immobile  
		wcma(1)=ratio(1)*WC(1);  
	else
		wcma(1)=WC(1);
    end
end

if immobile 
	wcma(1,2:ncs)=(wcmo(1,2:ncs)+wcmo(1,1:ncs-1))/2;
	wcma(ncs+1)=wcma(ncs);
else
	wcma(1,2:ncs)=WC(1,2:ncs);
	wcma(ncs+1)=wcma(ncs);
end
 
    wcob(1,1:ncs) = wco(1,1:ncs);
    wco(1,1:ncs)=WC(1,1:ncs);

%%     calculation of the pore velocity across the nodes 
 
pvelab(1,1:ncs+1) = pvela(1,1:ncs+1);
pvela(1,1:ncs+1)=wat_flxsa(1,1:ncs+1)./wcma(1,1:ncs+1); 
pvelob(1,1:ncs) = pvelo(1,1:ncs);
wat_flxso(1,1:ncs)=(wat_flxsa(1,1:ncs)+wat_flxsa(1,2:ncs+1))/2;  
pvelo(1,1:ncs)=wat_flxso(1,1:ncs)./wcmo(1,1:ncs);

%%     calculation of the fluxes and pore velocity across the nodes and   half ime

wat_flxsah(1,1:ncs+1)=(wat_flxsab(1,1:ncs+1)+wat_flxsa(1,1:ncs+1))/2;
wcmah(1,1:ncs+1)=(wcmab(1,1:ncs+1)+wcma(1,1:ncs+1))/2;
pvelah(1,1:ncs+1)=(pvelab(1,1:ncs+1)+pvela(1,1:ncs+1))/2;
pveloh(1,1:ncs)=(pvelob(1,1:ncs)+pvelo(1,1:ncs))/2;

%%     set the upper boundary condition
if wat_flxsa(1)>0 
	wat_flxsa(1)=0.0;
	wat_flxsab(1)=0.0;
	wat_flxsah(1)=0.0;
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%calculation of the apparent diffsivity and choice of the solute upper BC type

dm=lm.*abs(pvela');
%%%%%%%%%%%%%%%%%%%%
% if lm<=1.0
%    dm=max(0,(lm-0.139*dx)); % correction for numerical dispersion
% end
%dm(1:ncs) = max(0,(lm(1:ncs)-0.139*dx).*abs(pvela(1:ncs)'));
%%%%%%%%%%%%%%%%%%%%
diffusb=diffus;    % apparent diffusion coefficient at t-dt 
for j=1:nsol
    de(1:ncs,j)=(dif(1:ncs,j)).*ar(1:ncs,j).*exp(br(1:ncs,j).*(wcmo(1,1:ncs))'); % effective diffusion coefficient at time t
    diffus(1:ncs,j)=de(1:ncs,j)./(wcma(1,1:ncs))'+dm(1:ncs); 
    diffush(1:ncs,j)=(diffusb(1:ncs,j)+diffus(1:ncs,j))/2;            %%%%%%%%%% apparent diffusion coefficient at half time

%    [conc,solute_applic,solboco_top_type,solboco_input_type]=In_solute_boundary(t,nsol);
    if solboco_top_type == 2
    diffush(1,j)=0;   %%% set to 0 for flux BC or delete this line for concentration BC
    end
end

if  joachim == 1
    initsol =2;
end


