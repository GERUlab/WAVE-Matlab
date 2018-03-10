function [dt,hs,phsa,pot_surface_flux,flxa1,boco_top_type,boco_top,boco_bot_type,...
    boco_bot,pond,bctop_changed,bcbot_changed,first_time_bc,dt_changed,phsurf,esa,BOUNDARY_CONDITIONS_MATRIX]=...
    soil_boundary_conditions(var_changed,t,ph,dt,phsurf,pot_surface_flux,...
   	boco_top_type,boco_top,boco_bot_type,boco_bot,pond,runoff,simplant,dt_min,dt_start,...
    BOUNDARY_CONDITIONS_MATRIX,esa,hs,phsa);

% Calculate/find the current boundary conditions data
%
% IN :
%   General characteristics:t,ph,dt,pot_surface_flux,simplant
%   Characteristics of the bc in the last iteration: boco_top_type,boco_top,boco_bot_type,boco_bot,pond,runoff,phsurf
%   Changes during last iteration: var_changed (=bctop_changed,bcbot_changed,dt_changed_bc,first_time_bc)
% OUT :
%   dt,hs,phsa,pot_surface_flux,flxa1,boco_top_type,boco_top,boco_bot_type,
%   boco_bot,pond,bctop_changed,bcbot_changed,first_time_bc,dt_changed
% CALLS:
%   In_Boundary_conditions.m,In_ETsplit.m
% CALLED BY:
%   wavemat105.m
%------------------------------------------------
% Javaux M., Lambot S. & Vanclooster M. (2000)

%definitions of variables
bctop_changed=var_changed(1);
bcbot_changed=var_changed(2);
dt_changed=var_changed(3);
first_time_bc=var_changed(4);

%initialization
if first_time_bc
   dt_changed=0;
   first_time_bc=0;
end
%dt_changed,first_time_bc
% Searches for the appropriate boundary condition for the current time t 
i=max(find(BOUNDARY_CONDITIONS_MATRIX(:,1)<=t));
i2=max(find(BOUNDARY_CONDITIONS_MATRIX(:,1)<=(t+dt)));% at time t+dt

bctop_changed=0;
% if ponding appears
info_pond=0;

if dt_changed==1 & i==i2 
   dt_changed=0;
   if(BOUNDARY_CONDITIONS_MATRIX(i,2)~=boco_top_type)
   % dt has been optimized in the previous t to fall just in the change time
   % now, new user BC =>changes in the boundary condition
   if (boco_top_type==1 & boco_top>0)
      %bocotop is ponding=> the imposed flux is recorded
      imposed_flux=BOUNDARY_CONDITIONS_MATRIX(i,3);
      info_pond=1;
      %next BC will be ponding
      BOUNDARY_CONDITIONS_MATRIX(i,2)=1;
      BOUNDARY_CONDITIONS_MATRIX(i,3)=boco_top;      
      %if no changes of boco_top_type but in boco_top
      %=>Only user boundary condition is taken in account
   end
	%changes of boco_bot_type is prioritary on actual BC=> no adaptations are needed
	%reinitialization of dt_changed
   dt_changed=0;
   end
else
   %3 possibilities:
   %	a)dt_changed==1 & i~=i2 : dt has been reduced but re-changed in calc_dt.m.
   %	b)dt_changed==0 & i~=i2 : changes in future dt: reduce dt to stay in the same UBC
   %	c)dt_changed==0 & i==i2 : no changes in the current or future User BC
   if bctop_changed==1
      %Current BC are different of user BC
      if BOUNDARY_CONDITIONS_MATRIX(i,2)==2			
         % if user BC is flux (current is ph), then define the imposed flux
         imposed_flux=BOUNDARY_CONDITIONS_MATRIX(i,3);
         info_pond=1;
      end
      %definition of the future BC=current BC
      BOUNDARY_CONDITIONS_MATRIX(i,2)=boco_top_type;
      if boco_top_type==1									
         % new TBC is pression/succion
         BOUNDARY_CONDITIONS_MATRIX(i,3)=phsurf;
      else														
         % new TBC is flux
      	BOUNDARY_CONDITIONS_MATRIX(i,3)=pot_surface_flux ; 
      end
   end
   %check if the dt has to be reduce
   if i~=i2
      %dt encompass two or more changes of BC
      if i~=i2-1
         i2=1+i;
		end
      % reduction of the time step
      dt=BOUNDARY_CONDITIONS_MATRIX(i2,1)-t;
      dt_changed=1;
      %if the new dt is < than dtmin, do the change in this iteration
      if dt < dt_min
         i=i2;
         dt_changed=0;
         dt=dt_start;
      end
   else
      dt_changed=0;
   end
end

%Next Boundary condition for time t
   boco_top_type=BOUNDARY_CONDITIONS_MATRIX(i,2) ;
   boco_top=BOUNDARY_CONDITIONS_MATRIX(i,3);
   boco_bot_type=BOUNDARY_CONDITIONS_MATRIX(i,4);
   boco_bot=BOUNDARY_CONDITIONS_MATRIX(i,5);


% Upper boundary condition
if boco_top_type==1
   %succion/pression
   flxa1=esa;
   phsurf=boco_top;
elseif boco_top_type==2										
   % flux, pluie
   flxa1=esa+boco_top; 
end
if info_pond==1
   %ponding appears
   flxa1=esa+imposed_flux;
end

%potential flux on the top (=flxar)
pot_surface_flux=flxa1-pond/dt;
%that is the water which is imposed by the user and the water in the ponding

%if runoff, pot surf flux is reducted
if runoff~=0
   pot_surface_flux=pot_surface_flux-runoff;
end