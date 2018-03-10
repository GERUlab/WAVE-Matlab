function [hs,phsa,pot_surface_flux,flxa1,boco_top_type,boco_top,boco_bot_type,...
    boco_bot,pond,bctop_changed,bcbot_changed,first_time_bc,dt_changed_bc,next_dt_new_bc,...
    phsurf,esa,BOUNDARY_CONDITIONS_MATRIX,ncbot]=...
    soil_boundary_conditions(var_changed,t,ph,dt,phsurf,pot_surface_flux,...
   	boco_top_type,boco_top,boco_bot_type,boco_bot,pond,runoff,simplant,dt_min,dt_start,...
    BOUNDARY_CONDITIONS_MATRIX,esa,hs,phsa,ncbot,dx_inter,ncs,arel,brel,tprec);

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
% Javaux M., Lambot S. & Vanclooster M. (2008)
%modified by M.Sall 25/11/09


%%%%%%%%%%%%REMARQUE flxar est ici égale à la variable pot_surface_flux
%%%%%%%%%%%%MAMADOU


%definitions of variables
bctop_changed=var_changed(1);
bcbot_changed=var_changed(2);
dt_changed_bc=var_changed(3);
next_dt_new_bc =var_changed(4);
first_time_bc=var_changed(5);

%initialization
if first_time_bc
    dt_changed_bc=0;
    first_time_bc=0;
    next_dt_new_bc=0;
end
info_pond=0;

% Searches for the appropriate boundary condition for the current time t 

i=max(find(BOUNDARY_CONDITIONS_MATRIX(:,1)<=t));

%if dt_changed_bc==1
%    dt_changed_bc=0;  
    
%    next_dt_new_bc=1; 
%end    

%%%%%%THIS SHOULD BE IMPROVEMENT, BUT IT DOES NOT WORK
%if dt_changed_bc==1 | next_dt_new_bc==1 % In case the dt have been changed to fit the limits of the BC (see calcdt)
%    if dt_changed_bc==1
%        dt_changed_bc=0;  
%        next_dt_new_bc=1;
%    else %(dt_changed_bc==0)     
%        next_dt_new_bc=0;
%        if(BOUNDARY_CONDITIONS_MATRIX(i,2)==2 & boco_top_type==1 & boco_top>0)
            % New user BC is now flux (2) and  different from previous (1)
            % and there is water which has not been infiltrated (boco_top>0);
            % In this case, ponding will appear. 
            % In all other cases, no recalculations are needed
%            imposed_flux=BOUNDARY_CONDITIONS_MATRIX(i,3);
%            info_pond=1;
            %next BC will be ponding
%            BOUNDARY_CONDITIONS_MATRIX(i,2)=1;
%            BOUNDARY_CONDITIONS_MATRIX(i,3)=boco_top; % previous reserve is still considered     
%        end
%    end    
%end 
%%%%%%%%%%%%%%%

if dt_changed_bc==1 % In case the dt have been changed to fit the limits of the BC (see calcdt)
    dt_changed_bc=0;    
    if next_dt_new_bc==1
        next_dt_new_bc=0;
        if(BOUNDARY_CONDITIONS_MATRIX(i,2)==2 & boco_top_type==1 & boco_top>0)
        % New user BC is now flux (2) and  different from previous (1)
        % and there is water which has not been infiltrated (boco_top>0);
        % In this case, ponding will appear. 
        % In all other cases, no recalculations are needed
            imposed_flux=BOUNDARY_CONDITIONS_MATRIX(i,3);
            info_pond=1;
        %next BC will be ponding
            BOUNDARY_CONDITIONS_MATRIX(i,2)=1;
            BOUNDARY_CONDITIONS_MATRIX(i,3)=boco_top; %verifié % previous reserve is still considered     
        end
    end
    next_dt_new_bc=1;
end    



if bctop_changed==1
% Current BC are different of user BC
% The calculated flux by solving the flow equation (calc_boco, check_bc) 
% is incompatible with imposed BC
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
      	BOUNDARY_CONDITIONS_MATRIX(i,3)= pot_surface_flux; %%à revoir la valeur de pot_surface flux
    end
end

%Next Boundary condition for time t
boco_top_type=BOUNDARY_CONDITIONS_MATRIX(i,2) ;
boco_top=BOUNDARY_CONDITIONS_MATRIX(i,3);
boco_bot_type=BOUNDARY_CONDITIONS_MATRIX(i,4);
boco_bot=BOUNDARY_CONDITIONS_MATRIX(i,5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% In case of the presence of a ground water table
A= BOUNDARY_CONDITIONS_MATRIX(i,4);
B= BOUNDARY_CONDITIONS_MATRIX(i,5);
if A== 5 | A== 6 |A== 7;
 [ncbot,boco_bot,gwl,BOUNDARY_CONDITIONS_MATRIX]=...
    find_gwl(A,B,ph,dx_inter,ncs,BOUNDARY_CONDITIONS_MATRIX,phsurf,arel,brel);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Upper boundary condition
if boco_top_type==1
   %succion/pression
   flxa1=esa;
   phsurf=boco_top;
elseif boco_top_type==2										
   % flux, pluie
   %flxa1=esa+boco_top; %%changé par Mamadou
   i=max(find(BOUNDARY_CONDITIONS_MATRIX(:,1)<=t));
   prec=tprec(i,2);
   flxa1=esa+prec; 
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


