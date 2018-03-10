function [dt,dt_changed_bc]=calc_dt(t,dt,iter,dt_min,dt_max,BOUNDARY_CONDITIONS_MATRIX,nsol,...
    tprec,wat_flxs,simsol,compartiments_number,rtex,dx,wc_change_tol,sol_upper_bc_type,dt_start)

% CALC_DT   Calculates the optimal time step
%
%
% IN : 
%	dt : the actual time step 
%	iter : the number of "solve_flow" iterations
%   dt_min,dt_max
% OUT: 
% 	dt = the optimized time step (min)
% CALL:
%	none
% CALLED BY:
%   wavemat101.m
%-----------------------
% M. Vanclooster 2/2/2000
%modified by MSall Mars 2010





if (iter<=5)
   dt=dt*2;
end
if ((iter>5)&(iter<=10))
   dt=dt*1.6;
end
if ((iter>10)&(iter<=20))
   dt=dt*1.3;
end
if ((iter>20)&(iter<=25))
   dt=dt;
end
if (iter>25)
   dt=dt*0.7;
end

if dt<dt_min
   dt=dt_min;
end
if dt>dt_max
   dt=dt_max;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MAM
% define maximum timestep in terms of iteration in WAVE.FOR
%if((maxiter.gt.1).and.(ndts.gt.0))then
%fac=0.5+0.25*(7-iter);
%fac= max(0.5,fac);
%dt=dt*fac;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MAM


% Adapt dt in terms of boundary conditions
i=max(find(BOUNDARY_CONDITIONS_MATRIX(:,1)<=t));        % at time t
i2=max(find(BOUNDARY_CONDITIONS_MATRIX(:,1)<=(t+dt)));  % at time t+dt
if i~=i2            % dt encompass two or more changes of BC
    if i~=i2-1      % dt encompass moe changes of BC
        i2=1+i;
	end
    dt=BOUNDARY_CONDITIONS_MATRIX(i2,1)-t; % reduction of the time step, this time step can be smaller than dtmin
    dt_changed_bc=1;
    if dt < dt_min  %if the new dt<dtmin, do the change in this iteration
        i=i2;
        dt_changed_bc=0;
        dt=dt_start;
    end
else
    dt_changed_bc=0;
end

%%%AJOUT MSALL 08 04 10
%%     define maximum time step in terms of water content change
dtdt=1.0e-10;
flxsn=-(rtex(1)+(wat_flxs(1)-wat_flxs(2))/dx);
dtdt= max(dtdt, abs(flxsn));
for i=2:compartiments_number
    flxsn=-(rtex(i)+(wat_flxs(i)-wat_flxs(i+1))/dx);
	dtdt= max(dtdt, abs(flxsn));
end

%%     define maximum time step in terms of solute transport
if simsol
   [conc,solute_applic]=In_solute_boundary(t,nsol,tprec,wat_flxs,sol_upper_bc_type);
   i=max(find(tprec(:,1)<=t));
   dt_sol = dt;
    sp = 1: nsol;
   if conc(sp) > 0
        dt_sol= min(dt_sol,0.0099);
        if tprec(i+1,1)-tprec(i,1)<0.0099
            dt_sol =dt_min;
        end
    end
    
else
    dt_sol=100;
end
dt_sol=min(dt,dt_sol);
dt= max(dt_min, min(wc_change_tol/dtdt,dt_sol));