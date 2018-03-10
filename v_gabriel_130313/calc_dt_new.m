function [dt,dt_changed_bc]=calc_dt(t,dt,iter,dt_min,dt_max,BOUNDARY_CONDITIONS_MATRIX,rtex,wat_flxs,dx,ncs)

% Calculates the optimal time step
%
%
% IN : 
%	dt : the actual time step (min)
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


wc_change_tol = 0.002;


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

dt_start = dt;
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

% dtdt=1.0e-10;
% flxsn=-(rtex(1)+(wat_flxs(1)-wat_flxs(2))/dx);
% dtdt= max(dtdt, abs(flxsn));
% for i=2:ncs
%     flxsn=-(rtex(i)+(wat_flxs(i)-wat_flxs(i+1))/dx);
% 	dtdt= max(dtdt, abs(flxsn));
% end
% dt_sol = dt_max;
% dt= max(dt_min, min(wc_change_tol/dtdt,dt_sol));
