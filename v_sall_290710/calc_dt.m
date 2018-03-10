function dt = calc_dt(t,dt,iter,dt_min,dt_max,BOUNDARY_CONDITIONS_MATRIX)

% CALC_DT   Calculates the optimal time step
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

% Adapt dt in terms of boundary conditions
i=max(find(BOUNDARY_CONDITIONS_MATRIX(:,1)<=t));        % at time t
i2=max(find(BOUNDARY_CONDITIONS_MATRIX(:,1)<=(t+dt)));  % at time t+dt

if i~=i2            % dt encompass two or more changes of BC
    if i~=i2-1      % dt encompass moe changes of BC
        i2=1+i;
	end
    dt=BOUNDARY_CONDITIONS_MATRIX(i2,1)-t; % reduction of the time step
    dt_changed=1;
%if the new dt is < than dtmin, do the change in this iteration
    if dt < dt_min
        i=i2;
        dt_changed=0;
        dt=dt_start;
    end
end
