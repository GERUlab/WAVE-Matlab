function [dt_sol,applic_boolean,dt_sol_count] = calc_dt_sol(t,nsol,dt, applic_boolean,dt_min, units,cs,csp,dt_sol_count,dx,reservoir,fsol, conirsol,wdepsol,ddepsol,In_solute_boundary,theta_table)

%Adapt dt in terms of solute boundary conditions
% IN: 
% t : current time
% nsol : number of solutes
% dt: timestep
% applic_boolean: boolean indicating end of a solute application (1 = true, 0 = false)
% dt_min: minimum allowed timestep
% units: units of Mass, Time and Space
%OUT:
% dt_sol = timestep for the solute module
% applic_boolean: boolean indicating end of a solute application (1 = true, 0 = false)    

%Joachim Vansteenkiste, 11/03/2010

%Find the current application
%[conc,solute_applic,solboco_top_type,solboco_input_type]=In_solute_boundary(t,nsol);
solute_applic = In_solute_boundary.solute_applic;
solboco_top_type = In_solute_boundary.solboco_top_type;
solboco_input_type = In_solute_boundary.solboco_input_type;
i=max(find(solute_applic(:,1)<=t));
if numel(i)~=0
conc=(solute_applic(i,2:end));
else
    conc=zeros(1,nsol);
end
    

if solboco_input_type == 2 %input is defined in M L-2
%find current climatic conditions
climate = theta_table;
index = max(find(t>=climate(:,1)));
current_rain = climate(index, 3);

%Find current applications
index = max(find(t>=fsol(:,1)));
current_fert = fsol(index,:);
index = max(find(t>=conirsol(:,1)));
current_irr = conirsol(index,:);
%Wet and dry deposition
ressol=current_rain*dt.*wdepsol + ddepsol;
%Inorganic fertiliser + Fertigation
ressol = ressol + current_fert(2:end) + current_irr(2:end);
cs = ressol + reservoir;
else
cs = conc;
end


%If there is application decrease the timestep
if sum(cs) ~=0
    if strcmp(units{1,2},'day' )==1
    dt1 = 0.099;
    elseif strcmp(units{1,2},'hour' )==1
    dt1 = 0.0099*24;
    elseif strcmp(units{1,2},'min' )==1 
    dt1 = 0.0099*24*60;
    elseif strcmp(units{1,2},'s' )==1
    dt1 = 0.0099*24*60*60;
    end
else
    dt1 = dt;
end

if solboco_input_type == 1
%Adapt dt to the top BC
i=max(find(solute_applic(:,1)<=t));        % at time t
i2=max(find(solute_applic(:,1)<=(t+dt)));  % at time t+dt
if i~=i2            % dt encompass two or more changes of BC
    if i~=i2-1      % dt encompass moe changes of BC
        i2=1+i;
	end
    dt2=solute_applic(i2,1)-t; % reduction of the time step, this time step can be smaller than dtmin
else
    dt2 = dt;
end
else
    dt2 =dt;
end

previous_time_step = sum(csp);
next_time_step = sum(cs);
if applic_boolean ==0 && isnan(previous_time_step)==0;;
if previous_time_step ~=0 && next_time_step ==0 %application is finished, smaller time steps are required
    applic_boolean = 1;
    dt_sol_count = 0;
    dt3 = dt_min;
else
    dt_sol_count = 0;
    dt3 = dt;
end
else dt3 =dt;
end

if applic_boolean ==1
    dt_sol_count = dt_sol_count +1;
    dt3 = dt_min;
    if dt_sol_count == 10
        applic_boolean = 0;
    end
end
dt_sol = min(dt1, dt2);
dt_sol = min(dt_sol, dt3);

