function check_ini(soil_parameters,dt_start,dt_min,dt_max);
%check parameters set
%check initial values for dt
%Javaux, M.
%%%%%%%%%%%%%%%%%%%%%%%%%%

if dt_start<dt_min
    error('dt_start < dt_min');
elseif dt_start>dt_max
    error('dt_start > dt_max');
end

[validity] = lambda_control(soil_parameters);
if validity==0
      error('Jeu de parametres non valide');
end