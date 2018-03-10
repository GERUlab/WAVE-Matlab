function soil_parameters=In_soil_parameter_data(ncs);

% Give for each compartiment the different parameters of Van Genuchten
% equation parameters:
% 			=1 if no hysteresis
% 			>1 if hysteresis
% IN:
% 	ncs: nber of compartiments
% OUT:
%	soil_parameters
% CALLS:
%	none
% 	Column 1 Thr
% 	Column 2 Ths
% 	Column 3 alfa (1/cm)
% 	Column 4 n
% 	Column 5 Ks (cm/min)
% 	Column 6 lambda
% 	Column 7 alfa_w/alfa_d 



soil_parameters=ones(ncs,1)*[0.08 0.43 0.04 1.6 50 0.5 1];

%Initial parameter, silt loam 2, Lambot, WRR03