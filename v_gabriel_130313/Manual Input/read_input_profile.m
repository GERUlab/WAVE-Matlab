function [soil_parameters, solute_param1, solute_param2, om_param] =  read_input_profile(ncomp, scenario_file, scenario_number)

%aantal lagen uit vorige dataset
% ncomp = 85;
% 
% scenario_file= ('input_profiel.xlsx');
scenario_number = 1;
%Read datafile
input_dataset = xlsread(scenario_file)

%Select the correct line to be used as an input and reshape
soil_parameters = reshape (input_dataset(scenario_number,1:ncomp*7),ncomp, 7);
solute_param1 = reshape (input_dataset(scenario_number,ncomp*7+1:ncomp*12),ncomp, 5);
solute_param2 = reshape (input_dataset(scenario_number,ncomp*12+1:ncomp*27),ncomp, 15);
om_param = reshape(input_dataset(scenario_number,ncomp*27+1:ncomp*34),ncomp, 7);