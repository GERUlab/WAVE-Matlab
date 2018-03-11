function write_profile_data(scenario_number)
%Write inputfile for the profile data
ncs = 85; %Check with the previous inputfile!
nsol = 3

%Van Genuchten parameters for each soil compartment
%Calibrated
% soil_parameters(1:30,:) = ones(30,1)*[0.01 0.41 0.0047 1.4 50 0.5 1];
% soil_parameters(31:85,:) = ones(55,1)*[0.01 0.41 0.0047 1.15 50 0.5 1];

%test
soil_parameters(1:30,:) = ones(30,1)*[0.01 0.41 0.003 1.4 50 0.5 1];
soil_parameters(31:85,:) = ones(55,1)*[0.01 0.41 0.004 1.2 50 0.5 1];

%Measured
% soil_parameters(1:30,:) = ones(30,1)*[0.01 0.41 0.008 1.4 50 0.5 1];
% soil_parameters(31:85,:) = ones(55,1)*[0.01 0.41 0.0047 1.15 50 0.5 1];


%Solute parameters for each soil layer
%bd lm ratio alf f
bd = 1.5;
lm = 8.5;
ratio = 1;
alf = 0;
f = 1;

% Ureum
%Rates Kd diff aa bb
rates = 0;
kd = 0;
diff = 0.1^2; %Values according to Ducheyne (2001)
aa = 0.1;
bb = 10;
solute{1,1} = ones(ncs,1)*[rates kd diff aa bb];
% Ammonium
rates = 0;
kd = 0.0025;
diff =10* (0.1^2); %Values according to Ducheyne (2001)*
aa = 0.1;
bb = 10;
solute{2,1} = ones(ncs,1)*[rates kd diff aa bb];
% Nitrate
rates = 0;
kd = 0;
diff = 10*(0.1^2); %Values according to Ducheyne (2001)
aa = 0.1;
bb = 10;
solute{3,1} = ones(ncs,1)*[rates kd diff aa bb];

%Decaying constants (can be specified per compartment)
knitrif= 0.1;                   % nitrification constant   (T-1)
kdenit= 0.04 ;                     % denitrification constant (T-1)
khydro= 0.000;                     % hydrolyse constant      (T-1)
kvolat= 0.00  ;                    % volatilisation constant  (T-1)
k_lito = 00.05;%/24/60;     %Litter 
k_mano = 000.05;%/24/60;    %Manure
k_humo = 00.002*1.25;%/24/60;    %Humus 
om_param = ones(ncs,1)*[knitrif, kdenit, khydro, kvolat, k_lito, k_mano, k_humo];


solute_param1 = ones(ncs,1)*[bd lm ratio alf f];
for j=1:nsol
    if j==1
    solute_param2(:,1) = solute{j,1}(:,1);
    solute_param2(:,nsol+1) = solute{j,1}(:,2);
    solute_param2(:,2*nsol+1) = solute{j,1}(:,3);
    solute_param2(:,3*nsol+1) = solute{j,1}(:,4);
    solute_param2(:,4*nsol+1) = solute{j,1}(:,5);
    else
    solute_param2(:,j) = solute{j,1}(:,1);
    solute_param2(:,j+nsol) = solute{j,1}(:,2);
    solute_param2(:,j+2*nsol) = solute{j,1}(:,3);
    solute_param2(:,j+3*nsol) = solute{j,1}(:,4);
    solute_param2(:,j+4*nsol) = solute{j,1}(:,5);    
end
end

outputmatrix = [soil_parameters(:)' solute_param1(:)' solute_param2(:)' om_param(:)']

xlswrite('outputmatrix_profile.xls', outputmatrix)