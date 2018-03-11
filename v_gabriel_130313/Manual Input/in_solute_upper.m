function [ddepsol, wdepsol, fsol, conirsol,ppdepth]= in_solute_upper(nsol,dx,scenario_number)

%Solute Upper Boundary Condition as in Wave_Fortran
%When this boundary condition is needed, specify this in the file
%In_solute_boundary

% In:
% nsol : number of solutes (-)
% dx: space discretisation (L)
% OUT:
% ddepsol: average deposition of species (M L-2 T-1)
% wdepsol average load of species in precipitation (M L-3)
% fsol: Inorganic fertilization (M L-2)
% conirsol : irrigation (Fertigation) (M L-2)
% ppdepth: plowing depth (L)

%J. Vansteenkiste 11/03/2010
dose = 4 

%Deposition
%Wet deposition (WDEPSOL): average load of species in precipation (M L-3)
%Dry deposition (DDEPSOL): average daily deposition of species (M L-2)
ddepsol = zeros(1,nsol);
wdepsol = zeros(1,nsol);

%Inorganic fertilization (M L-2)
%[t ureum ammonium nitrate]


    if dose ==4
%Dose 4
fsol = [0 0 0 0
         68 0/10^5 90/10^5 90/10^5 
         114 0/10^5 9/10^5 9/10^5
         163 0 0 0 
         212 0/10^5 70/10^5 70/10^5
         366 0 0 0
         441 0/10^5 87/10^5 87/10^5
         484 0/10^5 70.5/10^5 70.5/10^5
         534 0 0 0
         585 0/10^5 71.5/10^5 71.5/10^5
         731 0 0 0
         798 0 82/10^5 82/10^5
         841 0 66.5/10^5 66.5/10^5
         874 0 0 0
         954 0 37.5/10^5 37.5/10^5];
    
    elseif dose ==3
%Dose 3     
fsol = [0 0 0 0
         68 0/10^5 65/10^5 65/10^5 
         114 0/10^5 12.5/10^5 12.5/10^5
         163 0/10^5 0 0 
         212 0/10^5 50/10^5 50/10^5
         366 0 0 0
         441 0/10^5 62/10^5 62/10^5
         484 0/10^5 54/10^5 54/10^5
         534 0 0 0
         585 0/10^5 53/10^5 53/10^5
         731 0 0 0
         798 0 52.5/10^5  52.5/10^5
         841 0 37.5/10^5 37.5/10^5
         874 0 0 0
         954 0 48/10^5 48/10^5];
     
elseif dose ==2     
%Dose 2
fsol = [0 0 0 0
         68 0/10^5 40/10^5 40/10^5 
         114 0/10^5 17.5/10^5 17.5/10^5
         163 0 0 0
         212 0/10^5 40/10^5 40/10^5
         366 0 0 0
         441 0/10^5 37/10^5 37/10^5
         484 0/10^5 32/10^5 32/10^5
         534 0 0 0
         585 0/10^5 37.5/10^5 37.5/10^5
         731 0 0 0
         798 0 25.5/10^5 25.5/10^5
         841 0 26.5/10^5 26.5/10^5
         874 0 0 0 
         954 0 25/10^5 25/10^5         ];
     
elseif dose ==1     
%Dose 1
fsol = [0 0 0 0
         68 0/10^5 15/10^5 15/10^5 
         114 0/10^5 21.5/10^5 21.5/10^5
         163 0 0 0
         212 0/10^5 20/10^5 20/10^5
         366 0 0 0
         441 0/10^5 12/10^5 12/10^5
         484 0/10^5 17/10^5 17/10^5
         534 0 75/10^5 0
         585 0/10^5 8.5/10^5 8.5/10^5
         874 0 0 0];
    end 
fsol(:,1) = fsol(:,1)+733774;
%Irrigation (Fertigation) (M L-2)
%[t urem ammonium nitrate]
conirsol = [0 0 0 0
            65 0 0 0];

%Plowing (L)
%[t plowing depth]
ppdepth = [0 0
            10 0];

%Calculate the compartments affected by plowing
ppdepth(:,3) = ceil(ppdepth(:,2)/dx);