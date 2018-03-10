function [BOUNDARY_CONDITIONS_MATRIX,hs,phsa,arel,brel,initial_gwl,tprec]=In_Boundary_conditions;

% Enter your top and bottom boundary conditions in a matrix called:
%   BOUNDARY_CONDITIONS_MATRIX
% This matrix is composed by:
% Col 1: time 
% Col 2: type of top boundary condition
%		1 = pressure head condition 
%		2 = flow condition : + for evaporation, - for drainage and/or irrigation
% Col 3: top boundary condition
% Col 4: type of bottom boundary condition
%		1 = pressure head condition
%		2 = flow condition : + for capillarity rise, - for drainage 
%		3 = seepage face (lysimeter)
%		4 = free drainage
%%%AJOUT MAMADOU
%		5 = a groundwater is present and the groundwater level (gwl) is given at each time
%		6 = a groundwater is present and the fonction between gwl and flux is available
%		7 = a groundwater is present and the flux is given at each time

% Col 5: bottom boundary condition
% 		if Col 4 = 3, then Col 5 must be equal to 0
%       if Col 4 = 5, then Col 5 = gwl  in absolute value
%       if Col 4 = 6, then Col 5 = initial gwl 
%       if Col 4 = 7, then Col 5 = flux 
%%%%%%%%%%%%%%%%%%%%%%
% IN :
%	none
% OUT :
% 	BOUNDARY_CONDITIONS_MATRIX
%	hs:maximum ponding depth (cm)
%	phsa:Minimum allowed pressure head at the surface (cm)
%CALLS:
%	none
%
%%%%%%%%%%%%%%%%%%%%%%%
%User boundary condition

BOUNDARY_CONDITIONS_MATRIX = [00,2,-1,4,0
                              20,2,0.,4,0
                              40,2,-1,4,0
                              50,2,0.,4,0];


tprec=[BOUNDARY_CONDITIONS_MATRIX(:,1) BOUNDARY_CONDITIONS_MATRIX(:,3)];


hs=0;								% Maximum ponding depth (cm)
phsa=-1E4; %origine				% Minimum allowed pressure head at the surface (cm)

%Specify the iniitial gwl
initial_gwl=95;
% for a groundwater flux relationship as bottom condition give the
% coefficient used in the relation Q = arel*exp(brel*gwl)
arel= -0.01;
brel= -0.005;

