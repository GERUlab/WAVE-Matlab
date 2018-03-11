function [cs,tcsolo_ini] = In_initial_solute(nsol,ncs);

%Initial solute present in profile
%IN: 
%   nsol:number of solute spieces
%   ncs: number of compartiments
%OUT:
%   tcsolo_ini: initial mass of solute in the profile (M/L2)
%   cs: Concentration at soil surface (M/L3)

% CALL:none¨
% CALLED BY: WAVE_MAT main programme
%---------------------------------------
% M. Sall  26/11/08

% MATRIX WITH INITIAL DATA SOLUTE FOR EACH SOLUTE AND EACH COMPARTIMENT 
% ONE COLUMN=ONE SOLUTE
%if nitrogen is modelled:
% solute 1 = ureum
% solute 2 = ammonium
% solute 3 = nitrate
 sol_data=ones(ncs,nsol)*0.0;
 sol_data(1:30,1) = 23/100000/30;
 sol_data(1:30,2) = 23/100000/30;
 sol_data(1:30,3) = 23/100000/30;
 
  sol_data(31:60,1) = 5/100000/30;
 sol_data(31:60,2) = 5/100000/30;
 sol_data(31:60,3) = 5/100000/30;
%
 sol_data(61:ncs,1) = 5/100000/25;
 sol_data(61:ncs,2) = 5/100000/25;
sol_data(61:ncs,3) = 5/100000/25;

sol_data(:,1) = 0 ;



%Calculations with input
tcsolo_ini=sol_data;
reservoir=0*ones(nsol,1); % reservoir initial for each solute(on soil surface)     
cs=ones.*(reservoir);



