function [cs,tcsolo_ini] = In_initial_solute(nsol,ncs);

%IN INITIAL_SOLUTE
%IN: 
% nsol:number of solute spieces
% ncs: number of compartiments
%OUT:
%sol_data: matrix of mass solute for each soil compartiment
% tcsolo_ini: initial mass of solute in the profile (mg/m²)
% reservoir: solute concentration on the surface soil 
%Initial solute in the soil (mg/m2) for each compartiment

% CALL:none¨
% CALLED BY: WAVE_MAT main programme
%---------------------------------------
% M. Sall  26/11/08

%%%%MATRIX WITH INITIAL DATA SOLUTE FOR EACH SOLUTE AND EACH COMPARTIMENT 
  %%%% ONE COLUMN=ONE SOLUTE
    %%%%%%%% FILL IN
    %%%%      SOL1    SOL2   SOL3 
    
  
         
             %% OR
    sol_data=[[0 0 0.]; ones(ncs-1,nsol)*0.0];
       
% ATTRIBUTION of initial CONCENTRATIONS TO THE COMPARTIMENTS
for j=[1:nsol];
    for i=1:ncs;
           tcsolo_ini(i,j)=sol_data(i,j);
     end
end
        reservoir=0*ones(nsol,1); % reservoir initial for each solute(on soil surface)
       
        cs=ones.*(reservoir);



