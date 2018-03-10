function [conc,solute_applic]=In_solute_boundary_sol(t,nsol,tprec,wat_flxs,sol_upper_bc_type);
%solute_boundary:find the actual flux concentration
%compartiment 
%IN:
%t: time information (see general data) 
% solute_applic_matrix
%nsol: number of solutes
%OUT
% conc: solute concentration (ML-3)
% CALL:none
%CALLED BY:solve_solute
%--------------------------------
% M. Sall 14/04/09
     
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% SOLUTE APPLICATION
%%% solute_applic1=Concentration in mg/cm³ (MG CM-3)=
%% concentration of infiltrating water for solute upper boundary condition type=Concentrationtion Flux boundary Condition 
% or concentration at the upper boundary for solute upper boundary condition type=Concentration Boundary Condition
%%column1= time, column2 = solute 1,  column3 = solute 2, column4 = solute3
           
%%% solute_applic2=Concentration in mg/cm² (MG CM-2)
% possible only if solute upper boundary condition type = Concentration Flux Condition 
%%as in WAVE FORTRAN else set concentration to zéro 




%% SOLUTE APPLICATION
%%% solute_applic1=Concentration in mg/cm³ (MG CM-3)
solute_applic1=...
[0	0	0	163
50 0 0 0]; 

%%% solute_applic2=Concentration in mg/cm² (MG CM-2)
if sol_upper_bc_type==1;  %%%'concentration_flux_bc'
solute_applic2=...
[0	0	0	0
50 0 0 0];
end


%%%%%%%%%%Calculation of the actual concentration
i=max(find(solute_applic1(:,1)<=t));
if numel(i)~=0
    concfo=(solute_applic1(i,2:end));
else
    concfo=zeros(1,nsol);
end

if sol_upper_bc_type==1;   %%%'concentration_flux_bc';
    
    j=max(find(solute_applic2(:,1)<=t));
    if numel(i)~=0
        concw=(solute_applic2(i,2:end));
    else
        concw=zeros(1,nsol);
    end
    u=max(find(tprec(:,1)<=t));
    prec=tprec(u,2);
    if prec ~= 0
        concf=concfo*abs(prec); %%% Dans ce cas conc est équivalent à solsur pour WAVEFOR
    else
        concf=zeros(1,nsol);
        concw=zeros(1,nsol);
    end  
conc =concf+concw;

elseif sol_upper_bc_type==2;   %%%'concentration_bc' 
    conc=concfo;
else
    error('solute upper boundary condition type is undefined')
end



solute_applic=[];



