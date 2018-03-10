function solute_balance = store_solute_balance(cberr_sol,tflsol,rleasa,tcsink,tsoli,tcsolo_ini,dsol,tcsolo,solute_balance,t,nsol,dt)

%Store the info calculated in sol_intgr in a variable named solute_balance.

% %IN:
% t: time
% tflsol : cumulative solute inflow (M L^-2)
% rleasa: cumulative netto flux at bottom of the profile (M L^-2)
% dsol: change in solute content for the entire profile (M L^-2)
% cberr_sol: error on the balance (M L^-2)
%tcsink : sink term (M L^-2)
% solute balance : summary of all above variables before

%OUT:
% solute balance : updated version of the solute balance

for i=1:nsol
solute_balance{i,1}(end+1,:)= [t+dt tflsol(i) rleasa(i) dsol(i) cberr_sol(i) tcsink(i)];
end

