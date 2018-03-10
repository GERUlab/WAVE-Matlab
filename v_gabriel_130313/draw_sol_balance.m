function draw_sol_balance(solute_balance,units,nsol,N_reaction_balance,sim_nitro)

%Draw the solute balance

 %solute_balance{1,1}(:,2:end) = solute_balance{1,1}(:,2:end)*100000;
 %solute_balance{2,1}(:,2:end) = solute_balance{2,1}(:,2:end)*100000;
 %solute_balance{3,1}(:,2:end) = solute_balance{3,1}(:,2:end)*100000;

for i=1:nsol  
t = solute_balance{i,1}(:,1);
tflsol=solute_balance{i,1}(:,2);
rleasa = solute_balance{i,1}(:,3);
dsol=solute_balance{i,1}(:,4);
cberr_sol = solute_balance{i,1}(:,5);
tcsink = solute_balance{i,1}(:,6);

if sim_nitro
figure
subplot(2,1,1)
plot(t, tflsol, t, rleasa, t,dsol,t, cberr_sol, t,tcsink  )
legend('cum solute inflow', 'cum netto flux at bottom', 'change of solute in profile','error on the balance', 'sink term')
xlabel('time')
ylabel([units{1,3} '/' units{1,1} '^2'])
title (['Solute Balance of solute ' num2str(i) ])

ptup(:,1)=N_reaction_balance(:,1);
ptup(:,2)=N_reaction_balance(:,2);
ptup(:,3)=N_reaction_balance(:,3);
ptmin(:,1)=N_reaction_balance(:,4);
ptmin(:,2)=N_reaction_balance(:,5);
ptmin(:,3)=N_reaction_balance(:,6);
pthyd=N_reaction_balance(:,7);
ptnit=N_reaction_balance(:,8);
ptvol=N_reaction_balance(:,9);
ptden=N_reaction_balance(:,10);
subplot(2,1,2)
if i == 1 ;%ureum
plot(t, -pthyd)
legend('hydrolyse ureum')
elseif i==2;%ammonium
plot(t, -ptup(:,2), t, ptmin(:,2),t,-ptnit,t,-ptvol,t,pthyd)%,t, -ptup(:,2)+ptmin(:,2)-ptnit-ptvol+pthyd)
legend('uptake ammonium', 'mineralisation ammonium', 'nitrification ammonium','volatilisation','+ hydrolyse ureum');%, 'sink term')
elseif i ==3 ;%nitrate
    plot(t, -ptup(:,3), t,ptmin(:,3), t, ptnit,t,-ptden)%,t, -ptup(:,3)+ptmin(:,3)+ptnit-ptden )
    legend('uptake nitrate', 'mineralisation nitrate', 'nitrification ', 'denitrification')%,'sink term')
end
xlabel('time')
ylabel([units{1,3} '/' units{1,1} '^2'])
if i==1
title ('Origin of sink term of urea')
elseif i==2
    title ('Origin of sink term of ammnium')
elseif i==3
   title ('Origin of sink term of nitrate')
end

% subplot(3,1,3)
% if i == 1 ;%ureum
% plot(t, -pthyd, t,tcsink)
% legend('sink term nit', 'sink term sol')
% elseif i==2;%ammonium
% plot(t, -ptup(:,2)+ptmin(:,2)-ptnit-ptvol+pthyd,t,tcsink,t, -ptup(:,2)+ptmin(:,2)-ptnit-ptvol+pthyd-tcsink)
% legend('sink term nit', 'sink term sol', 'error')
% elseif i ==3 ;%nitrate
%     plot(t, -ptup(:,3)+ptmin(:,3)+ptnit-ptden,t,tcsink)
%     legend('sink term nit', 'sink term sol')
% end
% xlabel('time')
% ylabel([units{1,3} '/' units{1,1} '^2'])
else
  

figure

plot(t, tflsol, t, rleasa, t,dsol,t, cberr_sol, t,tcsink  )
legend('cum solute inflow', 'cum netto flux at bottom', 'change of solute in profile','error on the balance', 'sink term')
xlabel('time')
ylabel([units{1,3} '/' units{1,1} '^2'])
title (['Sink term of solute ' num2str(i) ])
end

end