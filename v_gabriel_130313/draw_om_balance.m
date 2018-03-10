%function draw_om_balance(om_balance1, tnode1,units, N_reaction_balance1, om_balance2, tnode2, N_reaction_balance2)
%function draw_om_balance(om_balance1, tnode1,units, N_reaction_balance)

dcorg1 = om_balance(:,8); 
dnorg1 = om_balance(:,9);
cberr_corg1 = om_balance(:,10); 
cberr_norg1 = om_balance(:,11);
tflcorg1 = om_balance(:,12);
tflnorg1 = om_balance(:,13);

clear ptup ptmin pthyd ptnit ptvol ptden
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

totalMin = (-ptmin(1:end,1)-ptmin(1:end,2) -ptmin(1:end,3))*10^5;
totalUp = (-ptup(1:end,1)-ptup(1:end,2) -ptup(1:end,3))*10^5;
totalfirstorder = pthyd + ptnit + pthyd + ptden;

% %Nitrogen balance
figure
plot(tnode,totalMin(2:end), tnode,totalUp(2:end),tnode, totalfirstorder(2:end))
legend('total mineralisation', 'total uptake','firstorder reactions')
%  figure
%  plot( tnode1, tflnorg1-tflnorg1(147), tnode2, tflnorg2)
%  title ('input')
%  figure
%  plot(tnode1, dnorg1-dnorg1(147), tnode2, dnorg2)
%  title('delta N')
%  figure
%  plot(tnode1, cberr_norg1, tnode2, cberr_norg2)
% title ('error')
 
% legend('loss via miner/gain via immoob', 'input', 'delta N', 'error')
 xlabel ('time')
%  ylabel([units{1,3} '/' units{1,1} '^2'])
%  title ('Nitrogen Balance')



% %Carbon balance
% subplot(2,1,1)
% plot(tnode, -om_balance(:,7), tnode, tflcorg, tnode, dcorg, tnode, cberr_corg)
% legend('loss via CO_2', 'input', 'delta carbon', 'error')
% xlabel ('time')
% ylabel([units{1,3} '/' units{1,1} '^2'])
% title('Carbon Balance')
% 
% %Nitrogen balance
% subplot(2,1,2)
% plot(tnode,-ptmin(2:end,1)- ptmin(2:end,2) -ptmin(2:end,3), tnode, tflnorg, tnode, dnorg, tnode, cberr_norg)
% legend('loss via miner/gain via immoob', 'input', 'delta N', 'error')
% xlabel ('time')
% ylabel([units{1,3} '/' units{1,1} '^2'])
% title ('Nitrogen Balance')
% 
% %Detailed description of balance
% %Organic carbon balance
% figure
% subplot(2,1,1)
% plot(tnode, om_balance(:,1:2:5), tnode, -om_balance(:,7),tnode, dcorg, tnode, cberr_corg)
% legend('litter', 'manure', 'humus','co2', 'change in profile', 'error' )
% title('carbon balance in profile (change in profile is due to CO_2 production or input from top)')
% ylabel([units{1,3} '/' units{1,1} '^2'])
% xlabel(['time (' units{1,2} ')'])
% % 
% % %Organic nitrogen balance
% subplot(2,1,2)
plot(tnode, om_balance(:,2:2:6),tnode,-ptmin(2:end,1)- ptmin(2:end,2) -ptmin(2:end,3), tnode, dnorg, tnode, cberr_norg)
legend('litter', 'manure', 'humus', 'N due to mineral/immm' , 'change in profile', 'error')
title('organic nitrogen balance in profile (loss in profile is due to mineralisation, gain is from immobilisation or input from top)')
ylabel([units{1,3} '/' units{1,1} '^2'])
xlabel(['time (' units{1,2} ')'])



