%Control charts for crop cauli
%Labels of the leek output
Colheaders2 = {'DATE','DAP','GP','RMAINT','Mean Day Temp','Te','TeSum','TotDM','DMRoot','DMLeaf',...
    'LeafArea','DMStem','TotDMLeaf','TotDMStem','TotDMRoot','TotDmShoot','TotLeafArea','LAICov','LAI','NLeaf',...
    'NStem','NperLeaf','NperStem','TotNperLeaf','TotNperStem','TotRootDepth1','RCil','Transpiration','Evaporation','CumTranspiration',...
    'CumEvaporation','TotNperShoot','ShaftDi','ShaftLe'};

ColheadersNup2 ={'CumPlantNUptakeActualLeaf','CumPlantNUptakeActualStem','CumPlantNUptakeActualLeafStem',...
    'CumPlantNUptakeActualRoot','PlantNUptakeActualStem',...
    'PlantNUptakeActualLeaf', 'PlantNUptakeActualLeafStem','PlantNUptakeActualRoot','NFac','PercNShoot',...
    'PercNLeaftab', 'PercNStemtab','PlantNDemand'};

leek = output_crop{1,2};
leekNup = output_cropNup{1,2};
DAP = 1:1:size(leekNup,1);

figure
k = [1 2 4];
plot(DAP, leekNup(:,k))
legend(ColheadersNup2{k})
xlabel('DAP')
ylabel('Nuptake (gr/plant)')

figure
subplot(3,1,1) ;%NFac
k = [9];plot(DAP, leekNup(:,k));xlabel('DAP');ylabel(ColheadersNup2{k}); ylim([0 1])
xlim([0 140])
subplot(3,1,2)
k = [11 12]; %Nitrogen percentages
CropParamsfile = importdata(inputfile_leek,',',1);
LeekParams = CropParamsfile.data(1,:);
DAPN = LeekParams(50:57);
PercNShaft = LeekParams(76:83);
PercNLeaf = LeekParams(84:91);
TSumN = LeekParams(92:99);
hold on
plot(DAP, leekNup(:,k),DAPN, PercNShaft, DAPN, PercNLeaf)
xlabel('DAP');ylabel('perc');
xlim([0 140])
legend(ColheadersNup2{(11)},ColheadersNup2{(12)}, 'interpol shaft', 'interpol leaf')
subplot(3,1,3)%Nitrogen demand and uptake
plot(DAP(2:end), diff(leekNup(:,13)), DAP(2:end), Nuptake_time(176:308,2))
legend('demand', 'uptake')
xlim([0 140])

figure
subplot(2,2,1) ;%
k =[7 33];plot(leek(:,k(1)), leek(:,k(2)));xlabel(Colheaders2{k(1)});ylabel(Colheaders2{k(2)})
subplot(2,2,2);%TotLeafArea
k =[2 17];plot(leek(:,k(1)), leek(:,k(2)));xlabel(Colheaders2{k(1)});ylabel(Colheaders2{k(2)})
subplot(2,2,3)%LAI
plot(leek(:,2), leek(:,[18 19]))
xlabel(Colheaders2{1});ylabel('LAI')
legend(Colheaders2{18},Colheaders2{19})
subplot(2,2,4);%Rooting depth
k =[2 26];plot(leek(:,k(1)), leek(:,k(2)));xlabel(Colheaders2{k(1)});ylabel(Colheaders2{k(2)})

%Dry matter
figure
plot(leek(:,2), leek(:,[13 14 15 16 ]));
xlabel(Colheaders2{2}); ylabel('gram/plant')
legend(Colheaders2{13:16})

% %NDemand
% figure
% plot(leek(:,2), leek(:,34:38));
% legend(Colheaders{34:38})

