%Control charts for crop cauli
Colheaders = {'DATE','DAP','GP','RMAINT','Mean Day Temp',...
    'Te','TeSum','DVS','TotDM','DMRoot',...
    'DMLeaf','LeafArea','DMStem','DMCurd','TotDMLeaf',...
    'TotDMStem','TotDMCurd','TotDMRoot','TotDmShoot','TotLeafArea',...
    'LAICov','LAI','TotDMGreenLeaf','TotDMDeadLeaf','TotRootDepth1',...
    'RCil','Transpiration','Evaporation','CumTranspiration','CumEvaporation',...
    'LeafStemNDemand','CurdNDemand','RootNDemand','PlantNDemand','TotNDemandLeafStem',...
    'TotNDemandCurd','TotNDemandRoot','TotNDemandPlant'};

ColheadersNup ={'CumPlantNUptakeActualCurd','CumPlantNUptakeActualLeaf','CumPlantNUptakeActualStem',...
    ' CumPlantNUptakeActualLeafStem','CumPlantNUptakeActualRoot','PlantNUptakeActualCurd',...
    'PlantNUptakeActualLeaf','PlantNUptakeActualStem',' PlantNUptakeActualLeafStem',...
    'PlantNUptakeActualRoot','NFac','PercNShoot'};

cauli = output_crop{1,1};
cauliNup = output_cropNup{1,1};
DAP = 1:1:size(cauliNup,1);

figure
k = [1 2 3 5]
plot(DAP, cauliNup(:,k))
legend(ColheadersNup{k})
xlabel('DAP')
ylabel('Nuptake (gr/plant)')

% figure
% k = [11];plot(DAP, cauliNup(:,k));xlabel('DAP');ylabel(ColheadersNup{k}); ylim([0 1])
% figure
% k = [12];plot(DAP, cauliNup(:,k));xlabel('DAP');ylabel(ColheadersNup{k});


% figure
% subplot(2,2,1) ;%DVS
% k =[7 8];plot(cauli(:,k(1)), cauli(:,k(2)));xlabel(Colheaders{k(1)});ylabel(Colheaders{k(2)})
% subplot(2,2,2);%TotLeafArea
% k =[2 20];plot(cauli(:,k(1)), cauli(:,k(2)));xlabel(Colheaders{k(1)});ylabel(Colheaders{k(2)})
% subplot(2,2,3)%LAI
% plot(cauli(:,2), cauli(:,[21 22]))
% xlabel(Colheaders{1});ylabel('LAI')
% legend(Colheaders{21},Colheaders{22})
% subplot(2,2,4);%Rooting depth
% k =[2 25];plot(cauli(:,k(1)), cauli(:,k(2)));xlabel(Colheaders{k(1)});ylabel(Colheaders{k(2)})

%Dry matter
figure
plot(cauli(:,2), cauli(:,[16 17 18 19 23 24 ]));
xlabel(Colheaders{k(2)}); ylabel('gram/plant')
legend(Colheaders{16} ,Colheaders{17},Colheaders{18},Colheaders{19} ,Colheaders{23},Colheaders{24})

% %NDemand
% figure
% plot(cauli(:,2), cauli(:,34:38));
% legend(Colheaders{34:38})

