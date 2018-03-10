%Process Holding results
clear
load crop_model Holding
Holding = Holding(2:end,:)
%[RRDATE,RRDAP,i,RDMLa(i), RLengthLa(i),RVolumeLa(i),ERLa(i),RLDLa(i),RLDRCil];

% 1: data
% 2: dap
% 3 : aantal lagen
% 4 : droge stof
% 5 : lengte
% 6: volume
% 7 : erla
% 8 : rldl
% 9: rldrcil

for JD = 1:100
index = find(Holding(:,2) == JD);
numbers = size(index,1);
LLA(1:numbers,JD) = Holding(index,5);
LLA(numbers:85,JD)= zeros(85-numbers+1,1)
end

for DAP =1:100
 rootDAP = [0 20 40 68 120];
x_y_param = [2.88 2.88 9.92 13.02 13.02];
threshold = 0.05;
rootdensparamxy =interp1(rootDAP,x_y_param,DAP);
distance = [0:0.01:40];cumulative_dist = normcdf(distance,0,rootdensparamxy);
index = min(find(cumulative_dist > 1-threshold));
radius(DAP) = distance(index);
end

area = radius.^2*pi; %cm²

RLD = LLA(:,end)/area(end);