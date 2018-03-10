function draw_soilprop(soil_parameters, phsa,units)

%Plot the moisture retention curve and hydraulic conductivity of each
%unique soil soil layer

%find the unique soil layers 
layers = unique(soil_parameters, 'rows');
ph_range = phsa:0.1:0;
dt =1;
update_historic = 1;


for i=1:size(layers,1)
soil_parameters = ones(length(ph_range),1)*(layers(i,:)); 

%Calculate the complete moisture retention curve
wc = moist_ret(ph_range,soil_parameters,dt,update_historic);
kh = conduct(ph_range,soil_parameters);

figure;
semilogx(abs(ph_range), wc)
xlabel(['|pressure head| ' units{1,1}])
ylabel (['water content ' units{1,1} '^3 ' '/' units{1,1} '^3'])
title (['Soil Retention Curve of layer ' num2str(i)])

figure;
semilogx(abs(ph_range), kh)
xlabel(['|pressure head|' units{1,1}])
ylabel (['Hydraulic conductivity ' units{1,1} '/' units{1,2}])
title (['Hydraulic conductivity of layer  ' num2str(i)])

end