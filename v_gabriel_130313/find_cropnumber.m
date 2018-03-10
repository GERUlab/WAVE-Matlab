function ncrop = find_cropnumber (t,simplant,plant_date, harvest_date);

boolean_plant = (t>= plant_date);
boolean_harvest = (t<= harvest_date);
boolean_plant_harvest = boolean_plant + boolean_harvest;
index = find(boolean_plant_harvest == 2);

if isempty(index) ~= 1
    ncrop = index;
else
    ncrop = NaN;
end