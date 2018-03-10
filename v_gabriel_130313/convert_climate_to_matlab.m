function crop_climate = convert_climate_to_matlab(crop_climate)

%Convert the julian day in the inputfile to matlabtime
for i=1:size(crop_climate,1)
matlabtime(i,1) = jday2matlab(crop_climate(i,2), crop_climate(i,1));
end

crop_climate = [matlabtime crop_climate(:,3:7)];