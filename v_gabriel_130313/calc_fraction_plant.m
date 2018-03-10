function fraction_plant = calc_fraction_plant(crop_type,ncrop,CropParams, crop_result)

 if crop_type(ncrop)==1    
               rootDAP = [0 20 40 68 120];
                x_y_param = [2.88 2.88 9.92 13.02 13.02];
                threshold = 0.05;
               rootdensparamxy =interp1(rootDAP,x_y_param,crop_result(2));distance = [0:0.01:40];cumulative_dist = normcdf(distance,0,rootdensparamxy);
               index = min(find(cumulative_dist > 1-threshold));
               fraction_plant = min(((distance(index)/100)^2*pi)/(1/CropParams(5))*1.2,1);

               
 else crop_type(ncrop)==2    
               rootDAP = [0 32 68 102 200];
                x_y_param = [10 31 59 68 68]; %due to overlap, root coverage is here used as a direct input
                threshold = 0.05;
               leek_rootcoverage =interp1(rootDAP,x_y_param,crop_result(2));
               leek_rootcoverage = max(leek_rootcoverage*1.4,1);
               fraction_plant = leek_rootcoverage/100;
 end         
