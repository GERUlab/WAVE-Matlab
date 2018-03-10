clear
clc

ncomp = 100;
units = {'cm','day','g','°C';};
WC = ones(ncomp,1)*0.1;
wcs = linspace(0.45,0.45,100)'; %ones(ncomp,1)*0.4;
[part, alamo, geomq, geomos,geomom,bd_water,heatcapsolids_corr, heatcapwater_corr]= In_temperature_param(units,ncomp);
geomq(:,1) = linspace(0,0.5,101);
geomq(:,2) = linspace(0,0.5,101);
geomq(:,3) = 1 - geomq(:,1) - geomq(:,2);
dummy = WC> 0.2;
alama = dummy*alamo(1,5) +(1-dummy)*alamo(1,4);
alama(ncomp+1) = alama(ncomp);


for m=1:3;
    if m==1
        g = geomq;
    elseif m==2
        g=geomom;
    elseif m==3
        g=geomos;
    end
    sums1 = 1./(1+((alamo(:,m)./alama) -1).*g(:,1));
    sums2 = 1./(1+((alamo(:,m)./alama) -1).*g(:,2));
    sums3=  1./(1+((alamo(:,m)./alama) -1).*g(:,3));        
    aka(:,m)=(sums1+sums2+sums3)/3;
end

% 


%Calculate lambda_v
dummy = WC> 0.2;
alamv = dummy*alamo(5) +(1-dummy).*(WC*alamo(5)/0.2);
     %decrease from 0.176 mCal cm-1 sec-1 °C to 0 when WC drops from 0.2 to 0
     %in wave fortran different: WC(i)*0.2/alamo(5)
%Fortran: alamv = dummy*alamo(5) +(1-dummy).*(WC*0.2/alamo(5));
alamv(ncomp+1) = alamv(ncomp);

geoma(:,1) = linspace(0,0.5,101);
geoma(:,2) = linspace(0,0.5,101);
geoma(:,3) = 1 - geoma(:,1) - geoma(:,2);
%geoma = calc_temp_g_air(WC, wcs);
%geoma(ncomp+1,:) = geoma(ncomp,:);


%Calculate k-coefficient for air  
    sums = 1./(1+(((alamo(:,4)+alamv)./alama) -1).*geoma(:,1))+...
            1./(1+(((alamo(:,4)+alamv)./alama) -1).*geoma(:,2))+...
            1./(1+(((alamo(:,4)+alamv)./alama) -1).*geoma(:,3));        
    aka(:,4)=1/3*sums;
    
plot(geomq(:,1), aka(:,1), geoma(:,1), aka(:,4))
legend ('quarts', 'air')

 