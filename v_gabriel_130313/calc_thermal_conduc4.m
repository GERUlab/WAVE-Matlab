function [allam, beta] = calc_thermal_conduc4(units,ncomp,soil_parameters, solute_param1,WC,wcp);

%Calculates thermal conductivity based on pure De Vries equations

%In :
% units: specify which units are used in in_general_data
% ncomp: number of compartments
% soil_parameters(:,2) :saturated water conten
% solute_param1(:,1) : bulk density
% WC: water content of the current time step
% wcp: water content of the previous time step

%OUT:
% allam: thermal conductivity between the nodes, halfway the time step 
%in (J L-1 T-1 C-1)
%beta : heat capacity of the soil (J/C/L^3)

% CALLS:
% In_temperature_param
% calc_temp_k
% calc_temp_g_air
% CALLED BY:
% solve_temperature_JV


%% From other modules:
wcs = soil_parameters(:,2); %saturated water content
bd = solute_param1(:,1); %bulk density
bd(ncomp+1) = bd(ncomp);
WC = WC'; %Water content

%% Define extra water contents for calculation allam (ncomp +1) 
%Joachim V. (Not included in fortran)
WC(ncomp +1)= WC(ncomp);
wcp(ncomp +1) = wcp(ncomp);
wcs(ncomp+1) = wcs(ncomp);


    %%calculate water content halfway the node, halfway the time step
    i = [2:ncomp+1];
    tw(1) = (WC(1) + wcp(1))/2;   
    tw(i) = (WC(i) + wcp(i) + WC(i-1) + wcp(i-1))/4;
    tw = tw';

%User defined input of temperature parameters
[part,K_i,bd_water,heatcapsolids_corr, heatcapwater_corr,shape_par, shape_par2, K_vapour]= In_temperature_param4(units,ncomp);    

part(:,5) = tw;
part(:,4) = 1-sum(part,2);
K_in = K_i(:,4);
K_vapour_in = K_vapour;
%% Calculate volumetric heat capacity (eq 4-2 Wave manual)
i = [1:ncomp+1];
beta(i) = heatcapsolids_corr*bd(i) + heatcapwater_corr*tw(i)*bd_water; %in J/C/L^3 
beta = beta';

ph = -(10^4.2);
wc_wp = moist_ret(ph,soil_parameters,1,0)';
wc_wp(:) = 0.11;
ph = -(10^2);
wc_fc = moist_ret(ph,soil_parameters,1,0)';
wc_fc(:) = 0.27;
wc_wp(ncomp+1) = wc_wp(ncomp);
wc_fc(ncomp+1) = wc_fc(ncomp);

%%Calculate K_air corrected for vapour transfer
dummy4 = tw >= wc_wp;
K_vapour = dummy4.*K_vapour + (1-dummy4).*K_vapour.*tw./wc_wp;
K_i(:,4) = K_i(:,4) + K_vapour;
%% Determine the medium 
%Medium is water between saturation and (pF 2 + pF 4.2)/2
%Medium is air between residual water content and pF 4.2
%Between (pF 2 + pF 4.2)/2 and pF 4.2 no medium is assumed thermal cond is
%interpolated

dummy1 = tw> (wc_wp + wc_fc)/2; %medium water
dummy2 = tw< wc_wp ;%medium air
dummy3 = 1-dummy1-dummy2;
K_0 = dummy1.*K_i(:,5) + dummy2.*K_i(:,4)+dummy3*9999;
x_0 = dummy1.*part(:,5) + dummy2.*part(:,4)+dummy3*9999;


%% Calculate g-factors
index = find(shape_par ==1); %Spherical
alpha = shape_par2(index);
for i=1:length(index)

if alpha(i) <1
    %g1(i) = 1/2 * alpha^2/(1-alpha^2)*(1/(alpha^2) + 1/(2*(1-alpha^2)^1/2)*log(( 1-(1-alpha^2)^(1/2))/(1+(1+alpha^2)^(1/2))
elseif alpha(i) ==1
    g1(index(i),:) = 1/3;
elseif alpha(i) > 1
    g1(:,index(i)) = 1/2 * alpha(i)^2/(alpha(i)^2-1)*(pi/(2*(alpha(i)^2-1)^(1/2)) - 1/(alpha(i)^2) - 1/((alpha(i)^2-1)^(1/2))*atan((alpha(i)^2-1)^(-1/2)));
    g2(:,index(i)) = g1(:,index(i));
    g3(:,index(i)) = 1-2*g1(:,index(i));
end
end

clear index
index = find(shape_par ==0); %Cylinder
m = shape_par2(index);
for i=1:length(index)
g1(:,index(i)) =  1/(m(i)+1);
g2(:,index(i)) = m(i)/(m(i)+1);
g3(:,index(i)) = 0;
end

g1 = repmat(g1,ncomp+1,1);g2 = repmat(g2,ncomp+1,1);g3 = repmat(g3,ncomp+1,1);

%Value for air
g1(:,4) = 0.333 - (0.333-0.070)*part(:,4)./(wcs - wc_wp);
g2(:,4) = g1(:,4);
g3(:,4) = 1-2*g1(:,4);


%% Calculate f
i=[1:5];
f = 1/3*(1./(1+((K_i(:,i)./repmat(K_0,1,5))-1).*g1(:,i)) + 1./(1+((K_i(:,i)./repmat(K_0,1,5))-1).*g2(:,i)) + 1./(1+((K_i(:,i)./repmat(K_0,1,5))-1).*g3(:,i)));

[K_pf42, K_pf_upper]= upper_lower_Kh(g1,g2,g3,wc_wp, part,wcs,K_i,wc_fc,K_vapour_in,K_in);
%% Calculate conductivity
i=[i:4];k=[1 2 3 5];
A =(wc_wp + wc_fc)/2;
X = [wc_wp(1,1) A(1,1)];
Y = [K_pf42(1,1), K_pf_upper(1,1)];
K_h = dummy1.*((x_0 .* K_0 + sum(f(:,i).*part(:,i).*K_i(:,i),2))./(x_0 + sum(f(:,i).*part(:,i),2)))+...
     dummy2.*((x_0 .* K_0 + sum(f(:,k).*part(:,k).*K_i(:,k),2))./(x_0 + sum(f(:,k).*part(:,k),2)))+...
    dummy3.*spline(X,Y,tw);

allam = K_h;
