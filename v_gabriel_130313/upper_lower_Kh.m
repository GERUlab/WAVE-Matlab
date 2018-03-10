function [K_pf42, K_pf_upper]=upper_lower_Kh(g1,g2,g3,wc_wp, part,wcs,K_i,wc_fc,K_vapour_in, K_in)

%Calculates upper and lower limit of the De Vries calculation of thermal
%conductivity.

%For moisture contents between wilting point (pf 4.2) and (wilting point +
%field capacity)/2 thermal conductivity is interpolated based on the output
%of this m-file

%In:
%g1,g2,g3: g-factors for the solid materials
%wc_wp : water content at wilting point
%part: fractions of the solid materials
%wcs: saturated water content
% K_i : thermal conductivity of the solid materials
% wc_fc : water content at field capacity
% K_vapour_in : intial value for the conductivity of the vapour
% K_in : initial value for the conductivity

%Out
% K_pf42 : thermal conductivity at pF 4.2 (wilting point)
% K_pf_upper: thermal conductivity at (wilting point +
%field capacity)/2

%Called by:
% calc_thermal_conduc4

%% Lower boundary

%Initialize g-factors for the solid materials
g1_pf42 = g1;
g2_pf42 = g2;
g3_pf42= g3;

%Set water content to water content at wilting point
part(:,5) = wc_wp;
part(:,4) = 1-part(:,1)-part(:,2)-part(:,3)-part(:,5);

%Correct conductivity of air with vapour
K_i(:,4) = K_in + K_vapour_in;

%Calculate g-factors for air
g1_pf42(:,4) = 0.333 - (0.333-0.070)*part(:,4)./(wcs - wc_wp);
g2_pf42(:,4) = g1_pf42(:,4);
g3_pf42(:,4) = 1-2*g1_pf42(:,4);

%Calculate f
h=[1 2 3 4 5];
f_pf42 =...
1/3*(...
1./(1+((K_i(:,h)./repmat(K_i(:,4),1,5))-1).*g1_pf42(:,h))+ ...
1./(1+((K_i(:,h)./repmat(K_i(:,4),1,5))-1).*g2_pf42(:,h)) + ...
1./(1+((K_i(:,h)./repmat(K_i(:,4),1,5))-1).*g3_pf42(:,h)));

%Calculate conductivity of the soil (medium is air, thus leaving the fourth
%column out the sum)
h = [1 2 3 5];
K_pf42 = (part(:,4) .* K_i(:,4) + sum(f_pf42(:,h).*part(:,h).*K_i(:,h),2))./(part(:,4) + sum(f_pf42(:,h).*part(:,h),2));

%% Upper Boundary
%Initialize g-factors for the solid materials
g1_pf_upper = g1;
g2_pf_upper = g2;
g3_pf_upper= g3;

%Set water content to water content at (wilting point + field capacity)/2
wc_upper = (wc_wp + wc_fc)/2;
part(:,5) = wc_upper;
part(:,4) = 1-part(:,1)-part(:,2)-part(:,3)-part(:,5);

%Calculate g-factors for air
g1_pf_upper(:,4) = 0.333 - (0.333-0.070)*part(:,4)./(wcs - wc_wp);
g2_pf_upper(:,4) = g1_pf_upper(:,4);
g3_pf_upper(:,4) = 1-2*g1_pf_upper(:,4);

%Calculate f
h=[1 2 3 4 5];
f_pf_upper =...
1/3*(...
1./(1+((K_i(:,h)./repmat(K_i(:,5),1,5))-1).*g1_pf_upper(:,h))+ ...
1./(1+((K_i(:,h)./repmat(K_i(:,5),1,5))-1).*g2_pf_upper(:,h)) + ...
1./(1+((K_i(:,h)./repmat(K_i(:,5),1,5))-1).*g3_pf_upper(:,h)));

%Calculate conductivity of the soil (medium is water, thus leaving the
%fifth column out the sum)
h = [1 2 3 4];
K_pf_upper = (part(:,5) .* K_i(:,5) + sum(f_pf_upper(:,h).*part(:,h).*K_i(:,h),2))./(part(:,5) + sum(f_pf_upper(:,h).*part(:,h),2));