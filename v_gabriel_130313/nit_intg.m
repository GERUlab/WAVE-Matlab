function[ptsnorg, ptscorg ptup, cupt, pthyd, chyd, ptnit, cnit, ptvol, cvol, ptden, cden,...
 dnorg, dcorg, cberr_norg, cberr_corg,tflnorg,tflcorg,ptmin,cmin,om_balance] = ...
 nit_intg(ptsnorg,ptscorg, ptup, cupt, ptmin, cmin, pthyd, chyd, ptnit, cnit, ptvol,...
 cvol, ptden, cden, soil_om,tcorgs,tnorgs,carbman,carblit,rnitman,rnitlit,decnorg,deccorg,...
 uptakem,uptakei,minerm,mineri,nitrifm,nitrifi,denitm,deniti,hydro_uream,hydro_ureai,...
    volatm,volati,tflnorg,tflcorg,dx,dt,om_balance)

%IN
% ptsnorg : cumulative organic nitrogen sink
%ptscorg :  cumulative organic carbon sink
%ptup : profile total solute uptake
%cupt : cumulative plant uptake in each compartiment
%ptmin :  profile total mineralisation
%cmin : cumulative mineralisation of nitrogen in each compartiment
%pthyd : profile total hydrolysis
%chyd : cumulative hydrolysis in each compartiment
%ptnit : profile total nitrification
%cnit : cumulative nitrification in each compartiment
%ptvol : profile total volatilisation
%cvol : cumulative volatilisation in each compartiment
%ptden : profile total denitrification
%cden : cumulative denitrification in each compartiment
%tnorgs : initial organic nitrogen profile total
%tcorgs : initial organic carbon profile total
%soil_om : actual organic carbon and notrogen in the profile
%tflnorg
%tflcorg
%OUT
%cberr_norg : mass balance error of the organic nitrogen
%cberr_corg : mass balance error of the organic carbon
%terms of the mass balance 


c_manp = soil_om(:,1);
n_manp = soil_om(:,2);
c_litp = soil_om(:,3);
n_litp = soil_om(:,4);
c_hump = soil_om(:,5);
n_hump = soil_om(:,6);
co2 = soil_om(:,7);
%initialisation
ptclit = 0; ptnlit = 0; ptcman = 0; ptnman = 0; ptchum = 0; ptnhum = 0; 
volume=dx;
%profile total
ptclit = sum(c_litp);  %%actual profile total
ptnlit = sum(n_litp);  %%actual profile total
ptcman = sum(c_manp);  %%actual profile total
ptnman = sum(n_manp);  %%actual profile total
ptchum = sum(c_hump);  %%actual profile total
ptnhum = sum(n_hump);  %%actual profile total
ptco2 = sum(co2);
ptsnorg = ptsnorg + sum(decnorg)*volume*dt;  %% cumulative profile total sink (ptsnorg conservé)
ptscorg = ptscorg + sum(deccorg)*volume*dt;  %% cumulative profile total sink (ptscorg conservé)

ptnitorg=ptnlit+ptnman+ptnhum;      %% stock actuel d'azote
ptcarborg=ptclit+ptcman+ptchum;      %%stock actial de carbone   

%profile sink

%ammoniuum uptake
ptup(2) = ptup(2)+ (sum(uptakem(:,2))+ sum(uptakei(:,2)))*dx*dt;  % profile total uptake
cupt(:,2) = cupt(:,2)+ (uptakem(:,2) + uptakei(:,2))*dx*dt;       % cumulative uptake
%ammonium mineralisation
ptmin(2) = ptmin(2) + (sum(minerm(:,2))+ sum(mineri(:,2)))*dx*dt; % profile total mineralisation
cmin(:,2) = cmin(:,2) + (minerm(:,2)+mineri(:,2))*dx*dt;           % cumulative mineralisation
%nitrate uptake
ptup(3)= ptup(3)+ (sum(uptakem(:,3))+ sum(uptakei(:,3)))*dx*dt;  % profile total uptake
cupt(:,3) = cupt(:,3)+ (uptakem(:,3) + uptakei(:,3))*dx*dt;       % cumulative uptake
% nitrate mineralisation
ptmin(3) = ptmin(3) + (sum(minerm(:,3))+ sum(mineri(:,3)))*dx*dt; % profile total mineralisation
cmin(:,3) = cmin(:,3) + (minerm(:,3)+mineri(:,3))*dx*dt;           % cumulative mineralisation

% hydrolise
pthyd = pthyd + (sum(hydro_uream) + sum(hydro_ureai))*dx*dt;  % profilr total
chyd = chyd + (hydro_uream + hydro_ureai)*dx*dt;

% nitrification
ptnit = ptnit + (sum(nitrifm) + sum(nitrifi))*dx*dt;
cnit = cnit + (nitrifm + nitrifi)*dx*dt;
% volatilisation
ptvol = ptvol + (sum(volatm) + sum(volatm))*dx*dt;
cvol = cvol + (volatm + volati)*dx*dt;
% denitrification
ptden = ptden + (sum(denitm) + sum(deniti))*dx*dt;
cden = cden + (denitm + deniti)*dx*dt;

% organic balance term
dnorg=ptnitorg - tnorgs; %%(variation stock = actual-initial)
dcorg=ptcarborg - tcorgs; %% (variation stock = actual-initial)

% organic mass balance error
tflnorg = tflnorg+rnitlit+rnitman; %Input (defined in in_om_applic)
tflcorg = tflcorg+carblit+carbman; %Input (defined in in_om_applic)

cberr_norg= dnorg - (tflnorg - ptsnorg); %% (error = variation stock -(entrées-sorties))
cberr_corg= dcorg - (tflcorg - ptscorg);

om_balance(end+1,:) = [ptclit ptnlit ptcman ptnman ptchum ptnhum ptco2 dcorg dnorg cberr_corg cberr_norg tflcorg tflnorg];
%N_reaction_balance(end+1,:) = [ptup(1) ptup(2) ptup(3) ptmin(1) ptmin(2) ptmin(3) pthyd ptnit ptvol ptden];
