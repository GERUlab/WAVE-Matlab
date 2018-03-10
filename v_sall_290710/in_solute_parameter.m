function [nsol,solute_param1,solute_param2,ncs,sol_upper_bc_type]=in_solute_parameter(compartiments_number,immobile);

%IN_SOLUTE PARAMETER
%%%IN: 
% compartiments_number
%OUT: 
%nsol:number of solute spieces
%ncs: number of compartiments
%solute_param1= [bd,lm,ratio,alf,f]
%where:
%bd=bulk density, lm= hydrodynamic dispersivity, ratio= (mobile
%water content/total water content), alf=mass coefficient transfer betwween mobile and immobile zones of
%the soil, f=sorption sites fraction in the mobile zone of  soil
%solute_param2= [rates,kd,dif,ar,br];
%where:
%rates= solute decaying coefficients, kd= distribution coefficient, dif=
%chemical difffusion in pure water, ar and br cefficients for the
%calculation of dispersion coefficient

 % CALL: none
 % CALLED BY: WAVE_MAT main programme
%----------------------------------------
%M. Sall  23/03/09


%% SOLUTE BOUNDARY TYPE choose one of the two cases 

%sol_upper_bc_type = 1; %%%%{'concentration_flux_bc'}  %%%(first case)
sol_upper_bc_type = 2; %%%['concentration_bc'] (second case)


ncs=compartiments_number;
nsol=3;                              %%%3 nomber de solutes
nlay=5;                            %%%%%%%% nomber of soil layers
%%% number of compartiment for each soil layer
ncs1=3;
ncs2=3;
ncs3=6;
ncs4=14;


%%%%SOIL LAYER PARAMETERS 
%%%% For each layer
rates=ones(ncs,nsol)*0.00;
lm=[10;10;10;10;10];                        %soil solute dispersivity(L)
bd=[1.52;1.52;1.52;1.52;1.52];                % bulk density (ML^-3)
if immobile
ratio=[0.5;0.5;0.5;0.5;0.5];        %%%mobile/total moisture content ratio(-)
alf=[0.0;0.0;0.0;0.0;0.0];          %%% mass transfer coefficient (T^-1)
f=[1.0;1.0;1.0;1.0;;1.0];                %%%adsorbed site fraction in tne mobile zone (-) 
else
    ratio=ones(ncs,1);alf=zeros(ncs,1);f=ones(ncs,1);
end
%%%%%SOLUTE SPIECES PARAMETRES 
%%%%%%for each spieces, for each soil layer,specify:
KD=[0.0,0.0,0.0; 0.0,0.0,0.0; 0.0,0.0,0.0; 0.0,0.0,0.0; 0.0,0.0,0.0];       %%distribution coefficient(M^-1L^-3)
diff=[0.0,0.0,0.0; 0.0,0.0,0.0; 0.0,0.0,0.0; 0.0,0.0,0.0; 0.0,0.0,0.0];     %%chem dffusion in pure water (L²T^-1)
aa=  [0.01,0.01,0.01; 0.01,0.01,0.01; 0.01,0.01,0.01; 0.01,0.01,0.01 ;0.01,0.01,0.01];   %%empirical coefficient (-)
bb= [10,10,10; 10,10,10; 10,10,10; 10,10,10; 10,10,10];                                  %%empirical coefficient (-)

%%%%%% ATTRIBUTION OF SOIL LAYER PARAMETERS TO THE COMPARTIMENTS  
 i=[1:ncs1];
    lm(i)=lm(1);bd(i)=bd(1);
if immobile
    ratio(i)=ratio(1);alf(i)=alf(1);f(i)=f(1);
end
i=[ncs1+1: ncs1+ncs2];
lm(i)=lm(2);bd(i)=bd(2);
if immobile
    ratio(i)=ratio(2);alf(i)=alf(2);f(i)=f(2);
end
i=[ncs1+ncs2+1: ncs1+ncs2+ncs3];
lm(i)=lm(3);bd(i)=bd(3);
if immobile
    ratio(i)=ratio(3);alf(i)=alf(3);f(i)=f(3);
end
i=[ncs1+ncs2+ncs3+1 :ncs];
lm(i)=lm(4);bd(i)=bd(4);
if immobile
    ratio(i)=ratio(4);alf(i)=alf(4);f(i)=f(4);
end
i=[ncs1+ncs2+ncs3+ncs4+1 :ncs];
lm(i)=lm(5);bd(i)=bd(5);
if immobile
    ratio(i)=ratio(5);alf(i)=alf(5);f(i)=f(5);
end

%% ATTRIBUTION OF SOLUTE SPIECES PARAMETRES FOR EACH COMPARTIMENT
for j=[1:nsol];
    for i=[1:ncs];
        if i<=ncs1
        kd(i,j)=KD(1,j);dif(i,j)=diff(1,j);ar(i,j)=aa(1,j);br(i,j)=bb(1,j);
        end
        if ncs1<i & i<=ncs1+ncs2
            kd(i,j)=KD(2,j);dif(i,j)=diff(2,j);ar(i,j)=aa(2,j);br(i,j)=bb(2,j);
        end
        if ncs1+ncs2<i & i<=ncs1+ncs2+ncs3
           kd(i,j)=KD(3,j);dif(i,j)=diff(3,j);ar(i,j)=aa(3,j);br(i,j)=bb(3,j);
        end
        if ncs1+ncs2+ncs3<i &i<=ncs1+ncs2+ncs3+ncs4;
              kd(i,j)=KD(4,j);dif(i,j)=diff(4,j);ar(i,j)=aa(4,j);br(i,j)=bb(4,j);
        end
        if ncs1+ncs2+ncs3+ncs4<i &i<=ncs
              kd(i,j)=KD(5,j);dif(i,j)=diff(5,j);ar(i,j)=aa(5,j);br(i,j)=bb(5,j);
        end
    end
end
    solute_param1=[bd,lm,ratio,alf,f];
    solute_param2=[rates,kd,dif,ar,br];