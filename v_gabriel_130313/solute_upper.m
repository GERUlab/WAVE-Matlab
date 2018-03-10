function [fsol, conirsol, solsur] = solute_upper(t, fsol, conirsol, wdepsol, ddepsol,dt,fraction_time,theta_table);

global strategy
%Process the upper boundary conditions
% IN: 
% solute_param1= [bd,lm,ratio,alf,f]
%where:
%bd=bulk density,
%lm= hydrodynamic dispersivity, 
%ratio= mobilewater content/total water content)
%alf=mass coefficient transfer betwween mobile and immobile zones of
%       the soil, f=sorption sites fraction in the mobile zone of  soil
%solute_param2= [rates,kd,dif,ar,br];
%where:
%rates= solute decaying coefficients, 
%kd= distribution coefficient,
%dif =chemical difffusion in pure water
%ar and br cefficients for the calculation of dispersion coefficient
% t: current time
% fsol : inorganic fertilisation (M L-2)
% conirsol : irrigation (M L-2)
% ppdepth : plowing depth (L)
% wdepsol: wet deposition (M L-3)
% ddepsol: dry deposition (M L-2)
% tcsolo: total concentration
% csolmo
% csolio
% acsolmo
% acsolio
% wcio
% wco
% nsol : number of solute (-)
% immobile : boolean for immobile concept
% dx : space discretisation (L)
% dt: time step (T)
% wcmo: 

%find current climatic conditions
climate = theta_table;
index = max(find(t>=climate(:,1)));
current_rain = climate(index, 3);


%Find current applications and set to zero (one time application)
index = max(find(t>=fsol(:,1)));
current_fert = fsol(index,:);
if sum(fsol(index,2:end))==0
    current_fert(2:end) = fsol(index,2:end);
else
    current_fert(2:end) = fsol(index,2:end);
%     if strcmp(strategy, 'Row')
%     current_fert(2:end) = fsol(index,2:end)/fraction_time(end);
%     elseif strcmp(strategy, 'Broad')
%     current_fert(2:end) = fsol(index,2:end);
%     elseif strcmp(strategy, 'Fert')
%     current_fert(2:end) = fsol(index,2:end)/fraction_time(end);    
%     end
        
end
fsol(index,2:end) = 0;
index = max(find(t>=conirsol(:,1)));
current_irr = conirsol(index,:);
conirsol(index,2:end)=0;

%Wet and dry deposition
ressol=current_rain*dt.*wdepsol + ddepsol*dt;
%Inorganic fertiliser + Fertigation
ressol = ressol + current_fert(2:end) + current_irr(2:end);


% C       DEPLETION DURING LAST DAY
%         RESSOL(SP)=DMAX1(RESSOL(SP),MN_REAL)-
%      $  DMAX1(SOLINFL(SP),MN_REAL)


solsur = ressol;


% C     CALCULATE SOLSUR
%       DO 100 SP = 1,NR_OF_SOL
%         IF(PREC(NDAY).GT.0.E0.OR.RIRR(NDAY).GT.0.E0)THEN
%           SOLSUR(SP)  = RESSOL(SP)
%         ELSE
%           SOLSUR(SP) = 0.E0
%         ENDIF
% 100   CONTINUE
%       RETURN
%       END	 