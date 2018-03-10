function[temp,rocp,lamda,ttop,tbot]=solve_temperature(t,dt,dx,WC,wcs,tempp,ttopp,...
    tbotp,rocpp,lamdap,ncs);

%MODULE TEMPERATURE
%SOLVE_TEMPERATURE solve the 1D heat flow equation
%IN
%t=time(mn)
%dt=time increment(mn)
%dx=compartiment size
%WC=water content(-); wcs=saturated water content
%tempp=temperatute at time t-dt
%ttopp= soil surface temperatute at time t-dt
%TBOT= soil bottom temperatute at time t-dt
%ncs=compartiment number
%rocpp,lamdap:  parameters at time t-dt
%OUT
%temp=temperature at time t
%ttop= soil surface temperatute at time t
%TBOT= soil bottom temperatute at time t
%rocp,lamda:  parameters at time t

%CALL: In_temperature_param

%CALLED BY: WAVE_MAT main programme
%-------------------------------------------------------------------------
% Mamadou SALL 04/12/08

% lamdap= lamda at  t-1
% LAMDAPM= lamda at t-1/2
%rocpp=rocp  t-1
%ROCPPM=rocp  t-1/2
[ttop,tbot,rocp,lamda]=In_temperature_param(t,WC,wcs,...
   ncs,ttopp,tbotp);


for I=[1:ncs];  
lamdapm(I)=(lamdap(I)+lamda(I))/2;
rocppm(I)=(rocp(I)+rocpp(I))/2;
end

%%I=1;              %%%%%%%%(I-1) replaced by I

AA(1,1)=dt*(lamdapm(1)+lamdapm(1))/2/(2*dx*(2*((rocppm(1)+rocppm(1))/2)*dx));

BA(1,1)= dt*(lamdapm(2)+lamdapm(1))/2/(2*dx*(2*((rocppm(2)+rocppm(1))/2)*dx));

A(1,1)=BA(1,1);                                         %%%%%%%%%%%%%%%%% tempj(i+1)
B(1,1)=1-AA(1,1)+BA(1,1);                                    %%%%%%%%%%%%%%%%% tempj(i)                                 
E(1,1)=tempp(1)+AA(1,1).*(ttop+tempp(1))-BA(1,1).*(tempp(2)+tempp(1))...
    +AA(1,1).*ttop;     %%%%%%%%%% Cte

 for I=2:ncs-1;
AA(I,1)=dt*(lamdapm(I-1)+lamdapm(I))/2/(2*dx*(2*((rocppm(I-1)+rocppm(I))/2)*dx));

BA(I,1)= dt*(lamdapm(I+1)+lamdapm(I))/2/(2*dx*(2*((rocppm(I+1)+rocppm(I))/2)*dx));

A(I,1)=BA(I,1);                                         %%%%%%%%%%%%%%%%% Tj(i+1)
B(I,1)=1-AA(I,1)+BA(I,1);                                    %%%%%%%%%%%%%%%%% Tj(i)
D(I-1,1)=-AA(I,1);                                        %%%%%%%%%%%%%%%%% Tj(i-1)
E(I,1)=tempp(I)+AA(I,1).*(tempp(I-1)+tempp(I))-BA(I,1).*(tempp(I+1)+tempp(I));     %%%%%%%%%% Cte
 end
%%I=ncs;                    %%%%%%%%(I+1) remplcé par I
AA(ncs,1)=dt*(lamdapm(ncs-1)+lamdapm(ncs))/2/(2*dx*(2*((rocppm(ncs-1)+rocppm(ncs))/2)*dx));

BA(ncs,1)= dt*(lamdapm(ncs)+lamdapm(ncs))/2/(2*dx*(2*((rocppm(ncs)+rocppm(ncs))/2)*dx));

B(ncs,1)=1-AA(ncs,1)+BA(ncs,1);                                    %%%%%%%%%%%%%%%%% Tj(i)
D(ncs-1,1)=-AA(ncs,1);                                        %%%%%%%%%%%%%%%%% Tj(i-1)
E(ncs,1)=tbot+AA(ncs,1).*(tempp(ncs-1)+tbot)-BA(ncs,1).*(tbot+tbot)-...
    BA(ncs,1).*tbot;         

%SOLUTION OF THE TRIDIAGONAL MATRIX
temp=(inv(diag(B)+diag(A,1)+diag(D,-1))*E).';
temp=ones(1,ncs)*20;

