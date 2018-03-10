
dx_inter(1)=0.5; dx_inter(2:100)=1; dx_inter(101)=0.5;
A=5;ncs=100;
if A==5  
gwl= 95;
dxn=cumsum(dx_inter);
i=max(find (dxn <gwl ));


ncomp=i;
    if ncomp == ncs;
   error('ground water level reached bottom of profile')
    end
ph(ncomp+1)=-abs(gwl)+abs(dxn(ncomp+1)) %%determine ph(ncomp+1) assuming unit gradient
break
ncbot=ncomp;
boco_bot =0;   %%boco_bot=ph(ncomp+1);              %% Dirichlet condition

  %bocobot=0;ncbot=94;
end