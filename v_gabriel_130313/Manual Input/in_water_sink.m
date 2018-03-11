function [p0, p1, p2h, p2l, p3, smx, ilinr] = in_water_sink(ncomp,dx);

%Inputfile for water sink term variables

% IN:
% ncomp: number of compartments
% dx = space increment
% 
% OUT:
% p0   :below this pressure head start extraction (above = anaerobis)
% p1   :below this pressure head optimal extraction
% p2h  :below this pressure head no longer optimal at high evaporative demand
% p2l  :below this pressure head no longer optimal at low evaporative demand
% smx: maximal root extraction rate (1/T)
% ilinr: Boolean for linear relationship alpha - pressure head

%%%%%%%%%%%%%%
%Vansteenkiste J., 25/02/2010

%Linear relationschip between reduction factor of the root water uptake
%(alpha) and the pressure head? (yes = 1, no = 0=hyperbolic);
ilinr = 1;

%Specify critical pressure heads
p0 =-15;    %Below this pressure head start extraction (above = anaerobis)
p1 = -30 ; %below this pressure head optimal extraction
p2h = -600; %below this pressure head no longer optimal at high evaporative demand
p2l = -600 ;%below this pressure head no longer optimal at low evaporative demand
p3 =-8000 ;%below this pressure head no extraction (= wilting point)

%Linear relationship between potential root water uptake and depth
%(yes = 1, no = 1)
ismxnl = 1;

%next lines are not used, because these parameters are also input in
%In_ETsplit.
if ismxnl== 1 
    %Linear, specify arer and brer of the relationship:
    %root uptake = arer + brer * x where x is the depth of the node
        arer = 0.5/24/60;         
        brer = -0.0/24/60;
else
    %No linear relationship, specify maximal root water uptake for each
    %compartment
end

%%%%%%%%%%%%%%%%%END OF INPUT SECTION%%%%%%%%%%%%%%%%%%%%%%%

% Calculates maximal root water uptake and check if it is not negative
i=1:ncomp;
if ismxnl
rt_upt = arer + brer*abs(dx*i);
dummy = rt_upt>0;
index = find(dummy ==0);
if isempty(index)==0
error('root water uptake is negative, check in_water_sink')
end

smx = rt_upt;
end

%Check the critical pressure head
if p1 > p0| p2h>p1| p2l >p1| p3>p2h| p3>p2l
    error('check the input values of the critical pressure heads in in_water_sink')
end