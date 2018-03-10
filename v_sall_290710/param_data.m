function [PARAM,U,V,DIS1,Qp,NM,smax,nf,stop,iinit,local,gamma] = param_data()
%PARAM_DATA User defines fixed and fitted parameters and also fitted parameters space in
%	wich the optimization process will be carried out (FITTING_GLOB call) or for which the
%	objective function will be calculated (FITTING call).

%	S. Lambot (October 2000)


% Define fixed parameters and put nan for fitted parameters
p1 = 0.025;			%WCR
p2 = 0.40; 		%WCS [Sol100: 0.337] [HYPRES: 0.403, 0.439, 0.430, 0.520, 0.614] [S: 0.3611, O: 0.4017]
p3 = nan;			%ALFA
p4 = nan;			%N
p5 = 4.2;			%Ks
p6 = 0.5;			%Lambda
p7 = 1  ;			%ALFA_r


PARAM = [p1 p2 p3 p4 p5 p6 p7];

% Define the parameters space in which the optimization will be carried out
U = [0.005, 1.20]';	%U: column vector of MVG lower bounds (Art1)
V = [0.050, 2.00]';	%V: column vector of MVG upper bounds (Art1)

%U = [0.005, 1.05, 1e-2, -3]';	%U: column vector of MVG lower bounds (Sol100)
%V = [0.050, 7.00, 1e+2,  3]';	%V: column vector of MVG upper bounds (Sol100)

%U = [3.57, 0.1362, 1e-1, -3]';	%U: column vector of ASS lower bounds
%V = [269, 1.0243, 1e+2, +3]';	%V: column vector of ASS upper bounds

% MCS arguments definition
np = length(U);
smax = 5*np+10;	% (default Art1: 5*np+10)
	% smax: number of levels (default: 5*n+10).
	% smax governs the relative amount of global versus local search. By
	% increasing smax, more weight is given to global search.

nf = 150*np^2;	% (default Art1: 150*np^2)
	% nf: maximum number of function evaluations (default: 50*n^2).
	% Increasing nf allows to spend more time on boxes that have a chance 
	% on their level to contain better points. This may be important for
	% hard problems, or for problems with very wide bounds.

%stop = 3*np;
stop(1) = 0;
stop(2) = -inf;

% stop: integer defining stopping test (default: stop = 3*n).
	% Increasing or decreasing stop increases or decreases the amount 
	% of work allowed before concluding that nothing more is gained; 
	% the default choice is quite conservative, and may try to locate
	% a better point long after the global minimum has been found.
	% stop = 5 works for many easier problems, too, with much fewer
   % wasted function values.
   
   % stop(1) = 0: the user can specify a function value that should be reached.
	% stop(2) = function value that is to be achieved.


iinit = 2;	%(default Art1: 2)
	% iinit: parameter defining the initialization list.
	% iinit = 0: corners and midpoint.
	% iinit = 2: 5u/6 + v/6, u/6 + 5v/6 and midpoint.
   
local = 75;	%(default Art1: 75)
   % local: maximal number of steps in local search (default: 50).

gamma = 1e-10;	%(default Art1: 1e-10)
	% gamma: stopping criterion for local search (default: eps).
	% Acceptable relative accuracy for local search.
	% A tiny gamma (default) gives a quite accurate but in higher 
	% dimensions slow local search. Increasing gamma leads to less work 
	% per local search but a less accurately localized minimizer.
   
%FITTING-------------------------------------------------------------------------------------------
% Define the discretization of parameters space (for FITTING only).
DIS1 = [50, 50, 50]';
	% DIS2: number of values between lower and upper bounds.
   % The discretization is linear, except for parameter Ks for which it is logarithmic.
   % length(DIS1) = length(U).
   
% define values of fitted parameters to use as constants (for FITTING only).
Qp = [0.1426, 1.5000, 25.9691]';
	% length(Qp) = length(U).
	% These constants serve as splitting hyper plans of the 8D objective function.

% define the different parameters couple for which the objective function will be
% calculated and plotted (for FITTING only).
NM = [4 5; 4 6; 5 6];
	% NM: i*2 matrix containing the indices in PARAM of corresponding parameters.
   % For instance, NM = [3 5; 4 5] will lead to the calculation and plotting
   % of the objective function respectively as a function of alpha and Ks (p3 and p4),
   % and as a function of n and Ks (p4 and p5).