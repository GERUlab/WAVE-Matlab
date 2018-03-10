%function plot_of()
%PLOT_OF Calculate and plot the i objective functions OFi as a function of two 
%	variable parameters NM(i).
%	Description of the variables :
%		P1...Pk: vectors containing the different values taken by each parameters.
%		DIS2: k vector containing the length of the vectors P1...Pk.
%		Q: k vector containing the parameters values who are kept constant during
%			the different obj. func. calculations (= parameters values used to slice the 
%			hyper space 8D in order to plot 3D graphs (2D contour graphs)).
%		NM: i*2 vector containing the different variable parameters couples.
%		IND: indice vector who allows to select specific values in the parameters
%			vectors P1...Pk.
%		n, m: identifies the two variable parameters.
%		X0: soil parameters vector for which the objective function is calculated.
%		OF1...OFi: 2D matrices containing the objective functions values.

%	INPUTS: function to edit before running FITTING
%		.DATA_INV_SOL 	define observed data.
%		.PARAM_DATA 	define fixed and variable parameters.
%							define parametric space and its discretization for which the
%								objective function will be calculated and plotted.
%							define objective function slicing parameters.
%		.see inputs of function MAIN to define the direct problem.

% S. Lambot (October 2000).


% global variables shared with OBJ_FUN
global DATA_INV weigth_ph weigth_wc weigth_fl
global PARAM
global Ymodel Ymeas ERR WEIGHT
global cpt nf

% read observed data
disp('Loading of measured data ...');
[DATA_INV] = calc_obs_data;
Ymeas = DATA_INV(:,2);
[weigth_ph,weigth_wc,weigth_fl] = calc_weigth(DATA_INV);

% read parameters data
[PARAM,U,V,DIS1,Qp,NM,smax,nf,stop,iinit,local,gamma] = param_data;
l = 1;
for k = 1:length(PARAM)
   if ~isnan(PARAM(k))
      eval(['P' num2str(k) '=' num2str(PARAM(k)) ';']);
      Q(k) = PARAM(k);
   else
      if k == 5
         Pk = logspace(log10(U(l)),log10(V(l)),DIS1(l));
      else
         Pk = linspace(U(l),V(l),DIS1(l));
      end
      eval(['P' num2str(k) '=' 'Pk' ';']);
      Q(k) = Qp(l);
      l = l+1;
   end
   DIS2(k) = length(eval(['P' num2str(k)]));
end
PARAM = ones(1,length(PARAM)).*nan;

% preliminary calculations
warning off;
parameter_text = {'{\theta}_{r} [m³/m³]','{\theta}_{s} [m³/m³]',...
      '{\alpha} [cm^{-1}]','n [{\fontsize{8}-}]','Ks [cm/min])',...
      '\lambda [{\fontsize{8}-}]','{\alpha}_{w}/{\alpha}_{d} [{\fontsize{8}-}]'};

% objective functions calculation
disp('Objective functions calculation...');
for i = 1:size(NM,1)
   IND = ones(1,length(PARAM));
   cpt = 0;
   n = NM(i,1);
   m = NM(i,2);
   nf = DIS2(n)*DIS2(m);
   while cpt < nf
      % calculate soil parameters vector X0
      X0 = Q;
      X0(n) = eval(['P' num2str(n) '(' num2str(IND(n)) ')']);
      X0(m) = eval(['P' num2str(m) '(' num2str(IND(m)) ')']);

      % calculate objective function value at point of indice IND
	   [OF] = obj_fun_plot(X0);
	   	
	   % storing of the objective function
      eval(['OF' num2str(i) '(IND(' num2str(n) '),IND(' num2str(m) '))=OF' ';']);
      
     	% parameters updating
        if IND(n) < DIS2(n)
           IND(n) = IND(n)+1;
        else
           IND(n) = 1;
           if IND(m) < DIS2(m)
              IND(m) = IND(m)+1;
           end
        end
   
	end
   
   % plot i objective functions
   figure;
   	contourf(eval(['P' num2str(m)]),eval(['P' num2str(n)]),...
      	eval(['log10(OF' num2str(i) ')']),10);
		colormap(gray);
		colorbar;
		xlabel(parameter_text(m),'FontSize',[10],'Fontname','times');
		ylabel(parameter_text(n),'FontSize',[10],'Fontname','times');
		set(gca,'position',[0.10 0.15 0.7 0.7], ...
         'fontname','times','fontsize',[10]);
      if m == 5
         set(gca,'xscale','log');
      elseif n == 5
         set(gca,'yscale','log');
      end
		set(colorbar,'position',[0.8625 0.3 0.025 0.4], ...
		   'FontSize',[10],'Fontname','times');
		set(gca,'units','normalized');
		h = axes('position',[0 0 1 1],'visible','off');
		set(gcf,'currentaxes',h);
      text(0.83,0.25,'log \phi [{\fontsize{8}-}]',...
         'FontSize',[10],'Fontname','times');
	drawnow;

end
