function draw_fluxes(tnode, potential_surface_flux,potential_transp, top_inf,trans,bot_inf,...
    cum_potential_surface_flux,cum_pot_transp,top_flux, cum_trans ,bot_flux, cum_infiltr, cum_evap,units, simplant,water_storage)

%Plot potential and actual fluxes from the water module

figure
plot(tnode, potential_surface_flux)
title ('potential surface flux','FontSize',14,'FontWeight','bold')
xlabel ([ 'time (' units{1,2} ')' ] )
ylabel ([ units{1,1} '/' units{1,2}])

if simplant
figure
plot(tnode, potential_transp)
title ('potential root water uptake','FontSize',14,'FontWeight','bold')
xlabel ([ 'time (' units{1,2} ')' ] )
ylabel ([ units{1,1} '/' units{1,2}])
end

figure
plot(tnode, top_inf)
title ('Actual surface flux','FontSize',14,'FontWeight','bold')
xlabel ([ 'time (' units{1,2} ')' ] )
ylabel ([ units{1,1} '/' units{1,2}])

if simplant
figure
plot(tnode, trans)
title ('Actual transpiration','FontSize',14,'FontWeight','bold')
xlabel ([ 'time (' units{1,2} ')' ] )
ylabel ([ units{1,1} '/' units{1,2}])
end

figure
plot(tnode, bot_inf)
title ('Bottom flux','FontSize',14,'FontWeight','bold')
xlabel ([ 'time (' units{1,2} ')' ] )
ylabel ([ units{1,1} '/' units{1,2}])


figure
plot(tnode, cum_potential_surface_flux(2:end))
title ('cumulative potential surface flux','FontSize',14,'FontWeight','bold')
xlabel ([ 'time (' units{1,2} ')' ] )
ylabel (units{1,1} )

if simplant
figure
plot(tnode, cum_pot_transp(2:end))
title ('cumulative potential transpiration','FontSize',14,'FontWeight','bold')
xlabel ([ 'time (' units{1,2} ')' ] )
ylabel (units{1,1} )
end

figure
plot(tnode, -top_flux)
title ('cum actual surface flux','FontSize',14,'FontWeight','bold')
xlabel ([ 'time (' units{1,2} ')' ] )
ylabel (units{1,1})

if simplant
figure
plot(tnode,cum_trans(2:end) )
title ('cum actual transpiration','FontSize',14,'FontWeight','bold')
xlabel ([ 'time (' units{1,2} ')' ] )
ylabel (units{1,1})
end

figure
plot(tnode, bot_flux)
title ('Cumulative bottom flux','FontSize',14,'FontWeight','bold')
xlabel ([ 'time (' units{1,2} ')' ] )
ylabel (units{1,1})

figure
plot(tnode, cum_infiltr(2:end))
title ('cumulative infiltration','FontSize',14,'FontWeight','bold')
xlabel ([ 'time (' units{1,2} ')' ] )
ylabel (units{1,1})

figure
plot(tnode, cum_evap(2:end))
title ('cumulative evaporation','FontSize',14,'FontWeight','bold')
xlabel ([ 'time (' units{1,2} ')' ] )
ylabel ( units{1,1} )

figure;
plot(tnode,water_storage);
grid on;
title('Water storage','FontSize',14,'FontWeight','bold');
ylabel(['water storage in profile(', units{1,1}, ')']);xlabel(['Time (',units{1,2},')']);

figure
storage2 = water_storage - water_storage(1); %Change of water stock in profile since beginning
error = top_flux + bot_flux - cum_trans(2:end) - storage2; %mass balance error
plot(tnode, bot_flux, tnode, top_flux, tnode, cum_trans(2:end), tnode, storage2, tnode, error)
title ('Water balance', 'Fontsize', 14, 'FontWeight', 'bold');
legend('bottom flux', 'top flux', 'transpiration', 'delta water storage', 'mass balance error')
ylabel(['cumulative balance terms(', units{1,1}, ')']);xlabel(['Time (',units{1,2},')']);
