function draw_time(bc)

%Draw time information

figure
plot (bc(:,1), bc(:,9))
xlabel ('time')
ylabel('timestep')
title ('time step information')

figure
plot (bc(:,1), bc(:,10))
xlabel ('time')
ylabel('number of iterations')
title ('iteration information')
