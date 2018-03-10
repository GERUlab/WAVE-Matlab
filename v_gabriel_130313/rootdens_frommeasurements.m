rootDAP = [0 20 40 68 120]
z1_param = [2.98 2.98 1.77 1.68 1.68]
z2_param = [1.36 1.36 5.24 7.57 7.57]

for t=1:120
    rootdensparam1 =interp1(rootDAP,z1_param,t);
    rootdensparam2 =interp1(rootDAP,z2_param,t);
    depth = 0:0.01:85
    rootdens = gampdf(depth,rootdensparam1,rootdensparam2);
    plot(depth,rootdens)
    drawnow
end