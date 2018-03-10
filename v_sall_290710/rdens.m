



if isucr;
    rdens = 1000*rldf;
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [fractions]=calc_rt_distr (ira, irz, t, dx, rna,drz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
i=[1:ncs-1];
    if fraction(i) ~= 0;
        rdens(i+1) = w0_rdens* fractions(i+1)/fractions(i);
    end;
end


