
call calc_rt_distr(fractions)
dimension fractions(*)

      call calc_dim_roots(ira,irz,rna)
      sum_of_rt_distr = 0.d0
      if (ira.ne.irz) then
		sum_of_rt_distr = rt_distr(ira)* (rna - (x(ira)-dx/2.))/dx
		do i = ira+1,irz -1
			sum_of_rt_distr = sum_of_rt_distr + rt_distr(i)
		enddo
c		for last compartment correct for partial uptake
		sum_of_rt_distr = sum_of_rt_distr+rt_distr(irz)*
     $                    ((dx/2.)+x(irz)-drz(nday))/dx
      else
		sum_of_rt_distr = (rna - drz(nday))/dx
      endif

      if (sum_of_rt_distr.gt.0.d0) then
		do i = 1, ira-1
			fractions(i) = 0.d0
		enddo
		do i = ira,irz
			fractions(i) = rt_distr(i)/sum_of_rt_distr
		enddo
c		for first & last compartment correct for partial uptake
		if (ira.ne.irz) then
			fractions(ira) = (rt_distr(ira)*(rna - (x(ira)-dx/2.))/dx)
     $         /sum_of_rt_distr
			fractions(irz) = (rt_distr(irz)*((dx/2.)+x(irz)-drz(nday))
     $                    /dx)/sum_of_rt_distr
		else
			fractions(ira) = 1.d0
		endif
      else
		do i = 1, irz
			fractions(i) = 0.d0
		enddo
      endif
      return
      end  