%%Sucros

function [dvs,rtdep,slai,slaig,wcrn,wkob,wiv,wrt,wso,wst]=sucros(amx,asrqso,dx,eai,ear,eff,harvest_date,hsh,ncrop,ncs,t,nsi,ntabel,plant_date,...
    rgr,rkdf,rlaicr,rlat,rldf,rmainso,rmatr,rmrd,scp,sla,specweig,ssl,t,table,tbase,temp,tmax,tmin)

%c%###################################################################################
      subroutine sucros
c     in   : amx, asrqso, dx, eai, ear, eff, harvest_date, hsh, ncrop, ncs,
c            nday, nsl, ntabel, plant_date, rgr, rkdf, rlaicr, rlat, rldf,
c            rmainso, rmatr, rmrd, scp, sla, specweig, ssl, t, table,
c            tbase, temp, tmax, tmin
c     out  : dvs, rtdep, slai, slaig, wcrn, wkob, wlv, wrt, wso, wst
c     calls: afgen, calc_nitreductgrow, calc_watreductgrow, report_err, stop_simulation
%c###################################################################################
%dimension gsdst(3),gswt(3)
%dimension fractions(kt_comps)
%parameter (pi = 3.1416d0, rd = pi/180.d0)
starts = 200.
%save gsdst, gswt, slaid,   tadrw
%save  cgphot,  wlvd,  cnphot
%data init /.true./
gsdst=[0.112702 0.5 0.887298d0];
gswt= [0.277778, 0.444444d0, 0.277778];
slaid = 0.0;
tadrw=0; tsum=0; cgphot=0;
wlvd =0;
cnphot=0; % rldf /0.0d0, kt_comps*0.0d0/
sucr_red_wat = calc_watreductgrow()
sucr_red_nit = calc_nitreductgrow()
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t == harvest_date | t > harvest_date
    slaid = 0; plai = 0; slaig = 0 slai = 0; rtdep = 0; wlvg= 0; wlvd= 0; wst= 0 ;wso= 0;
    wrt= 0; 	wlv= 0; wcrn= 0; wkob= 0; tadrw= 0; dvs = 0; sucr_red_wat=0; sucr_red_nit=0;
end

%if(dvs.gt.2.0d0.and.ncrop.lt.4)return

%     calculate average and effective temperature for growth (dteff)
tmina=tmin(nday)
tmaxa=tmax(nday)
tempa=0.5*(tmaxa+tmina)
tminr=tmina
tmaxr=tmaxa
if tminr< 0
   tminr=2 +tminr*0.805;
end
if tmaxr< 0
   tmaxr=2.d0+tmaxr*0.805
end
   ddtmp=tmax(nday)-(tmax(nday)-tmin(nday))/4.0
   dteff= dmax1(tempa-tbase, 0.0d0)

%	temperature sums
%   winter wheat(emergence day and temperature sum since the first of january)
if ncrop == 1
    if(idint(t).gt.365)then
			tsum=tsum+dteff
			if (init) then
				init = .false.
				plai = nsl * ssl
				slaig = plai
				slai = plai
				sla=afgen(tab10,0.d0,ntabel(10))
				wlvg = slaig/sla
				wlv = wlvg
				rtdep = 3.5d0
			endif
			if(tsum.lt.starts) return
		else
			return
		endif      
c	other crops (temperature sum since planting day)
      else
		if(idint(t).lt.plant_date)then
			return
		else
			if (init)then
				init =.false.
				plai = nsl * ssl
				if (plai.eq.0.d0) then
					call report_err 
     $					('the initial lai is zero for the crop growth model')
					call stop_simulation 
     $				('programme stopped: check the err_file')
				endif
				slaig = plai
				slai = plai
				if(ncrop.eq.3) sla=0.035
				wlvg = slaig/sla
				wlv = wlvg
				rtdep = 3.0d0
			endif
			tsum=tsum+dteff
		endif
      endif
c     amax as a function of dvs and temp
      if(ncrop.eq.4)then
		amdvs=1.0d0
      else 
		amdvs=afgen(tab1,dvs,ntabel(1))
      endif
      amtmp=afgen(tab2,ddtmp,ntabel(2))
      amax=amx*amdvs*amtmp
c     daylength
      dec = -23.4d0* dcos(2.d0*pi*(t+10.d0)/365.d0)
      sinld = dsin(dec*rd)* dsin(rlat*rd)
      cosld = dcos(dec*rd)* dcos(rlat*rd)
      dl = 12.d0*(pi+2.d0* dsin(sinld/cosld))/pi
c     daily integral of sine of solar inclination(sinb)
      dsinb = 3600.d0*(dl*sinld+24.d0*cosld*sqrt(1.d0-(sinld/cosld)**2)/pi)
      dsinbe = 3600.d0*(dl*(sinld+0.4d0*(sinld*sinld+0.5d0*cosld*cosld))+
     $	12.d0*cosld*(2.d0+3.d0*0.4d0*sinld)*
     $	sqrt(1.d0-(sinld/cosld)**2)/pi)
c	corrected solar constant (j/m2/d)
      sc=1370.d0*(1.d0+0.033d0* dcos(2.d0*pi*t/365.d0))
c     daily extra terrestrial radiation (j/m2/d) 
      dso=sc*dsinb
c     daily photosynthetically active radiation (j/m2/d)
      dpar=0.50*hsh(nday)*10000.d0
c	atmospheric transmission
      atmtr=hsh(nday)*10000.d0/dso
c     fraction diffusive radiation (frdf) 
      if(atmtr.le.0.07d0)then
		frdf=1.d0
      else if(atmtr.le.0.35d0)then
		frdf=1.d0-2.3d0*(atmtr-0.07d0)**2
      else if(atmtr.le.0.75d0) then
		frdf=1.333d0 -1.46d0*atmtr
      else  
		frdf=0.23d0
      endif
c     daily gross assimilation of the canopy (dtga, kg co2/ha/d)
      dtga=0.d0
      do i = 1,3
		hour=12.d0+dl*0.5d0*gsdst(i)
c		instantaneous radiation above the canopy 
		sinb= dmax1(0.d0,sinld+cosld*cos(2.*pi*(hour+12.)/24.))
c		diffus par (pardf) and direct par (pardr) (j/m2/s)
		par=dpar*sinb*(1.0d0+0.4d0*sinb)/dsinbe
		pardf=dmin1(par,frdf*dpar*sinb/dsinb)
		pardr=par-pardf
c		instantaneous gross assimilation at different canopy depths
		fgros=0.0d0
c		radiation profile in the canopy and instantaneous absorbed radiation for 
c         succesive leaf layers
		do j = 1,3
			rlaic =plai*gsdst(j)
c			canopy refelection coefficient (refl) as a function of the leaf scattering coefficient (sc)
			refl=(1.-sqrt(1.-scp))/(1.+sqrt(1.-scp))         
c			extinction coefficient for direct component (rkbl) and total direct flux 
c			(akdrt) / cluster factor ad ratio between emperical and theoretical value of
c			(rkdf)
			clustf=rkdf/(0.8d0*sqrt(1.d0-scp))
			rkbl=(0.5d0/sinb)*clustf
			akdrt=rkbl*sqrt(1.d0-scp)
c			diffus flux,total direct flux and direct component of direct flux
			parldf=(1.d0-refl)*pardf*rkdf* dexp(-rkdf*rlaic)
			parlt=(1.d0-refl)*pardr*akdrt* dexp(-akdrt*rlaic)
			parldr=(1.d0-scp)*pardr*rkbl* dexp(-rkbl*rlaic)
c			absorbed fluxes (j/m2 leaf/ s ) for shaded and sunlit lleaves
			parlsh=parldf+(parlt-parldr)
c			direct par absorbed by leaves perpendicular on direct beam
			parlpp=pardr*(1.d0-scp)/sinb
c			fraction of the sunlit leaf area
			fslla= dexp(-rkbl*rlaic)*clustf
c			assimilation of shaded leaf area (kg co2/ha leaf/hour)
			asssh=amax*(1.d0- dexp(-eff*parlsh/amax))
c			assimilation of sunlit leaf area (kg co2/ha leaf/hour)
			asssl=0.0d0
			do k=1,3
				parlsl=parlsh+parlpp*gsdst(k)
				asssl=asssl+amax*(1.d0- dexp(-parlsl*eff/amax))*gswt(k)
			enddo
c			hourly total gross assimilation (kg co2/ ha soil/h) 
			fgros=fgros+((1.d0-fslla)*asssh+fslla*asssl)*plai*gswt(j)
		enddo
c		integration of instantaneous assimilation to a daily total (dtga)
		dtga=dtga+fgros*dl*gswt(i)
      enddo
c	water and nutrient stress
      sucr_red_wat = calc_watreductgrow()
      sucr_red_nit = calc_nitreductgrow()
c	glucose production 
      gphot=dtga*(30.0d0/44.0d0)* dmin1(sucr_red_wat, sucr_red_nit)

c     calculation of the crop development rate
c     sugar beets
      if(ncrop.eq.5)then
		ratdvs= dmin1(19.0d0,( dmax1(tempa-2.0d0,0.d0)))
c     potatoes     
      else if(ncrop.eq.4)then
		if(tempa.lt.13.0d0)then
			term=tempa-2.0d0
		else
			term=29.0d0-tempa
		endif
		ratdvs= dmin1(11.d0, dmax1(term,0.d0))
c     spring,winter wheat and mais
      else if(ncrop.lt.4)then
		if(dvs.lt.1)then
			ratdvs=afgen(tab3,tempa,ntabel(3))
		else
			ratdvs=afgen(tab4,tempa,ntabel(4))
		endif                            
      endif

c     maintenance requirements
      teff=2.d0**((tempa-25.0d0)/10.d0)
      if(wlv.le.0.0d0)then
		rmndvs=0.0d0
      else
		rmndvs=wlvg/wlv
      endif
c	winter whear, spring wheat, potatoes
      rmaints=0.03d0*wlv+0.015d0*wst+0.015d0*wrt+rmainso*wso
c     sugar beets
      if(ncrop.eq.5) then
		rmaints=rmaints+0.015*wcrn
c	mais	
      else if(ncrop.eq.3) then
		rmaints=rmaints+0.015*wkob 
      endif
      rmaint= dmin1 (rmaints*teff*rmndvs, gphot)

c     growth fraction of crop components
      cgphot=cgphot+gphot
      cnphot=cnphot+gphot-rmaint
c     potatoe
      if(ncrop.eq.4)then
		ind=1.0d0/(0.0015d0+0.00079d0*rmatr)
		opt=(dvs-ind)/430.0d0
		fsh= dmin1(1.0d0, dmax1(0.8d0 + 0.2d0*opt,0.8d0))
		frt=1.0d0-fsh
		flv= dmin1(0.75d0, dmax1(0.75d0 -opt,0.d0))
		fso= dmin1(1.0d0, dmax1(opt,0.d0))
		fst=1.0-flv-fso
c     spring wheat, winter whear, mais, sugar beet 
      else
		fsh=afgen(tab5,dvs,ntabel(5))
		flv=afgen(tab6,dvs,ntabel(6))
		fst=afgen(tab7,dvs,ntabel(7))
c		spring, winter wheat 
		if(ncrop.le.2)then
			fso=1.0d0-flv-fst
			frt=1.0d0-fsh
c		sugarbeets
		else if(ncrop.eq.5)then
			fcrn=1.0d0-flv-fst
			frt=afgen(tab8,dvs,ntabel(8))*(1.0d0-fsh)
			fso=1.0d0-fsh-frt
c         maize
		else if(ncrop.eq.3)then
			fcob=afgen(tab8,dvs,ntabel(8))
			fso=1.0d0-flv-fst-fcob
			frt=1.0d0-fsh
		endif
      endif

c     dry matter production
c     maize
      if(ncrop.eq.3) then
		asrq = fsh*(1.46*flv+1.51*fst+1.51*fcob+asrqso*fso)
     $		+1.44*frt
c     sugarbeets
      else if(ncrop.eq.5) then
		asrq=fsh*(1.46d0*flv+1.51d0*fst+1.51d0*fcrn)
     $		+asrqso*fso+1.44d0*frt
c     winter wheat, spring wheat, potatoes 
      else
		asrq=fsh*(1.46d0*flv+1.51d0*fst+asrqso*fso)+1.44d0*frt
      endif

c     dry matter growth rates
      ratwtot=(gphot-rmaint)/asrq
      ratwrt=ratwtot*frt
      ratwsh=ratwtot*fsh
      ratwlvg=ratwsh*flv
      ratwst=ratwsh*fst
c     winter whear, spring wheat, potatoes 
      ratwso=ratwsh*fso
c     sugarbeets
      if(ncrop.eq.5)then 
		ratwso=ratwtot*fso
		ratwcrn=ratwsh*fcrn
c     maize
      else if(ncrop.eq.3)then
		ratwcob=ratwsh*fcob
      endif

c     lai and leaves  death rates

c	rdrdv: senescence death rate factor
c     spring and winter wheat
      if(ncrop.le.2)then
		if(dvs.lt.1)  then
			rdrdv=0.0d0
		else
			rdrdv=afgen(tab11,tempa,ntabel(11))
		endif
		if(ncrop.eq.1) rdrdv=rdrdv*afgen(tab9,dvs,ntabel(9))
c     sugar beets
      else if(ncrop.eq.5) then 
		rdrdv=afgen(tab11,dvs,ntabel(11))*ratdvs
c     maize
      else if(ncrop.eq.3) then
		if(dvs.lt.1.35d0) then
			rdrdv=0.0005d0
		else
			rdrdv=0.003d0
		endif
		rdrdv=rdrdv*(tempa-8.0d0)
c     potato 
      else if(ncrop.eq.4)then
		tssnc=tsum-725.0d0
		if(tssnc.lt.0.0d0)then
			rdrdv=0.0d0
		else 
			term= dmax1(tempa-2.0d0,8.0d0)
			rdrdv=term* dexp(-11.7d0+0.68d0*rmatr)* dexp(tssnc*(0.0068d0-
     $		0.0006d0*rmatr))
		endif
      endif
c     rdrsh: self shading death rate factor
      term=0.03d0*(plai-rlaicr)/rlaicr
      rdrsh = dmax1 (0.0d0, dmin1(term,0.03d0))
c     rdrlt: chilling temperatures death rate factor (maize only)
      if(ncrop.eq.3) then
		if(dvs.lt.1.25d0) then
			rdrlt=0.d0
		else
			term=(6.0d0-tempa)/6.0d0
			rdrlt =  dmin1(1.0d0, dmax1(0.0d0,term))
		endif
      endif
c     relative death rates of leaves
      rdrwlvd = dmax1(rdrsh,rdrdv)
c	maize
      if(ncrop.eq.3)then
		if(dvs.lt.1)then 
			rdrwlvd=0.0d0
		else 
			rdrwlvd= dmax1(rdrlt,rdrwlvd)
			rdrwlvd= dmax1(0.001d0, rdrwlvd)
		endif
      endif
c     leaf area death rate 
      ratlaid=slaig*( dexp(rdrwlvd)-1.d0)
      if(plai.le.0.0d0.or.wlvg.le.0.0d0)then
		ratwlvd=0.d0
      else
		ratwlvd=wlvg*ratlaid/slaig
      endif

c     ear area growth rates
c	maize, potato, sugar beet  
      if(ncrop.gt.2)then
		rateai=0.0d0
c	winter wheat, spring wheat
      else 
		if(dvs.lt.0.8d0)then
			rateai=0.0
		else
			rateai=ear*tadrw
		endif
		if (dvs.ge.1.3d0) rateai = rateai - rdrwlvd*eai
      endif

c	leaf area growth rate
      sla=afgen(tab10,dvs,ntabel(10))
c	potatoe, sugar beet
      if(ncrop.ge.4)then
		if(tsum.le.450.0d0.and.plai.le.0.75d0)then
			ratlaig=slaig*( dexp(rgr*dteff)-1.0d0)
		else
			ratlaig=sla*ratwlvg
		endif
c	winter wheat, spring wheat, maize
      else
		if(dvs.lt.0.3d0.and.plai.lt.0.75d0)then
			ratlaig=slaig*( dexp(rgr*dteff)-1.0d0)
		else
			ratlaig=sla*ratwlvg
		endif
      endif

c     root penetration rates and root density growth rates
      rlnew=ratwrt*specweig*1000.d0
      if (rlnew.ne.0.d0) then
c		effective temperature for root growth
		tempcr=(tmaxr+tminr)/2.0
		dtt=tempcr
		if(tmaxr.lt.0.0d0) then
			dtt=0.d0
		elseif(tminr.le.0.0d0) then
			tcor= dmax1(0.d0,0.1583d0*(tmaxr-tminr)- dabs(dtt)*0.4043d0)
			dtt=tempcr+tcor
		endif
c		rooting depth (rtdep = cm)
		daincr =  dmin1(dtt *0.22d0 ,1.8d0)
		rtdep= dmin1(rtdep+daincr, dabs(rmrd/10.d0))
c		calculate the  root density (rldf)
		dx_cm=dx/10.d0
		call rt_distr_sucr(fractions,rtdep)
		do l = 1, ncs
			rldf(l)=rldf(l)+ fractions(l)*rlnew/(100000000.d0*dx_cm)
		enddo
      endif

c     calculation of the integrals

c     plant development stage
      dvs=dvs+ratdvs

c     dry mass in the plant organs
      wlvg=wlvg+(ratwlvg-ratwlvd)
      wlvd=wlvd+ratwlvd
      wst=wst+ratwst
      wso=wso+ratwso
      wrt=wrt+ratwrt
      wlv=wlvg+wlvd
      tadrw=wlv+wst+wso
c	sugar beet: crown weight
      if(ncrop.eq.5) then
		wcrn=wcrn+ratwcrn
		tadrw=tadrw+wcrn
c	maize: cob weight
      else if(ncrop.eq.3) then
		wkob=wkob+ratwcob
		tadrw=tadrw+wkob 
      endif

c     leaf and ear area charcteristics 
      slaig= dmax1(slaig+(ratlaig-ratlaid),0.d0)
      slaid=slaid+ratlaid
      eai =  dmax1(eai + rateai,0.d0)
      plai=slaig+eai/2
      slai=slaig+slaid
      return
      end


c###################################################################################
      double precision function  root_sucros ()
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'crop.com'
      root_sucros = -rtdep*10.d0
      return
      end


c###################################################################################
      double precision function rlai_sucros ()
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'crop.com'
      rlai_sucros = slaig
      return
      end


c###################################################################################
      double precision function dvs_sucros ()
c###################################################################################
      implicit double precision (a-h,o-z)
      include   'constant'
      include   'gen.com'
      include   'crop.com'
      dvs_sucros = dvs
      return
      end

c###################################################################################
      subroutine rt_distr_sucr(fractions,rtdep)
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'wat.com'
      dimension fractions(*)
      drz(nday) = -rtdep*10.d0
      call calc_rt_distr (fractions)
      end

c###################################################################################
      subroutine rdens_sucros(rdens)
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'crop.com'
      dimension rdens(kt_comps)
      do i=1,ncs
        rdens(i) = 1000.d0*rldf(i)
	enddo
      end

c###################################################################################
      double precision function calc_xncle ()
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'crop.com'
      dimension tab12(2,kt_cropinfo)
      equivalence (tab12(1,1), table(1,1,12))

      calc_xncle = afgen(tab12,dvs,ntabel(12))
      return
      end

c###################################################################################
      double precision function calc_xncst ()
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'crop.com'
      dimension tab13(2,kt_cropinfo)
      equivalence (tab13(1,1), table(1,1,13))

      calc_xncst = afgen(tab13,dvs,ntabel(13))
      return
      end

c###################################################################################
      double precision function calc_xncrt ()
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'crop.com'
      dimension tab14(2,kt_cropinfo)
      equivalence (tab14(1,1), table(1,1,14))

      calc_xncrt = afgen(tab14,dvs,ntabel(14))
      return
      end

c###################################################################################
      double precision function calc_xncso ()
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'crop.com'
      dimension tab15(2,kt_cropinfo)
      equivalence (tab15(1,1), table(1,1,15))

      calc_xncso = afgen(tab15,dvs,ntabel(15))
      return
      end

c###################################################################################
      double precision function sucr_anlv ()
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'crop.com'
      dimension s_ancl(kt_crop)
      data s_ancl/0.05,0.05,0.05,0.05,0.05/

      sucr_anlv = s_ancl(ncrop)*ssl*nsl/sla
      return
      end

c###################################################################################
      subroutine plant_weights(wso_out,wlv_out,wrt_out,wst_out,wlvg_out)
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include  'crop.com'

      wso_out = wso
      wlv_out = wlv
      wlvg_out = wlvg
      wrt_out = wrt
      wst_out = wst
      end

c###################################################################################
      double precision function calc_nitreductgrow ()
c###################################################################################
      implicit double precision (a-h,o-z)
      include  'constant'
      include  'gen.com'
      include 'crop.com'
      common /local_nit_upt/ ancl, rlncl,rmncl

      if (idint(t).lt.plant_date.or.idint(t).gt.harvest_date) then
		calc_nitreductgrow=0.d0
		return
      endif
      if (simnit) then
c		potatoes and sugar beet
		if(ncrop.ge.4) then
			if(tsum.le.450.d0.and.plai.lt.0.75d0) then
				calc_nitreductgrow=1.d0
			else
				calc_nitreductgrow=
     $			dmin1(1.d0, dmax1(0.d0,((ancl-rlncl)/(rmncl-rlncl))))
			endif
c		winter wheat, spring wheat and maize
		else
			if(dvs.le.0.3d0.and.plai.lt.0.75d0) then
				calc_nitreductgrow=1.d0
			else
				calc_nitreductgrow=
     $			dmin1(1.d0, dmax1(0.d0,((ancl-rlncl)/(rmncl-rlncl))))
			endif
		endif
      else
		calc_nitreductgrow= 1.d0
      endif 
      return
      end

c###################################################################################
      function afgen(table,var,ndp)
c###################################################################################
c     x values must be in increasing order.
c     ndp = number of data pairs in the interpolation table
      implicit double precision (a-h,o-z)
      dimension table(2,ndp)
      if(var.lt.table(1,1)) then
c		x coordinate below range. assume gradient below 1st point = 0
		afgen=table(2,1)
      else if(var.gt.table(1,ndp)) then
c		x coordinate above range. assume gradient above last point = 0
		afgen=table(2,ndp)
      else
		n=1
10		n=n+1 
		if(var.gt.table(1,n)) goto 10
c		endvar lies between table(1,n-1) and table(1,n)
c		linear interpolation
		afgen=((var-table(1,n-1))*((table(2,n)-table(2,n-1))
     $	/(table(1,n)-table(1,n-1))))+table(2,n-1)
      endif
      return
