*72*********************************************************************
      subroutine radgen_init(CTARGET,ebeam,LST40,ixytst)
      implicit none 
c
#include "radinclude/tailcom.inc"
#include "radinclude/mcpolrad.inc"
#include "radinclude/cmpcom.inc"
#include "radinclude/radgenkeys.inc"
      integer ige,ixytst
      real ebeam
      INTEGER LST40
      CHARACTER*2 CTARGET
      CHARACTER*13 tname
c
*	open(734,file='kill_elas_res.dat',status='old')
*	read(734,*)kill_elas_res
*	close(734)

*     kill_elas_res =2, no rad corrections
*     kill_elas_res =0, all rad corrections
      kill_elas_res=0

*     initialise the rad. correction part !
*     =2 : fint cernlib extrapolation routine
*     =1 : 2dim spline (is not supported in this version)
*     =0 : make lookuptable
*     = -1:	do not lookup table , calculate all events
       ixytest=ixytst
*----------------------------
*    =2 : polrad 2.0
      ige= 2
*----------------------------
*     ire is target 1,2,3
      IF(CTARGET.eq.'P') THEN
	 ire = 1
      ELSEIF(CTARGET.eq.'D') THEN
	 ire = 2
      ELSEIF(CTARGET.eq.'HR'.or.CTARGET.eq.'HI') THEN
	 ire = 3
      ELSE
         print *,'UNKNOWN TARGET!!!',CTARGET
         stop
      ENDIF
*----------------------------
*     plrun : beam polarisation    [-1.,1.]
*     pnrun : nucleon polarisation     [-1.,1.]

      IF (LST40.gt.0)  THEN
*  larger cross section state (anti parallel )
	 plrun = -1.
	 pnrun = 1.
	 tname='xytab0ant.dat'
      ELSEIF (LST40.lt.0) THEN
*  smaller cross section state ( parallel )
	 plrun = -1.
	 pnrun = -1.
	 tname='xytab0par.dat'
      ELSE
	 plrun = 0.
	 pnrun = 0.
	 tname='xytab0unp.dat'
      ENDIF
      write(tname(6:6),'(i1)')ire

*----------------------------
* grid of important regions in theta (7*ntk)
      ntk = 35
*----------------------------
* photonic energy grid
      nrr = 100
*----------------------------
* min. energy in the calo (resolution parameter)
      demin=0.1
*----------------------------
       ap=2.*amp
       amp2=amp**2
       ap2=2.*amp**2


c	print *,' kill...', kill_elas_res,amc2
       if(kill_elas_res.eq.1)amc2=4d0
c	print *,' kill...', kill_elas_res,amc2

      if(ire.eq.1)then
       amt=.938272d0
       tara=1d0
       tarz=1d0
       fermom=0d0
      elseif(ire.eq.2)then
       amt=1.87561d0
       tara=2d0
       tarz=1d0
       fermom=.07d0
      elseif(ire.eq.3)then
       amt=2.80833d0
       tara=3d0
       tarz=2d0
       fermom=.164d0
      endif
       if(ixytest.ge.0)call
     .	xytabl(tname,ebeam,plrun,pnrun,ixytest,ire)
      end
*********************************************************




C-----------------------------------
      SUBROUTINE RADGEN(e1,VPGEN,VPRAD,PhRAD,q2tr,anutr,WEIGHT)
      IMPLICIT NONE
***********************************************************
*  INPUT :
*	 VPGEN(4) : virtual photon 4 vector generated
*  OUTPUT :
*	 VPRAD(4) : virtual photon 4 vector   _TRUE_
*	 PhRAD(4) : real rad. photon 4 vector
*	 q2tr	  : true Q2
*	 Utr	  : true U
***********************************************************
#include "radinclude/mcpolrad.inc"
#include "radinclude/phiout.inc"
#include "radinclude/tailcom.inc"
#include "radinclude/radgenkeys.inc"
      REAL VPGEN(4),VPRAD(4),PhRAD(4),q2tr,anutr
      REAL x,y,e1,q2,xs,ys,phil,weight
      INTEGER itagen
*

*     calculate -y,x from VPGEN(4)

      Q2=vpgen(1)**2+vpgen(2)**2+vpgen(3)**2-vpgen(4)**2
      ys=vpgen(4)/e1
      xs=q2/(2.*0.938272*ys*e1)

      phil=asin(-vpgen(2)/sqrt(vpgen(1)**2+vpgen(2)**2))
c      print *,'in radgen:e1,Q2,ys,xs',e1,Q2,ys,xs,vpgen
c      stop

      if(ixytest.ge.0)then
*     output		 itagen
       if(kill_elas_res.ne.2)then
	call itafun(ys,xs,plrun,pnrun,ixytest,itagen)
       else
	itagen=0
	sigcor=1.
       endif
      else
	itagen=-1
      endif
c
      if(itagen.ne.0.and.abs(ixytest).ne.2)then
	iphi=0
	call mpolrad(e1,ys,xs,1.,plrun,pnrun,itagen)
	if(ixytest.eq.1)ixytest=-2
      endif
      call radgam_pol(e1,ys,xs,phil,ixytest,itagen,q2tr,anutr)


* Perform calculation of VPRAD(*),PhRAD(*)

      phrad(1)=sngl(dplabg(1))
      phrad(2)=sngl(dplabg(2))
      phrad(3)=sngl(dplabg(3))
      phrad(4)=ys*e1-anutr

      vprad(1)=vpgen(1)-phrad(1)
      vprad(2)=vpgen(2)-phrad(2)
      vprad(3)=vpgen(3)-phrad(3)
      vprad(4)=vpgen(4)-phrad(4)

* set the weight
      if(kill_elas_res.eq.2)then
	  weight=1.
      else
	  weight = sigcor
      endif

      return
      end




      subroutine phidst_pol
*     -----------------
*
*     calculate phi-distribution of bremsstrahlung photons
*     needed for correction of hadronic distributions
*
#include "radinclude/double.inc"
c#include "radinclude/concom.inc"
c#include "radinclude/kincom.inc"
#include "radinclude/phiout.inc"
c#include "radinclude/gamcom.inc"
#include "radinclude/mcpolrad.inc"
      common/p/dpi,dpi2,dalfa,i1(8),i2(8)
*
      dimension db(3),dx(31),dy(31),dz(31)
      external dfufi_pol
      data ndim /31/
*
*     boundaries for 2 integrations
*
      db(1)=0.d0
      db(2)=dpi/10.d0
      db(3)=dpi
*
      k=0
      dsum=0.d0
      do i=1,2
	 dbnd1=db(i)
	 dbnd2=db(i+1)
	 call radquad_pol(dfufi_pol,dbnd1,dbnd2,dx,dy,dz,dh,ndim)
	 do j=1,ndim
	    if(i.eq.1.or.j.ne.1) k=k+1
	    dphi(k)=dx(j)
	    ddephi(k)=dh
	    dsumph(k)=dsum+dz(j)
	 enddo
	 dsum=dsum+dz(ndim)
      enddo
      kmp=k
      return
      end




      subroutine radgam_pol(erad,yrad,xrad,genphi,
     . ixytest,itagen,q2tr,anutr)
c      implicit none
#include "radinclude/double.inc"
*#include "radinclude/concom.inc"
*#include "radinclude/kincom.inc"
*#include "radinclude/radcom.inc"
#include "radinclude/phiout.inc"
*#include "radinclude/gamcom.inc"
#include "radinclude/mcpolrad.inc"
#include "radinclude/cmpcom.inc"
#include "radinclude/sxycom.inc"
#include "radinclude/tailcom.inc"
       parameter(nt=35)
       real*4 denstk
       common/density/denstk(nt,nt,245,3),width(nt,nt,7,3)
       real*8 taa,tm,sfm0,podinl
      common/amf2/taa,tm(8,6),sfm0(8)
*
      dimension rnd(1)
C	double precision dreal
*
*     selection of mass dmj and mass index iadmj
*
      call grndm(rnd,1)
      r1=rnd(1)
      if(ixytest.lt.0)then
      rr1=r1*(tbor+tine+tpro+tnuc)
c     print *,tine,tpro,tnuc
	 if (rr1.gt.(tine+tpro+tnuc)) then
	       ita=0
	 elseif(rr1.gt.(tpro+tnuc)) then
	       ita=1
	       sicurr=rr1-tpro-tnuc
C	       call conk2(dreal(erad),dreal(xrad),dreal(yrad),1)
	       call conk2(dble(erad),dble(xrad),dble(yrad),1)
	 elseif(rr1.gt.tnuc)then
	       ita=3
	       sicurr=rr1-tnuc
C	       call conk2(dreal(erad),dreal(xrad),dreal(yrad),1)
	       call conk2(dble(erad),dble(xrad),dble(yrad),1)
	 else
	       ita=2
	       sicurr=rr1
C	       call conk2(dreal(erad),dreal(xrad),dreal(yrad),2)
	       call conk2(dble(erad),dble(xrad),dble(yrad),2)
	 endif
       else
	 xs=dble(xrad)
	 ys=dble(yrad)
	 ita=itagen
	 if(ita.eq.1)sicurr=r1*tine
	 if(ita.eq.2)sicurr=r1*tnuc
	 if(ita.eq.3)sicurr=r1*tpro
	 if(ita.eq.2)then
	 call conk2(dble(erad),dble(xrad),dble(yrad),2)
	 else
	 call conk2(dble(erad),dble(xrad),dble(yrad),1)
c	 print *,an,xs,ys,ita,s,x
	 endif
	 isf1=1
	 isf2=isf20
	 isf3=1

c	 print *,an,xs,ys,erad,xrad,yrad
c	  sicurr=scgen
       endif

	 if(ita.ne.0)then
        call grndm(rnd,1)
        r1=rnd(1)*dsitkm(ndxtkm(ita),ita)
	 do ktk=1,ndxtkm(ita)
c	     print *,dsitkm(ktk,ita),r1,ndxtkm(ita)
	    if (r1.le.dsitkm(ktk,ita)) then
	       if (ktk.eq.1) then
		  dtk1=r1/dsitkm(ktk,ita)
	       else
		  dtk1=(r1-dsitkm(ktk-1,ita))/
     +			      (dsitkm(ktk,ita)-dsitkm(ktk-1,ita))
	       endif
	       dtk=dcvtkm(ktk,ita)+(dtk1-0.5d0)*ddetkm(ktk,ita)
	       goto 30
	    endif
	 enddo
	 dtk=dcvtkm(ndxtkm(ita),ita)
   30	 continue
	 taout=(sx-sqly*cos(dtk))/ap2
c	  print *,taout,dtk,ita,ktk,r1

c	  print *,sx,sqly,ap2
	 if(ita.eq.1)then
	    taa=taout
	    iphi=0
	    call tails(taout,tm)
	    call grndm(rnd,1)
	    rcal=ap*demin
	    rmax=(w2-amc2)/(1.+taout)
	    do krr=1,nrr
	       ddeler(krr,ktk)=(rmax-rcal)/nrr
	       drcurr(krr,ktk)=rcal+ddeler(krr,ktk)*(krr-0.5)
c		print *,krr,ktk,drcurr(krr,ktk)
	      if(krr.eq.1)then
	      dsigmr(krr,ktk)=podinl(drcurr(krr,ktk))
	      else
	      dsigmr(krr,ktk)=dsigmr(krr-1,ktk)+podinl(drcurr(krr,ktk))
	      endif
c	       print *,dsigmr(krr,ktk),drcurr(krr,ktk),taout
	    enddo
	    r2=rnd(1)
	    dsircu=r2*dsigmr(nrr,ktk)
	    do krr=1,nrr
	       if (dsircu.le.dsigmr(krr,ktk)) then
		  if(krr.eq.1)then
		     drr1=dsircu/dsigmr(krr,ktk)
		  else
		     drr1=(dsircu-dsigmr(krr-1,ktk))
     .			/(dsigmr(krr,ktk)-dsigmr(krr-1,ktk))
		  endif
		  drr=drcurr(krr,ktk)+(drr1-0.5d0)*ddeler(krr,ktk)
		  goto 20
	       endif
	    enddo
	    drr=drcurr(nrr,ktk)
   20	    continue
	    dom=drr/ap
	 else
	    dom=(sx-y)/ap/(1.+taout)
	endif
	rrout=ap*dom

*
*---	 selection of phi angle
*
*
c	  print *,rrout,taout
c	 pause
	 iphi=1
	 call phidst_pol

*
	 call grndm(rnd,1)
	 r3=rnd(1)
	 dphpoi=r3*dsumph(kmp)
	 do kph=2,kmp
	    if (dphpoi.le.dsumph(kph)) then
	       dphk1=(dphpoi-dsumph(kph-1))/
     +			   (dsumph(kph)-dsumph(kph-1))
	       dphk=dphi(kph)+(dphk1-1.0d0)*ddephi(kph-1)
	       goto 40
	    endif
	 enddo
	 dphk=dphi(kmp)
   40	 continue
	 call grndm(rnd,1)
	 r4=rnd(1)
	 if (r4.gt.0.5) dphk=-dphk

*
*     radiative photon
*
      deg=dom
      dthg=dtk
      dphig=dphk
      sigma_total=sngl(tbor+tine+tpro+tnuc)

	 q2tr=y+rrout*taout
	 anutr=sx/ap-dom

	   dgpz=deg*dcos(dthg)
	   dgpxy=deg*dsin(dthg)
	   dgpx=dgpxy*dcos(dphig)
	   dgpy=dgpxy*dsin(dphig)
*
*---	   momentum components in the LAB-system:
*---	   two rotations needed - first within the scattering plane
*---	   around the y-axis and second around the new z-axis (beam
*---	   direction) by Phi of the scattering plane
*
	   dgplx=-dgpz*dsts+dgpx*dcts
	   dgply=dgpy
	   dgplz=dgpz*dcts+dgpx*dsts
	   dcphi=dcos(dble(genphi))
	   dsphi=dsin(dble(genphi))
	   dplabg(1)=dgplx*dcphi-dgply*dsphi
	   dplabg(2)=dgplx*dsphi+dgply*dcphi
	   dplabg(3)=dgplz
       else
	 q2tr=2d0*amh*erad*xrad*yrad
	 anutr=yrad*erad
c	      write(*,'(i5,5f8.3)')ita,xs,ys,q2tr,anutr,sigcor

	   dplabg(1)=0.0d0
	   dplabg(2)=0.0d0
	   dplabg(3)=0.0d0


      endif



      end

       double precision function dfufi_pol(dx)
*     -----------------------------------
*     needed for correction of hadronic distributions
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire
      common/amf2/taa,tm(8,6),sfm0(8)
      common/p/pi,pi2,alfa,i1(8),i2(8)
#include "radinclude/mcpolrad.inc"
*

      phipoi=dx
      taa=taout
       call tails(taout,tm)
       dfufi_pol=an*alfa/pi*podinl(rrout)

      return
      end




      subroutine radquad_pol(dfunct,dlow,dup,dx,dy,dz,dh,ndim)
*     -------------------------------------------------
*
#include "radinclude/double.inc"
      dimension dx(31),dy(31),dz(31)
      external dfunct
*
      dsum2=0.d0
      if (ndim.gt.1) then
	 dh=(dup-dlow)/float(ndim-1)
	 do i=1,ndim
	    daux=dlow+dh*float(i-1)
	    dx(i)=daux
	    dy(i)=dfunct(daux)
	 enddo
	 do i=2,ndim
	    dsum1=dsum2
	    dsum2=dsum2+.5d0*(dx(i)-dx(i-1))*(dy(i)+dy(i-1))
	    dz(i-1)=dsum1
	 enddo
	 dz(ndim)=dsum2
      elseif (ndim.eq.1) then
	 dz(ndim)=dsum2
      endif
*
      return
      end

*****************************************************
******************   MPOLRAD ************************
*****************************************************
      subroutine mpolrad(e1curr,yscurr,xscurr
     .		      ,uncurr,plcurr,pncurr,itagen)
      implicit real*8(a-h,o-z)
      real e1curr,yscurr,xscurr,uncurr,plcurr,pncurr
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire
#include "radinclude/mcpolrad.inc"
#include "radinclude/radgenkeys.inc"
      common/delta5/delta5
      common/tmp/itmp
      dimension tai(5),si(2,3),si0(2,3),tls(2,3,4)

      itmp=1

	 e1=e1curr
	 xs=xscurr
	 ys=yscurr
	 pl=plcurr
	 pn=-pncurr
	 un=uncurr
	 qn=0.

	 call conk2(e1,xs,ys,1)
c
c
c delta is factorizing part of virtual and real leptonic bremsstrahlung
c
      call deltas(delta,delinf,tr)
c      call deltas

      il=1
      in=1

       si(il,in)=0d0
       si0(il,in)=0d0
       tls(il,in,1)=0d0
       tls(il,in,2)=0d0
       tls(il,in,3)=0d0
       tls(il,in,4)=0d0

      isf1=1
      isf2=isf20
      isf3=1

      do ii=1,4
      tai(ii)=0d0
      enddo
      sib=0d0
      sia=0d0

      if(itagen.eq.-1)then
	ita1=1
	ita4=4
      else
	ita1=itagen
	ita4=itagen
      endif

      do 30 ita=ita1,ita4
c
c     sib is born cross section with polarized initial
c     lepton and proton
c     sia is contribution of anomalous magnetic moment.
c
      if(ita.eq.2.and.kill_elas_res.eq.1)goto30
      if(ita.eq.3.and.ire.eq.1)then
	 tai(3)=0d0
c	  write(9,'('' tai = 0.0  '')')
	 goto 30
      end if
      if(ita.eq.1)then
	 call bornin(sib,sia)
      endif
c
c     tai(1),tai(2),tai(3) are contributions of radiative tails:
c		 1 - inelastic
c		 2 - elastic
c		 3 - quasielastic
c
      if(ita.eq.2) call conk2(e1,xs,ys,ita)

	 if(kill_elas_res.ne.2) call qqt(sib,tai(ita))

      if(ita.eq.2) call conk2(e1,xs,ys,1)
  30  continue
c      extai1=exp(alfa/pi*delinf)
c      extai2=((sx-y/tara)**2/s/(s-y/tara))**tr
c      extai3=((sx-y)**2/s/(s-y))**tr

      extai1=1.
      extai2=1.
      extai3=1.
      delinf=0.

c      sig=sib*extai1*(1.+alfa/pi*(delta-delinf))+sia
      sig=sib*redfac*(1.+vertex+vac+small)+sia
c     . 	   +tai(4)+tai(5)
     .		   +tai(4)
     . +tai(1)+(tai(2)*extai2+tai(3)*extai3)/tara

      si(il,in)=si(il,in)+sig
      si0(il,in)=si0(il,in)+sib
      tls(il,in,1)=tls(il,in,1)+tai(1)
      tls(il,in,2)=tls(il,in,2)+tai(2)*extai2/tara
      tls(il,in,3)=tls(il,in,3)+tai(3)*extai3/tara
c      tls(il,in,4)=sib*extai1*(1.+alfa/pi*(delta-delinf))+sia
c     . 	   +tai(4)+tai(5)
      tls(il,in,4)=sib*redfac*(1.+vertex+vac+small)+sia
     .		   +tai(4)
 20   continue
 10   continue

c     write(*,'(1x,8g11.4)')xs,ys,si0(il,in),si(il,in)
c    .	  ,tls(il,in,1),tls(il,in,2),tls(il,in,3)
c    .	  ,tls(il,in,4)
      if(itagen.eq.-1.or.itagen.eq.1)then
c	 if(abs(tine-tls(il,in,1))/tine.gt.5d-2)
c     .  write(*,*)' tine',tine,tls(il,in,1),itagen
	 tine=tls(il,in,1)
      endif
      if(itagen.eq.-1.or.itagen.eq.2)then
c	 if(abs(tnuc-tls(il,in,2))/tnuc.gt.5d-2)
c     .  write(*,*)' tnuc',tnuc,tls(il,in,2),itagen
	 tnuc=tls(il,in,2)
      endif
      if(itagen.eq.-1.or.itagen.eq.3)then
c	 if(abs(tpro-tls(il,in,3))/tpro.gt.5d-2)
c     .  write(*,*)' tpro',tnuc,tls(il,in,3),itagen
	 tpro=tls(il,in,3)
      endif
      if(itagen.eq.-1)then
       tbor=tls(il,in,4)
       sigrad=si(il,in)
       sig1g=si0(il,in)
      endif
c      write(*,'(1x,6g11.4)')sig1g,sigrad,tine,tnuc,tpro,tbor
      sigcor=sigrad/sig1g
       if(kill_elas_res.eq.2)sigcor=1.

c      print *,sngl(s),sngl(xs),sngl(ys)
c      print *,' kill_elas_res=',kill_elas_res,amc2
c      print *,' sig1g=',sig1g,sib
c      print *,' sigrad=',sigrad,sig
c      print *,' tine =',tine
c      print *,' tai(4) =',tai(4)
c      print *,' tnuc =',tnuc
c      print *,' tpro =',tpro
c      print *,' tbor =',tbor
c      print *,' sig1g=',sig1g
c      print *,' sigcor=',sigcor
c      print *,' delta=',alfa/pi*delta
c      print *,' vac  =',vac
c      print *,' vertex=',vertex
c      print *,' small=',small
c      print *,' tai(4)= ',sib*log(redfac)
cc	sigt1=sib*alfa/pi*delta+tai(5)
cc	sigt2=sib*(log(redfac)+vertex+vac+small)
cc	sigt3=sib*alfa/pi*(delta+delta5)
c      print *,' sigt=',sigt1,sigt2,sigt3
1000  end




****************** conk2 **************************************

      subroutine conk2(e1,xs,ys,iittaa)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      if(iittaa.eq.2)then
	 amp=amt
      else
	 amp=amh
      endif
c      print *,amp,amh,iittaa
      call conkin(e1,xs*amh/amp,ys)
      return
      end

****************** conkin *************************************

      subroutine conkin(e1,xs,ys)
c set of kinematical constants
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/pol/as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xss,yss,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
#include "radinclude/mcpolrad.inc"
      ap=2.*amp
      amp2=amp**2
      ap2=2.*amp**2
      s=ap*e1
      x=s*(1-ys)
      sx=s-x
      sxp=s+x
      y=s*xs*ys
      ym=y+al2
      tpl=s**2+x**2
      tmi=s**2-x**2
      w2=amp2+s-y-x
      als=s*s-al2*ap2
      alx=x*x-al2*ap2
      alm=y*y+4.*aml2*y
      aly=sx**2+4.*amp2*y
      sqls=dsqrt(als)
      sqlx=dsqrt(alx)
      sqly=dsqrt(aly)
      sqlm=dsqrt(alm)
      allm=dlog((sqlm+y)/(sqlm-y))/sqlm

      coe=xss/e1/1d3

      axy=pi*(s-x) * coe
      an=2.*alfa**2/sqls*axy*barn*amh/amp
      tamin=(sx-sqly)/ap2
      tamax=(sx+sqly)/ap2

      dcts=(s*sx + ap2*y)/sqly/sqls
      dsts=sin( acos(dcts) )

      as=s/2./aml/sqls
      bs=0.
      cs=-aml/sqls
	 ae=amp/sqls
	 be=0.
	 ce=-s/ap/sqls
      apq=-y*(ae-be)+ce*sx
      apn=(y+4.*aml2)*(ae+be)+ce*sxp
      dk2ks=as*ym+al2*bs+cs*x
      dksp1=as*s+bs*x+cs*ap2
      dapks=2.*(al2*(as*ae+bs*be)+ap2*cs*ce+ym*(as*be+bs*ae)
     .+s*(as*ce+cs*ae)+x*(bs*ce+cs*be))
      return
      end

****************** bornin *************************************

      subroutine bornin(sibor,siamm)
c
c     sibor is born cross section with polarized initial
c     lepton and polarized target
c     siamm is contribution of anomalous magnetic moment.
c
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/pol/as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire
       common/print/ipri1
      dimension sfm0(8),tm(8)
       ipri1=1
      call strf(0d0,0d0,sfm0)
       ipri1=0
      tm(1)=-(2.*aml2-y)
      tm(2)=(-(amp2*y-s*x))/(2.*amp2)
      tm(3)=(2.*(apq*dk2ks-dapks*y)*aml)/amp
      tm(4)=apq/amp*(-(dk2ks*sx-2.*dksp1*y)*aml)/amp2
      tm(7)=(-(4.*aml2+3.*apn**2-3.*apq**2+y))/2.
      tm(8)=apq/amp*(-3.*(apn*sxp-apq*sx))/(2.*ap)
      ek=(3.*apq**2-y)/amp2
      tm(5)=-ek*tm(1)
      tm(6)=-ek*tm(2)
      ssum=0.
      do 1 isf=isf1,isf2,isf3
	ppol=1.
	if(isf.eq.3.or.isf.eq.4)ppol=-pn
	if(isf.ge.5)ppol=qn/6
      ssum=ssum+tm(isf)*sfm0(isf)*ppol
    1 continue
      sibor=ssum*an/y/y*2.
c
c  formula (4) of kukhto and shumeiko paper
c
cc    res1=amp*ww1*(y+4.*aml2)-ww2*(s+x)**2/4./amp
cc    siamm=alfa/pi*al2*allm*(sibor+an*res1/y**2)
      siamm=0.
      return
      end

****************** deltas *************************************

      subroutine deltas(delta,delinf,tr)
*      subroutine deltas
c
c delta is factorizing part of virtual and real leptonic bremsstrahlung
c
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/delta5/delta5
      dimension am2(3)
c
c    am2 : squared masses of charge leptons
c
#include "radinclude/mcpolrad.inc"

      data am2/.26110d-6,.111637d-1,3.18301d0/
      del1=-ym*(alm*allm**2/2.+2.*fspen(2d0*sqlm/(y+sqlm))-pi2/2.)/sqlm
      del2=(3.*y/2.+4.*aml2)*allm-2.

      suml=0.
      do 10 i=1,3
	 a2=2.*am2(i)
	 sqlmi=dsqrt(y*y+2.*a2*y)
	 allmi=dlog((sqlmi+y)/(sqlmi-y))/sqlmi
  10  suml=suml+2.*(y+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./y
      if(y.lt.1.d0)then
	aaa = -1.345d-9
	bbb = -2.302d-3
	ccc = 4.091
      elseif(y.lt.64d0)then
	aaa = -1.512d-3
	bbb =  -2.822d-3
	ccc = 1.218
      else
	aaa = -1.1344d-3
	bbb = -3.0680d-3
	ccc = 9.9992d-1
      endif
      sumh = -(aaa+bbb*log(1.+ccc*y)) *2*pi/alfa
      sum=suml+sumh
c      print *,' vacl_my_old=',alfa/pi*suml
c      print *,' vach_my_old=',alfa/pi*sumh
c      print *,' vac_my_old=',alfa/pi*sum
c      print *,'in-radgen-deltas call vacpol',y
      sum=vacpol(y)
c      print *,' vac_my_new=',alfa/pi*sum

      aj0=2.*(ym*allm-1.)
      deltai=aj0*dlog((w2-amc2)/aml/dsqrt(w2))

      ss=x+y
      xx=s-y
      alss=ss**2-2.*w2*al2
      alxx=xx**2-2.*w2*al2
      sqlss=dsqrt(alss)
      sqlxx=dsqrt(alxx)
      allss=dlog((sqlss+ss)/(-sqlss+ss))/sqlss
      allxx=dlog((sqlxx+xx)/(-sqlxx+xx))/sqlxx
      dlm=dlog(y/aml2)
	 sfpr=dlm**2/2.-dlm*dlog(ss*xx/aml2/w2)
     . -(dlog(ss/xx))**2/2.+fspen((s*x-y*amp2)/ss/xx)-pi2/3.
      delta0=(ss*allss+xx*allxx)/2.+sfpr
      delta=deltai+delta0+del1+del2+sum
      delinf=(dlm-1.)*dlog((w2-amc2)**2/ss/xx)
      tr=alfa/pi*(dlm-1.)

      vac=alfa/pi*sum
      vertex=alfa/pi*del2
      small_old=alfa/pi*(pi2/6.-fspen(1.-amp2*y/s/x)
     .			+ fspen(1.-s/x) + fspen(1.-x/s))
	small=alfa/pi*(-pi2/6.-log(ss*xx/s/x)*log(sx/y)
     .	+log(ss/s)*log(xx/x)
     . -1.5*log(s/x)**2+2.*log(xx/x)*log(s/y)+2.*log(ss/s)*log(x/y)
     . -2.*fspen(-y/x)-2.*fspen(y/s)+2.*fspen(y/sx)+fspen(s*x/ss/xx)
     . -fspen(s*y/ss/sx)-fspen(x*y/sx/xx))
c      print *,' small_old=',small_old
c      print *,' small_new=',small
      redfac=exp(-alfa/pi*(dlm-1.)*log(s*x/(4.*amp2*demin**2) ))
      delta5=-(dlm-1.)*log((w2-amc2)**2*s*x/(4.*amp2*demin**2*ss*xx))
     .	-log(xx/x)*log(amp2*y/s**2)-log(ss/s)*log(amp2*y/x**2)
     .	-2.*fspen(-y/x)-2.*fspen(y/s)+2.*fspen(-tamin)+2.*fspen(-tamax)
     .	-fspen((-y-s*tamin)/xx)
     .	-fspen((-y-s*tamax)/xx)
     .	-fspen(( y-x*tamin)/ss)
     .	-fspen(( y-x*tamax)/ss)
      return
      end

****************** vacpol *************************************

      double precision function vacpol(y)
      implicit real*8(a-h,o-z)
      common/p/pi,pi2,alfa,i1(8),i2(8)
      dimension am2(3),a(5),b(5),c(5)
c
c    am2 : squared masses of charge leptons
c

      data am2/.26110d-6,.111637d-1,3.18301d0/
      data a/0d0,0d0,0d0,1.2227d-3,1.64178d-3/
      data b/2.2877d-3,2.51507d-3,2.79328d-3,2.96694d-3,2.92051d-3/
      data c/4.08041425d0,3.09624477d0,2.07463133d0,1d0,1d0/

      suml=0.
      do 10 i=1,3
	 a2=2.*am2(i)
	 sqlmi=dsqrt(y*y+2.*a2*y)
	 allmi=dlog((sqlmi+y)/(sqlmi-y))/sqlmi
  10  suml=suml+2.*(y+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./y

      if(y .lt. 4d0)then
	k=1
      elseif(y .lt. 16d0)then
	k=2
      elseif(y .lt. 100d0)then
	k=3
      elseif(y .lt. 8317.44d0)then
	k=4
      elseif(y .ge. 8317.44d0)then
	k=5
      else
        print *,'y is',y,k,suml
	stop ' Y<0 in VACPOL'
      endif

      sumh = (a(k)+b(k)*log(1.+c(k)*y)) *2.*pi/alfa
      vacpol=suml+sumh

      end


****************** fspens *************************************

      double precision function fspens(x)
c
c    spence function
c
      implicit real*8(a-h,o-z)
      f=0.d0
      a=1.d0
      an=0.d0
      tch=1.d-16
  1   an=an+1.d0
      a=a*x
      b=a/an**2
      f=f+b
      if(b-tch)2,2,1
  2   fspens=f
      return
      end
****************** fspen **************************************

      double precision function fspen(x)
      implicit real*8(a-h,o-z)
      data f1/1.644934d0/
      if(x)8,1,1
  1   if(x-.5d0)2,2,3
    2 fspen=fspens(x)
      return
    3 if(x-1d0)4,4,5
    4 fspen=f1-dlog(x)*dlog(1d0-x+1d-10)-fspens(1d0-x)
      return
    5 if(x-2d0)6,6,7
    6 fspen=f1-.5*dlog(x)*dlog((x-1d0)**2/x)+fspens(1d0-1d0/x)
      return
    7 fspen=2d0*f1-.5d0*dlog(x)**2-fspens(1d0/x)
      return
    8 if(x+1d0)10,9,9
   9  fspen=-.5d0*dlog(1d0-x)**2-fspens(x/(x-1d0))
      return
  10  fspen=-.5*dlog(1.-x)*dlog(x**2/(1d0-x))-f1+fspens(1d0/(1d0-x))
      return
      end

****************** qqt ****************************************

      subroutine qqt(bo,tai)
      implicit real*8(a-h,o-z)
      external rv2
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire
#include "radinclude/mcpolrad.inc"
      dimension dbtk(8)
      data ep/1d-8/
       dsumtk=0.d0
       derrtk=0.d0
       isumtk=0
       if(ita.eq.1 .or. ita.eq.5)then
	  tade=(w2-amc2)/(ap*demin) -1.d0
	  costkm=(sx-ap2*tade)/sqly
	  if(costkm.ge.1d0)then
	     tai=0.
	     return
	  endif
	  dtkmax=acos( max(-1d0,costkm ))
       else
	 dtkmax=pi
       endif
       call intbtk2(dbtk,nbtk,dtkmax)
       do 10 itk=1,nbtk
c      integrate each bin by ntk subbins
c      ntk=number of bins within the big bin dbtk(i)...dbtk(i+1)
       call inttk2(isumtk,dbtk(itk),dbtk(itk+1),dsumtk,derrtk)
c	write(*,*)itk,isumtk,dsumtk
c
10     continue
      if(ita.le.3)ndxtkm(ita)=isumtk
      tai=dsumtk


      end

****************** tails **************************************

       subroutine tails(ta,tm)
       implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/pol/as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/tmp/itmp
#include "radinclude/mcpolrad.inc"
       common/bseo/ois,oir,oi12,eeis,eeir,eei12,
     . eei1i2,eb,eeb,tm3(6,4,3)
       dimension tm(8,6),ajm2(2),ajm3(3),ii(8)
      data ii/1,2,3,4,1,2,5,6/
      if(iphi.eq.0)then
       b2=(-aly*ta+sxp*sx*ta+2.*sxp*y)/2.
       b1=(-aly*ta-sxp*sx*ta-2.*sxp*y)/2.
       c1=-(4.*(amp2*ta**2-sx*ta-y)*aml2-(s*ta+y)**2)
       c2=-(4.*(amp2*ta**2-sx*ta-y)*aml2-(ta*x-y)**2)
       bb=1./sqly
       sc1=dsqrt(c1)
       sc2=dsqrt(c2)
       bi12=(sxp*(sx*ta+2.*y))/(sc1*sc2*(sc1+sc2))
       bi1pi2=1./sc2+1./sc1
       bis=-b1/sc1/c1+b2/sc2/c2
       bir=b2/sc2/c2+b1/sc1/c1
       b1i=-b1/aly/sqly
       b11i=(3.*b1**2-aly*c1)/2./aly**2/sqly
      else
       sqrtmb=sqrt((ta-tamin)*(tamax-ta)*(s*x*y-y**2*amp2-aml2*aly))
       z1=(y*sxp+ta*(s*sx+ap2*y)-ap*cos(phipoi)*sqrtmb)/aly
       z2=(y*sxp+ta*(x*sx-ap2*y)-ap*cos(phipoi)*sqrtmb)/aly
       bb=1./sqly/pi
       bi12=bb/(z1*z2)
       bi1pi2=bb/z2+bb/z1
       bis=bb/z2**2+bb/z1**2
       bir=bb/z2**2-bb/z1**2
       b1i=bb*z1
       b11i=bb*z1**2
      endif
	 sps=as+bs
	 spe=ae+be
	 ccpe=(ae-be)*ta+2.*ce
	 ccps=(as-bs)*ta+2.*cs
      sis=(2.*bi1pi2*sps+bir*sps*ta+bis*ccps)/2.
      sir=( (2.*bi12*sps*ta+bir*ccps+bis*sps*ta))/2.
      si12=(bi12*ccps+bi1pi2*sps)/2.
      eis=(2.*bi1pi2*spe+bir*spe*ta+bis*ccpe)/2.
      eir=( (2.*bi12*spe*ta+bir*ccpe+bis*spe*ta))/2.
      ei12=(bi12*ccpe+bi1pi2*spe)/2.
      ois=((2.*bi1pi2+bir*ta)*(ccpe*sps+ccps*spe)+(ccpe*ccps+
     . spe*sps*ta**2)*bis+8.*bb*spe*sps+4.*bi12*spe*sps*ta**2)/
     . 4.
      oir=( ((2.*bi12+bis)*(ccpe*sps+ccps*spe)*ta+(ccpe*ccps+
     . spe*sps*ta**2)*bir+4.*bi1pi2*spe*sps*ta))/4.
      oi12=((ccpe*ccps+spe*sps*ta**2)*bi12+(ccpe*sps+ccps*spe)*
     . bi1pi2+4.*bb*spe*sps)/4.
      eeis=((ccpe**2+spe**2*ta**2)*bis+8.*bb*spe**2+4.*bi12*spe
     . **2*ta**2+4.*bi1pi2*ccpe*spe+2.*bir*ccpe*spe*ta)/4.
      eeir=( ((ccpe**2+spe**2*ta**2)*bir+4.*bi12*ccpe*spe*ta+4.
     . *bi1pi2*spe**2*ta+2.*bis*ccpe*spe*ta))/4.
      eei12=((ccpe**2+spe**2*ta**2)*bi12+4.*bb*spe**2+2.*bi1pi2
     . *ccpe*spe)/4.
      ei1pi2=(4.*bb*spe+bi12*spe*ta**2+bi1pi2*ccpe)/2.
      eei1i2=((ccpe**2+spe**2*ta**2)*bi1pi2+4.*(2.*ccpe-spe*ta)
     . *bb*spe+8.*b1i*spe**2+2.*bi12*ccpe*spe*ta**2)/4.
      eb=((ccpe-spe*ta)*bb+2.*b1i*spe)/2.
      eeb=((ccpe-spe*ta)**2*bb+4.*(ccpe-spe*ta)*b1i*spe+4.*b11i
     . *spe**2)/4.
       call ffu(1,bb,bis,bir,bi12,bi1pi2,sir,sis,si12
     .,eis,eir,ei12,ei1pi2,ta)
       call ffu(2,eb,eis,eir,ei12,ei1pi2,oir,ois,oi12
     .,eeis,eeir,eei12,eei1i2,ta)
       call ffu(3,eeb,eeis,eeir,eei12,eei1i2,0d0,0d0,0d0
     .,0d0,0d0,0d0,0d0,ta)
       ajm2(1)=apq/amp
       ajm2(2)=-1./amp
       ajm3(1)=(y-3.*apq**2)/amp2
       ajm3(2)=6.*apq/amp2
       ajm3(3)=-3./amp2
       do 15 i=1,8
       do 13 l=1,6
   13  tm(i,l)=0
       do 10 k=1,i2(i)
       ajk=1.
       if(i.eq.4.or.i.eq.8)ajk=ajm2(k)
       if(i.eq.5.or.i.eq.6)ajk=ajm3(k)
       do 10 j=k,i1(i)+k-1
       tm(i,j)=tm(i,j)+tm3(ii(i),j-k+1,k)*ajk
       if((i.eq.5.or.i.eq.6).and.k.eq.2)
     . tm(i,j)=tm(i,j)+tm3(ii(i),j-k+1,1)*ta/amp2
  10   continue
  15   continue
       return
       end

****************** ffu ****************************************

       subroutine ffu(n,bb,bis,bir,bi12,bi1pi2,sir,sis,si12
     .	      ,eis,eir,ei12,ei1pi2,ta)
       implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/pol/as,bs,cs,ae,be,ce,apn,apq,dk2ks,dksp1,dapks
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/p/pi,pi2,alfa,i1(8),i2(8)
       common/bseo/ois,oir,oi12,eeis,eeir,eei12,
     . eei1i2,eb,eeb,tm3(6,4,3)
      hi2=aml2*bis-ym*bi12
      shi2=aml2*sis-ym*si12
      ehi2=aml2*eis-ym*ei12
      ohi2=aml2*ois-ym*oi12
       goto(10,20,30)n
  10   continue
      tm3(3,1,n)=(8.*(apq*dk2ks-dapks*y)*aml*hi2)/amp
      tm3(3,2,n)=(-2.*((2.*(bi12*dk2ks*ta-2.*shi2)*apq+(2.*shi2-
     . sir*y+sis*ym)*apn+4.*dapks*hi2*ta)-4.*((2.*ei12-eis)*
     . dk2ks-(si12-sis)*apn)*aml2)*aml)/amp
      tm3(3,3,n)=(2.*(((2.*si12+sir-sis)*apn*ta-2.*dk2ks*ei12*ta
     . -6.*ohi2-oir*y+ois*ym)-4.*aml2*oi12)*aml)/amp
      tm3(3,4,n)=(2.*(2.*oi12-oir+ois)*aml*ta)/amp
      tm3(5,1,n)=-2.*(4.*aml2+3.*apn**2-3.*apq**2+y)*hi2
      tm3(5,2,n)=-2.*(6.*aml2*apn*eir-3.*apn**2*bi12*ta+3.*apn*
     . apq*bi1pi2+6.*apq*ehi2+hi2*ta)
      tm3(5,3,n)=-(24.*aml2*eei12-6.*apn*ei1pi2-6.*apq*ei12*ta-
     . 2.*bb-bi12*ta**2)
  20   continue
      tm3(4,1,n)=(-4.*(dk2ks*sx-2.*dksp1*y)*aml*hi2)/amp2
      tm3(4,2,n)=(((2.*(sxp-2.*sx)*shi2+2.*bi12*dk2ks*sx*ta+8.*
     . dksp1*hi2*ta-sir*sxp*y+sis*sxp*ym)-4.*(2.*bi12*dk2ks-bis*
     . dk2ks-si12*sxp+sis*sxp)*aml2)*aml)/amp2
      tm3(4,3,n)=((((sxp*ta-ym)*sis-(sxp*ta-y)*sir+2.*bi12*dk2ks
     . *ta+6.*shi2-2.*si12*sxp*ta)+4.*aml2*si12)*aml)/amp2
      tm3(4,4,n)=(-(2.*si12-sir+sis)*aml*ta)/amp2
      tm3(6,1,n)=(-3.*(apn*sxp-apq*sx)*hi2)/amp
      tm3(6,2,n)=(-3.*(2.*(apn*bir+eir*sxp)*aml2-(2.*bi12*sxp*ta
     . -bi1pi2*sx)*apn+(bi1pi2*sxp+2.*hi2)*apq+2.*ehi2*sx))/(2.*
     . amp)
      tm3(6,3,n)=(-3.*(8.*aml2*ei12-apn*bi1pi2-apq*bi12*ta-ei12*
     . sx*ta-ei1pi2*sxp))/(2.*amp)
  30   continue
      tm3(1,1,n)=-4.*(2.*aml2-y)*hi2
      tm3(1,2,n)=4.*hi2*ta
      tm3(1,3,n)=-2.*(2.*bb+bi12*ta**2)
      tm3(2,1,n)=(((sxp**2-sx**2)-4.*amp2*y)*hi2)/(2.*amp2)
      tm3(2,2,n)=(2.*aml2*bir*sxp-4.*amp2*hi2*ta-bi12*sxp**2*ta+
     . bi1pi2*sxp*sx+2.*hi2*sx)/(2.*amp2)
      tm3(2,3,n)=(2.*(2.*bb+bi12*ta**2)*amp2+4.*aml2*bi12-bi12*
     . sx*ta-bi1pi2*sxp)/(2.*amp2)
       return
       end

****************** podinl *************************************

      double precision function podinl(r)
c
c     integrand (over r )
c
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire
      common/p/pi,pi2,alfa,i1(8),i2(8)
      common/amf2/ta,tm(8,6),sfm0(8)
c      common/tmp/itmp
      dimension sfm(8),iel(8)
      data iel/0,2,1,2,2,4,2,3/
       if(ita.ne.5)call strf(ta ,r,sfm)
      podinl=0.
      do 11 isf=isf1,isf2,isf3
	ppol=1.
	if(isf.eq.3.or.isf.eq.4)ppol=-pn
	if(isf.ge.5)ppol=qn/6
	if(ita.eq.2)ppol=ppol*(amt/amp)**iel(isf)
	 irrend=i1(isf)+i2(isf)-1
	 if(ita.eq.5)irrend=1
      do 1 irr=1,irrend
      if(ita.eq.5)then
	 pres=-sfm0(isf)/2./y**2/r
      else
	 pp=sfm(isf)
	 if(irr.eq.1.and.ita.eq.4)pp=pp-sfm0(isf)*(1.+r*ta/y)**2
	 pres=pp*r**(irr-2)/(y+r*ta)**2/2.
      endif
      podinl=podinl-tm(isf,irr)*pres*ppol
c      print *,' isf=',isf,' irr=',irr,podinl
c      write(*,*)' isf=',isf,' irr=',irr,podinl,sfm(isf)
    1 continue
   11 continue
      podinl=podinl*(1.+alfa/pi*vacpol(Y+R*ta))
c      if(itmp.eq.1)then
c	write(*,*)podinl,tm(1,1),amp,ap,ap2,aml2,aml
c	itmp=0
c      endif
       return
      end

****************** rv2 ****************************************

      double precision function rv2(ta)
c
c     integrand (over ta )
c
      implicit real*8(a-h,o-z)
      external podinl
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire
      common/amf2/taa,tm(8,6),sfm0(8)


      common/p/pi,pi2,alfa,i1(8),i2(8)
#include "radinclude/mcpolrad.inc"


      taa=ta
       call strf(0d0,0d0,sfm0)
       call tails(ta,tm)
      rmin=1d-8
       rcal=ap*demin
      rmax=(w2-amc2)/(1.+ta)
      if(ita.eq.1)then
c	 call dqn32(rmin,rmax,podinl,res)
       dsumtk=0.d0
       derrtk=0.d0
c
	dbmj2=rmax
	dbmj1=rcal
	nmj=nrr
       dr=(rmax-rcal)/nrr
       rcur=rcal-dr*.5d0
       dsum1=0.d0
c      loop over all bins
       do irr=1,nrr
	  rcur=rcur+dr
	  drcurr(irr,itkcur) = rcur
	  ddeler(irr,itkcur) = dr
	  dsum1=dsum1+podinl(rcur)
	  dsigmr(irr,itkcur) = dr*dsum1
       enddo

       res=dr*dsum1

      elseif(ita.eq.2 .or. ita.eq.3)then
	 aa=amt/amp
	 if(ita.eq.3)aa=1.
	 res=podinl((sx-y/aa)/(1d0+ta/aa))/(1.+ta/aa) /aa**2
      elseif(ita.eq.4)then
	  rend=min(rcal,rmax)
	 call dqn32(rmin,rend,podinl,res)
c	 call qunc8(podinl,rmin,rend,1d-4,1d-9,res,er,nn,fl,3500)
      elseif(ita.eq.5)then
	 res=podinl(1.d0/log(rmax/rcal))
      endif
      rv2=res
c      write(9,*)res,dr,rcur
      return
      end



      subroutine strf(ta,rr,sfm)
c
c     the programm calculates deep inelastic (ita=1),
c     elastic (ita=2), quasielastic (ita=3) structure functions
c     in kinematical point (ta,rr).
c	   rr=sx-tt,
c	   ta=(t-y)/rr,
c    where tt=t+amf2-amp2, amf2 is invarint mass of final hadrons
c
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/sxy/s,x,sx,sxp,y,ym,w2,als,alx,alm,aly,
     .sqls,sqlx,sqly,sqlm,allm,an,tamin,tamax,xs,ys,tpl,tmi
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire
       common/print/ipri1
      DIMENSION DHEPAR(3)
      dimension sfm(8)
      DATA DHEPAR /-0.01589,-0.663, -5.96/
      t=y+rr*ta
      if(ita.eq.1 .or. ita.eq.4 .or. ita.eq.5)then

      tt=sx-rr
      amf2=tt-t+amp2
      aks=t/tt
      anu=tt/ap
      epsi=ap2/tt
      g2=0.d0


      if(ire.eq.1)f2=df2h8(t,aks)
      if(ire.eq.2)f2=df2d8(t,aks)
      if(ire.eq.3)f2=(df2h8(t,aks)+2d0*df2d8(t,aks))/3d0


      f1=amp*(1.+anu**2/t)*f2/anu/(1.+r1990f(aks,t))


      if(ire.eq.1)then
	Asym1=aks**0.725
       elseif(ire.eq.2)then
       adeu=8.2885d0
	bdeu=3.23589d-2
	cdeu=.142777d0
	asym1=(ddexp(-adeu*aks)-1d0)*(bdeu**(cdeu)-aks**(cdeu))

       elseif(ire.eq.3)then
	   dz= 1./2.*log(1.+exp(2.0-1000.*aks))
	   dfn=1.
	   df2nf2p=0.67225*(1.0-aks)**1.6254-0.15436*aks**0.048301
     1		  +(.41979+.047331*dz-0.17816*dz**2)
	   df=dfn*(1./((2./df2nf2p)+1))

	 a1nue=0.00024-.00463*(aks**0.1+aks**0.5)
     .		-3.48645*aks+1.59218*aks**1.5
     .		+8.59393*aks**2-5.74029*aks**3
	 asym1=a1nue*df
      endif
      g1=f1*asym1

       b1=0.d0
       b2=0.d0
       b3=0.d0
       b4=0.d0


      goto 10
      endif
c
c   tarn,tarz,tara are n,z,a of the nucleus
c
      epsi=ap2/t
      tarn=tara-tarz
c
c     tf is t in fermi**(-2)
c
      tf=t/chbar**2
c
c   gep,gmp are electric and magnetic form factors of proton
c   s.i.bilenkaya et al pisma zhetf 19(1974)613
c
      call ffpro(t,gep,gmp)
      if(ita.eq.2)then
      tau=t/4./amt**2
      tau1=1.+tau
      if(ire.eq.1)then
	   f1=2.*amp2*tau*gmp**2
	   f2=4.*amp2*tau*(gep**2+tau*gmp**2)/tau1
	   g1=amp2*tau*2.*gmp*(gep+tau*gmp)/tau1
	   g2=2.*amp2*tau**2*gmp*(gep-gmp)/tau1
      elseif(ire.eq.2)then
	  call ffdeu(t,fcdeu,fmdeu,fqdeu)
	   fc=fcdeu*amt
	   fm=fmdeu*amt
	   fq=fqdeu*amt
	   f1=4./3.*tau*tau1*fm**2
	   f2=4./9.*tau*(8.*tau**2*fq**2+6.*tau*fm**2+9.*fc**2)
	   g1=-1./3.*tau*fm*(2.*tau*fq+6.*fc+3.*tau*fm)
	   g2=-1./3.*tau**2*fm*(2.*tau*fq+6*fc-3.*fm)
	   b1=2.*tau**2*fm**2
	   b2=4.*fm**2*tau**2+
     .	      16./3.*tau**3/tau1*(tau*fq+3.*fc-3.*fm)*fq
	   b3=-4./3.*(3.*tau+2.)*fm**2*tau**2+
     .	      16./9.*tau**3/tau1*(tau*fq+3.*fc-3.*fm)*fq
	   b4=4./3.*(6.*tau+1.)*fm**2*tau**2+
     .	      16./9.*tau**3/tau1*(tau*fq+3.*fc+3.*(3.*tau+2.)*fm)*fq
      elseif(ire.eq.3)then
	  call ffhe3(t,ge,gm)
	   f1=2.*amt**2*tau*gm**2
	   f2=4.*amt**2*tau*(ge**2+tau*gm**2)/tau1
	   g1=amt**2*tau*2.*gm*(ge+tau*gm)/tau1
	   g2=2.*amt**2*tau**2*gm*(ge-gm)/tau1
       endif
       else if(ita.eq.3)then
	   tau=t/4./amp**2
	   tau1=1.+tau
	   call ffquas(t,geun,gmun,gepo,gmpo)
	   f1=2.*amp2*tau*gmun**2
	   f2=4.*amp2*tau*(geun**2+tau*gmun**2)/tau1
	   g1=amp2*tau*2.*gmpo*(gepo+tau*gmpo)/tau1
	   g2=2.*amp2*tau**2*gmpo*(gepo-gmpo)/tau1
      endif
  10  continue

	   b1=0.
	   b2=0.
	   b3=0.
	   b4=0.
      sfm(1)=f1+qn/6.*b1
      sfm(2)=epsi*(f2+qn/6.*b2)
      sfm(3)=epsi*(g1+g2)
      sfm(4)=epsi**2*g2
      sfm(5)=epsi**2*b1
      sfm(6)=epsi**3*(b2/3.+b3+b4)
      sfm(7)=epsi*(b2/3.-b3)
      sfm(8)=epsi**2*(b2/3.-b4)
      return
      end

********************** f1sfun ***********************************
      double precision function f1sfun(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      f2=f2sfun(aks,t)
      anu=t/ap/aks
      f1sfun=amp*(1.+anu**2/t)*f2/anu/(1.+r1990f(aks,t))
      end

CDECK  ID>, F2SFUN.
********************** f2sfun ***********************************
      double precision function f2sfun(aks,t)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
#include "radinclude/tailcom.inc"
      if(ire.eq.1)f2sfun=df2h8(t,aks)
      if(ire.eq.2)f2sfun=df2d8(t,aks)
      if(ire.eq.3)f2sfun=(df2h8(t,aks)+2d0*df2d8(t,aks))/3d0
      end


CDECK  ID>, R1990.
********************** r1990 ************************************
      double precision function r1990f(dx,dq2)
      implicit none
#include "options.inc"
      double precision dq2,dx,r1990fbosted,r1990fradgen
c
      if(cl_radsfs.eq.1) then ! standard radgen
       r1990f=r1990fradgen(dx,dq2)
      else if(cl_radsfs.eq.2) then ! Bosted SFs
       r1990f=r1990fbosted(dx,dq2)
      else
       print *,'ERR: Unknown input for RADGEN' 
       stop
      endif
      return
      end
c
      double precision function r1990fbosted(dx,dq2)
      implicit none
#include "options.inc"
            double precision dq2,dx
            double precision Z,A,qsq,wsq,F1,F2,rc,g2
********************** r1990 ************************************
            if(cl_target.eq.1) then ! hydrogen
             Z=1.0
             A=1.0
           else if(cl_target.eq.2) then ! deuteron 
             Z=1.0
             A=2.0
           else if(cl_target.eq.3) then ! helium 
             Z=2.0
             A=3.0
           else 
             print *,'FIX Targets!!!'
             stop
           endif
            qsq=dq2
            Wsq=dq2*(1.D0-dx)/dx+0.88035
      call F1F2IN09(Z, A, QSQ, Wsq, F1, F2, rc) 
            g2=4*0.88035*dx*dx/dq2
       r1990fbosted=(1.D0+g2)*F2/2.0/dx/F1-1.D0
      return
      end 
c
      double precision function r1990fradgen(aks,tc)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      t=tc
      if(tc.lt..35d0)t=0.35
      teta=1.+12.*t/(1.+t)*(.125**2/(aks**2+.125**2))
      zn=teta/log(t/.04)
      ra=.0672*zn+.4671/(t**4+1.8979**4)**(.25d0)
      rb=.0635*zn+.5747/t-.3534/(t**2+.09)
      rc=.0599*zn+.5088/sqrt((t-5.*(1.-aks)**5)**2+2.1081**2)
      rrr=(ra+rb+rc)/3.
      r1990fradgen=rrr
c	call r1990(tc,aks,r1990f)
      return
      end


CDECK  ID>, G1SFUN.
********************** g1sfun ***********************************
      double precision function g1sfun(aks,t)
      implicit real*8(a-h,o-z)
#include "radinclude/tailcom.inc"

      if(ire.eq.1)then
c	a1p=1.90202d-2+aks**(-1.16312d-3)*(1d0-ddexp(-1.84517*aks))
c	a1p=1.1092*aks**(-.13658)*(1d0-ddexp(-1.3258*aks))
	a1p=aks**0.7
	g1sfun=a1p*f1sfun(aks,t)

       elseif(ire.eq.2)then
       adeu=8.2885d0
	bdeu=3.23589d-2
	cdeu=.142777d0
	a1d=(ddexp(-adeu*aks)-1d0)*(bdeu**(cdeu)-aks**(cdeu))
c	 a1d=0.15225*aks**(.92104)*(1d0-ddexp(-15.272*aks))
	g1sfun=a1d*f1sfun(aks,t)

       elseif(ire.eq.3)then
c	 a1he3=-0.01589*aks**(-.663)*(1d0-ddexp(-5.96*aks))
	   dz= 1./2.*log(1.+exp(2.0-1000.*aks))
	   dfn=1.
	   df2nf2p=0.67225*(1.0-aks)**1.6254-0.15436*aks**0.048301
     1		  +(.41979+.047331*dz-0.17816*dz**2)
	   df=dfn*(1./((2./df2nf2p)+1))

	 a1nue=0.00024-.00463*(aks**0.1+aks**0.5)
     .		-3.48645*aks+1.59218*aks**1.5
     .		+8.59393*aks**2-5.74029*aks**3
	 a1he3=a1nue*df
	g1sfun=a1he3*f1sfun(aks,t)
      endif

      end


CDECK  ID>, G2SFUN.
********************** g2sfun ***********************************
      double precision function g2sfun(aks,t)
      implicit real*8(a-h,o-z)
      g2sfun=0.
      end




      double precision function df2d8(dq2,dx)
      implicit none
#include "options.inc"
      double precision dq2,dx,df2d8bosted,df2d8radgen
c
      if(cl_radsfs.eq.1) then ! standard radgen
       df2d8=df2d8radgen(dq2,dx)
      else if(cl_radsfs.eq.2) then ! Bosted SFs
       df2d8=df2d8bosted(dq2,dx)
      else
       print *,'ERR: Unknown input for RADGEN' 
       stop
      endif
      return
      end

      double precision function df2d8bosted(dq2,dx)
      implicit none
*:=====================================================================:
*
            double precision dq2,dx
            double precision Z,A,qsq,wsq,F1,F2,rc
            Z=1.0
            A=2.0
            qsq=dq2
            Wsq=dq2*(1.D0-dx)/dx+0.88035
      call F1F2IN09(Z, A, QSQ, Wsq, F1, F2, rc) 
        print *,'new-df2',Wsq,qsq,dx, F1, F2, rc
        df2d8bosted=F2
        return
        end


      double precision function df2d8radgen(dq2,dx)
*:=====================================================================:
*:								       :
*:	author:    m.dueren	   last update: 06.03.1991	       :
*:				   tested: yes			       :
*:								       :
*:	arguments: dq2,dx: double prec. input xbj,q2		       :
*:		   df2d8* double prec f2  output		       :
*:								       :
*:	called by: mkf2 					       :
*:								       :
*:	action:    calculate f2 structure function of the deuteron     :
*:		   nmc fit of dis-region with 8 parameters	       :
*:		   data of nmc (1992), slac,bcdms		       :
*:								       :
*:		   parametrized with a1,bi (i=1...4) as 	       :
*:								       :
*:		   f2_dis(x,q2) ~prop.				       :
*:		     [1/b(n1,n2)*x**n1*(1-x)**n2 + 1/3*n3*(1-x)**n4 ]  :
*:		     *s(x,q2)					       :
*:			  with	x = (q2+m_a)/(2m*nu + m_b**2)	       :
*:				ni= ai+bi*s			       :
*:				s = ln(ln(q2+m_a**2)/lambda)/ln(..q2_0):
*:		   reference:					       :
*:		   the new muon collaboration			       :
*:		   nuclear physics b 371 (1992) 3-31		       :
*:=====================================================================:
c
      implicit double precision (d)
c
c
c *** d1,...,d8 = 8 param of nmc, slac, bcdms (92)
c *** d9,...,d10 = 2 parameters: (1 for resonance) + (1 for background)
c *** daw,dbw =  weizmann variable in bodek's d2 fit
c	     values: daw=1.512(gev2), dbw=0.351(gev2)
c	     ref:  bodek et al., p.r.d20(1979)1427.
c	     see p.1495, eq(5.1) and table viii
c
c *** dl2 = lamda**2 = 0.2**2 = 0.04 (gev2)
c *** q0**2 = 2 gev2 ... (2+0.351)/0.04 = 58.771
c *** fit by y.m.(25-nov-88 19h43m14s)
c
      data d1,d2,d3,d4,d5,d6,d7,d8
     :	   ,d9,d10
     :	   ,daw,dbw
c
c     f2 from nmc, slac, bcdms - data (92)
     :	   /.764053,-.171774,3.48979,.611064,.946086
     :	   ,1.06818,13.8472,-2.40967
c     resonance-region
     :	   ,.89456,.16452
     :	   ,1.512,  .351 /
c
c
      df2d8radgen=1.d-30
      dw2 = .8803686078d0+dq2*(1.d0/dx-1.d0)
      dwmas = dsqrt(dw2)
      ddw = (dwmas-1.03d0)
c
c *** ddw = w - (resonance threshold - smearing 0f 0.05 gev)
c *** lamda(qcd) = 0.2 gev
c
      if(ddw.le.0.d0) return
c
c *** use weizmann variable for low q2
c *** values: daweiz=1.512(gev2), dbweiz=0.351(gev2)
c *** ref:    bodek et al., p.r.d20(1979)1427.
c ***	      see p.1495, eq(5.1) and table viii
c
      dq2w=dq2+dbw
      dxw=dq2w/(dq2/dx+daw)
c
      dsbar = dlog(dlog(dq2w/.04d0)) - 1.404555751d0
c
      detav1 = d1+d2*dsbar
      detav2 = d3+d4*dsbar
      detas1 = d5+d6*dsbar
      detas2 = d7+d8*dsbar
c
      dxw1=1.d0-dxw
      de1=detav1
      de2=detav2+1.d0
c
c *** supression bue to "quark sum rule"
c ***	    sup.fact.= 1 - gd**2, with gd=nucleon dipole f.f.
c *** further supression due to low w phase space volume
c
      den1=(1.d0+dq2/.71d0)
      dgd2=1.d0/den1**4
      dssum = (1.d0-dgd2)
      ds = dssum * (1.d0-dexp(-4.177d0*ddw))
c
      df2d8radgen =
     :	   ( .83333333333333*dgamma(de1+de2)/dgamma(de1)/dgamma(de2)
     :	   *dxw**detav1*dxw1**detav2
     :	   + .33333333333333*detas1*dxw1**detas2  ) * ds
c
c *** resonance contribution
c
      dres = 0.d0
c      dres2=0.d0
c      dres3=0.d0
c
c *** lorentzian resonance ( small fermi-smearing effect)
c ***			     gamma(fermi)=0.0883 gev
c *** gamma(d2) = sqrt ( gamma(h2)**2 + gamma(fermi)**2 )
c 1.232**2 = 1.518
c 1.520**2 = 2.310
c 1.681**2 = 2.826
c 1.232**2 * 0.15**2 = 0.0342
c 1.520**2 * 0.14**2 = 0.0453
c 1.681**2 * 0.14**2 = 0.0554
c
      if (dwmas .le. 2.3d0) then
	 dres1 = d9**2*dexp(-(dwmas-1.232d0)**2/.0053d0)
     :	      /den1**3
c      dres2 = d10**2/( (dw2-2.310d0)**2 + 0.0453d0 )
c    :	      * dgd2
c      dres3 = d11**2/( (dw2-2.826d0)**2 + 0.0554d0 )
c    :	      * dgd2
c      endif
c
c *** background under resonances
c
c    mp**2 = 0.8803686078
c    mp**2-m(pi)**2=0.8608892416
c *** dqs = momentum of one pion, decaying from the resonance, in cm
c frame
c
	 dw2m = (dwmas+.05d0)**2
	 dqs=dsqrt((dw2m+0.8608892416d0)**2/4.d0/dw2m-0.8803686078d0)
	 dbg = (d10**2*dqs )
     :	      * dexp(-0.5d0*ddw**2) / den1
	 dres=(dres1 + dbg) * dssum
      endif
c
c *** total f2 of d2
c
      df2d8radgen = df2d8radgen + dres
c
      if (df2d8radgen.gt.0d0) return
      df2d8radgen = 1.d-30
      return
      end

      double precision function df2h8(dq2,dx)
      implicit none
#include "options.inc"
      double precision dq2,dx,df2h8radgen,df2h8bosted
c
      if(cl_radsfs.eq.1) then ! standard radgen
       df2h8=df2h8radgen(dq2,dx)
      else if(cl_radsfs.eq.2) then ! Bosted SFs
       df2h8=df2h8bosted(dq2,dx)
      else
       print *,'ERR: Unknown input for RADGEN' 
       stop
      endif
      return
      end
c
      double precision function df2h8bosted(dq2,dx)
      implicit none
#include "options.inc"
*:=====================================================================:
            double precision dq2,dx
            double precision Z,A,qsq,wsq,F1,F2,rc,g2
            if(cl_target.eq.1) then ! hydrogen
             Z=1.0
             A=1.0
           else if(cl_target.eq.2) then ! deuteron 
             Z=1.0
             A=2.0
           else if(cl_target.eq.3) then ! helium 
             Z=2.0
             A=3.0
           else 
             print *,'FIX Targets!!!'
             stop
           endif
            qsq=dq2
            Wsq=dq2*(1.D0-dx)/dx+0.88035
      call F1F2IN09(Z, A, QSQ, Wsq, F1, F2, rc) 
c        print *,'new-hf2',Wsq,qsq,dx, F1, F2, rc
            df2h8bosted=F2
c      print *,'xxq2xF2',dx,dq2,F2
            return
            end
c
      double precision function df2h8radgen(dq2,dx)
*:=====================================================================:
*:								       :
*:	author:    m.dueren	   last update: 06.03.1991	       :
*:				   tested: yes			       :
*:								       :
*:	arguments: dq2,dx: double prec. input xbj,q2		       :
*:		   df2h8* double prec f2  output		       :
*:								       :
*:	called by: mkf2 					       :
*:								       :
*:	action:    calculate f2 structure function of the proton       :
*:		   nmc fit of dis-region with 8 parameters	       :
*:		   data of nmc (1992)slac,bcdms 		       :
*:								       :
*:		   parametrized with a1,bi (i=1...4) as 	       :
*:								       :
*:		   f2_dis(x,q2) ~prop.				       :
*:		     [1/b(n1,n2)*x**n1*(1-x)**n2 + 1/3*n3*(1-x)**n4 ]  :
*:		     *s(x,q2)					       :
*:			  with	x = (q2+m_a)/(2m*nu + m_b**2)	       :
*:				ni= ai+bi*s			       :
*:				s = ln(ln(q2+m_a**2)/lambda)/ln(..q2_0):
*:		   reference:					       :
*:		   the new muon collaboration			       :
*:		   nuclear physics b 371 (1992) 3-31		       :
*:=====================================================================:
c
      implicit double precision (d)
*
c *** d1,...,d8 = 8 param of nmc, slac, bcdms (92)
c *** d9,...,d10 = 2 parameters: (1 for resonance) + (1 for background)
c *** daw,dbw =  weizmann variable in bodek's d2 fit
c	     values: daw=1.512(gev2), dbw=0.351(gev2)
c	     ref:  bodek et al., p.r.d20(1979)1427.
c	     see p.1495, eq(5.1) and table viii
c
c *** dl2 = lamda**2 = 0.2**2 = 0.04 (gev2)
c *** q0**2 = 2 gev2 ... (2+0.351)/0.04 = 58.771
c *** fit by y.m.(25-nov-88 19h43m14s)
*
      data d1,d2,d3,d4,d5,d6,d7,d8
     :	   ,d9,d10,d11,d12,d13,d14
     :	   ,d15,d16
     :	   ,daw,dbw
c     f2 from nmc, slac, bcdms	data '92 (final)
     :	   /.886627,-.11191,3.3951,1.04064,1.02702,1.40335,12.4577,-
     :	   .100622
c     resonance-region:
     :	   ,.1179, .044735, .038445, .27921, 8.8228d-5, 6.2099d-5
     :	   ,1.421,1.2582
     :	   ,1.642, .376/
c
      df2h8radgen =1.d-30
      dw2 = .8803686078d0+dq2*(1.d0/dx-1.d0)
      dwmas = dsqrt(dw2)
      ddw = (dwmas-1.08d0)
c
c *** ddw = w - (resonance threshold)
c *** lamda(qcd) = d2: 0.20 gev
c *** lamda(qcd) = h2: 0.15 gev
c
      if(ddw.le.0.d0) return
c
c *** use weizmann variable for low q2
c *** values = d2 : daweiz=1.512(gev2), dbweiz=0.351(gev2)
c *** values = h2 : daweiz=1.642(gev2), dbweiz=0.376(gev2)
c *** ref:    bodek et al., p.r.d20(1979)1427.
c ***	      see p.1495, eq(5.1) and table viii
c
      dq2w=dq2+dbw
      dxw=dq2w/(dq2/dx+daw)
      dsbar=dlog(dlog(dq2w/.0225d0)) - 1.538942135d0
c      dsbar=dlog(dlog(dq2w/.0225d0)/dlog((2.d0+dbw)/.0225d0))
c      dsbar = dlog(dlog(dq2w/.04d0)) - 1.404555751d0
      detav1 = d1+d2*dsbar
      detav2 = d3+d4*dsbar
      detas1 = d5+d6*dsbar
      detas2 = d7+d8*dsbar
c
      dxw1=1.d0-dxw
      de1=detav1
      de2=detav2+1.d0
c
c *** supression due to "quark sum rule"
c ***	    sup.fact.= 1 - gd**2, with gd=nucleon dipole f.f.
c *** further supression due to low w phase space volume
c
      den1=(1.d0+dq2/.71d0)
      dgd2=1.d0/den1**4
      dssum = (1.d0-dgd2)
      dsthr = 1.d0
      if( ddw .le. 5.0d0) then
	 dtemp = dexp(ddw*3.98213222d0)
	 dsthr = (dtemp-1.d0)/(dtemp+1.d0)
      endif
      ds = dssum * dsthr
c
      df2h8radgen =
     :	   ( .83333333333333*dgamma(de1+de2)/dgamma(de1)/dgamma(de2)
     :	   *dxw**detav1*dxw1**detav2
     :	   + .33333333333333*detas1*dxw1**detas2  ) * ds
c
c *** resonance region
c
      dres = 0.d0

      if(ddw .le. 5.0d0) then
c
c *** >>> + background under the resonance
c
c *** appropriate kinematic variables
c     ... dqs = momentum of one pion in pi-nucleon c.m.-system
c		 in the case of single pion production
c *** n.b.  dqs = 0  at  w = 1.08gev (inelastic threshold)
c	 mp**2 = 0.8803686078
c	 mp**2-m(pi)**2=0.8608892416
c
	 dqs2 = (dw2+0.8608892416d0)**2/4.d0/dw2 - 0.8803686078d0
	 dqs = dsqrt(dqs2)
c
c *** >>> + resonance shape
c
c *** lorentzian resonance
c *** this includes the w**2-dependence of the res.width.
c
c *** appropriate kinematic variables
c  1) correction to res.width due to resonance threshold
c     ... dqs = momentum of one pion in pi-nucleon c.m.-system
c		 in the case of single pion production
c  2) correction to res.width due to the q2-dependence
c     ... dks = momentum of virtual photon in pi-n c.m.-system
c
	 if(ddw .le. 1.d0) then
c
	    dks2 =
     :		 (dw2+0.8803686078d0+dq2)**2/4.d0/dw2 - 0.8803686078d0
	    dks = dsqrt(dks2)
c
c *** resonance form factor (effective guess!)
c
	    dresff = 1. / den1**(d15**2)
c
c *** 1236
c     wres**2	      = 1.232**2	    = 1.518
c     (wres*gamma)**2 = 1.232**2 * 0.119**2 = 0.02149
c
	    dw2r = 1.518d0
	    dqsr2 =
     :		 (dw2r+0.8608892416d0)**2/4.d0/dw2r - 0.8803686078d0
	    dqsr = dsqrt(dqsr2)
	    dksr2 =
     :		 (dw2r+0.8803686078+dq2)**2/4./dw2r - 0.8803686078
	    ddcorr = (dqs/dqsr) * (1.+.16/dqsr2)/(1.+.16/dqs2)
	    dncorr = ddcorr * (dksr2+.16)/(dks2+.16)
	    ddcorr = ddcorr**2
	    dres1 = d9**2 * dncorr
     :		 /( (dw2-1.518d0)**2 + 0.02149d0*ddcorr )
c
c *** 1520
c     wres**2	      = 1.520**2	    = 2.310
c     (wres*gamma)**2 = 1.520**2 * 0.097**2 = 0.02127
c     n.b. q2-dependence of the width is neglected
c
	    dw2r = 2.310d0
	    dqsr2 =
     :		 (dw2r+0.8608892416d0)**2/4.d0/dw2r - 0.8803686078d0
	    ddcorr = dqs2/dqsr2
	    dres2 = d10**2 * ddcorr
     :		 / ( (dw2-2.310d0)**2 + 0.02127d0*ddcorr )
c
c *** 1681
c     wres**2	      = 1.681**2	    = 2.826
c     (wres*gamma)**2 = 1.681**2 * 0.105**2 = 0.03115
c     n.b. q2-dependence of the width is neglected
c
	    dw2r = 2.826d0
	    dqsr2 =
     :		 (dw2r+0.8608892416d0)**2/4.d0/dw2r - 0.8803686078d0
	    dqsr = dsqrt(dqsr2)
	    ddcorr = ( dqs/dqsr )**3
	    dres3 = d11**2 * ddcorr
     :		 / ( (dw2-2.826d0)**2 + 0.03115d0*ddcorr )
c
c
c *** sum of all resonances
c	  * resonance form factor(q2-dependence)
c
	    dres = (dres1 +dres2 +dres3)*dresff
c
c *** end of resonance calculation (only if   ddw < 1.0 gev)
c
	 endif
c
c *** background under the resonances
c     n.b. exp(-0.92**2*3.5393817) = 0.05
c	 dbgff = dq2val/dxval /den1**(dp(16)**2)
c
	 dbgff = 1. / den1**(d16**2)
	 dbg = (d12**2*dsqrt(dqs) +d13**2*dqs +d14**2*dqs2 )
     :	      * dbgff * dexp(-0.5d0*ddw**2)
c
c *** (resonance region) = ( (resonance) + (background) )
c			   * dssum(=suppression)
c
	 dres = (dres + dbg) * dssum
c
c *** end of resonance calculation
c
      endif
c
c *** (total) = (qcd part) + (resonance region)
c
      df2h8radgen = df2h8radgen + dres
c
      if(df2h8radgen .gt. 0.d0) return
c
      df2h8radgen=1.d-30
      return
      end

      double precision function gammln(xx)
      implicit real*8(a-h,o-z)
      dimension cof(6)
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     *	 -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 11 j=1,6
       x=x+one
       ser=ser+cof(j)/x
11    continue
      gammln=tmp+log(stp*ser)
      return
      end
      double precision function dgamma(xx)
      implicit real*8(a-h,o-z)
	dgamma=exp(gammln(xx))
      end

****************** ffpro **************************************
      subroutine ffpro(t,gep,gmp)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      gep=1.2742/(1.+t/0.6394**2)-.2742/(1.+t/1.582**2)
      gmp=(1.3262/(1.+t/0.6397**2)-.3262/(1.+t/1.3137**2))*amm
c     gep=1./((1.+.61*t)*(1.+2.31*t)*(1.+.04*t))
c     gmp=amm*gep
      end


      subroutine ffdeu(t,gc,gm,gq)
      implicit real*8(a-h,o-z)
      parameter (c2i3 = 6.6666666666666667d-1)	! 2/3
      dimension a(4),b(4),c(4)
      dimension al2ar(4),be2ar(4),ga2ar(4)
      data amd/1.8756280D0/chbar/0.19732858d0/
      data amp/0.9382796D0/
      data dmu/0.857406d0/dqu/25.84d0/
      gd2=1d0/(1.d0+t/4./0.8952**2)**4
      eta=t/4.d0/amd**2
      gd2e=gd2/(2d0*eta+1d0)
      sq2e=sqrt(2d0*eta)

      al2ar(1)=1.8591*chbar**2
      al2ar(4)=2d0*amd*0.58327d0
      be2ar(1)=19.586*chbar**2
      be2ar(4)=2d0*amd*0.1d0
      ga2ar(1)=1.0203*chbar**2
      ga2ar(4)=2d0*amd*0.17338d0
      do i=2,3
	 al2ar(i)=al2ar(1)+(al2ar(4)-al2ar(1))/3d0*(i-1)
	 be2ar(i)=be2ar(1)+(be2ar(4)-be2ar(1))/3d0*(i-1)
	 ga2ar(i)=ga2ar(1)+(ga2ar(4)-ga2ar(1))/3d0*(i-1)
      enddo

      a(1)=2.4818*chbar**2
      a(2)=-10.850*chbar**2
      a(3)=6.4416*chbar**2
      a(4)=al2ar(4)*(1d0-a(2)/al2ar(2)-a(3)/al2ar(3)-a(1)/al2ar(1))
      b(1)=-1.7654*chbar
      b(2)=6.7874*chbar
      bzn=1d0/be2ar(4)-1d0/be2ar(3)
      bbb=(2d0-dmu*amd/amp)/2d0/sqrt(2d0)/amd
      b(3)=(b(1)/be2ar(1)+b(2)/be2ar(2)-b(1)/be2ar(4)-b(2)/be2ar(4)
     . -bbb)/bzn
      b(4)=-(b(1)/be2ar(1)+b(2)/be2ar(2)-b(1)/be2ar(3)-b(2)/be2ar(3)
     . -bbb)/bzn
      ccc=(1d0-dmu*amd/amp-dqu)/4./amd**2
      c(1)=-0.053830d0
      znc2=ga2ar(1)*(ga2ar(3)*ga2ar(4)-ga2ar(2)*ga2ar(3)
     .	+ga2ar(2)**2-ga2ar(2)*ga2ar(4))
      c(2)=-ga2ar(2)/znc2*(c(1)*(
     .ga2ar(3)*ga2ar(4)-ga2ar(1)*ga2ar(3)+ga2ar(1)**2-ga2ar(1)*ga2ar(4)
     .)-ccc*ga2ar(3)*ga2ar(4)*ga2ar(1) )
      znc3=ga2ar(1)*(ga2ar(3)-ga2ar(2))*(ga2ar(4)-ga2ar(3))
      c(3)=ga2ar(3)/znc3*(c(1)*(
     .ga2ar(2)*ga2ar(4)-ga2ar(1)*ga2ar(4)+ga2ar(1)**2-ga2ar(1)*ga2ar(2)
     .)-ccc*ga2ar(2)*ga2ar(4)*ga2ar(1) )
      znc4=ga2ar(1)*(ga2ar(4)-ga2ar(2))*(ga2ar(4)-ga2ar(3))
      c(4)=-ga2ar(4)/znc4*(c(1)*(
     .ga2ar(2)*ga2ar(3)-ga2ar(1)*ga2ar(3)+ga2ar(1)**2-ga2ar(1)*ga2ar(2)
     .)-ccc*ga2ar(2)*ga2ar(3)*ga2ar(1) )

      g00=0d0
      gp0=0d0
      gpm=0d0
      sqt=sqrt(t)
      do i=1,4
	g00=g00+a(i)/(al2ar(i)+t)
	gp0=gp0+sqt*b(i)/(be2ar(i)+t)
	gpm=gpm+t*c(i)/(ga2ar(i)+t)
      enddo

      gc=gd2e*( (1d0-c2i3*eta)*g00+4d0*c2i3*sq2e*gp0
     .			       +c2i3*(2d0*eta-1d0)*gpm)
      gm=gd2e*(2d0*g00+2d0*(2d0*eta-1d0)/sq2e*gp0-2d0*gpm)
      gq=gd2e*(-g00+2d0/sq2e*gp0-(1d0+1d0/eta)*gpm)

      end

****************** ffhe3 **************************************
      subroutine ffhe3(t,ge,gm)
c
c   l.i.shiff phys. rev. 133(1964)3b,802
c
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      tf=t/chbar**2
       qf=sqrt(tf)
      a=.675
      b=.366
      c=.836
      am=.654
      bm=.456
      cm=.821
      d=-6.78d-3
      p=.9
      q0=3.98
      f0=ddexp(-a**2*qf**2) - b**2*qf**2*ddexp(-c**2*qf**2)
      fm=ddexp(-am**2*qf**2) - bm**2*qf**2*ddexp(-cm**2*qf**2)
      df=d*ddexp(-((qf-q0)/p)**2)
      ge=f0+df
      gm=fm*tara/tarz * (-2.13)
      end

****************** ffquas **************************************

      subroutine ffquas(t,geun,gmun,gepo,gmpo)
      implicit real*8(a-h,o-z)
      common/cmp/amp,amp2,ap,ap2,aml,aml2,al2,amc2,amh,
     .amt,tara,tarz,fermom,amm,amn,chbar,barn,isf20
      common/tail/un,pl,pn,qn,ita,isf1,isf2,isf3,ire
      call ffpro(t,gep,gmp)
      tf=t/chbar**2
      tau=t/4./amp**2
      tau1=1.+tau
c
c   t. de forest and d.j. valecka adv. in phys. 15(1966) no.57
c
      if(ire.eq.2)then
             supele=supst(t)
c            qbold=sqrt(t*tau1)
c            supele=1.-dsbern(qbold)**2
             supmag=supele
      else       
	supele=1.
	supmag=1.
	     if(t.le.(2.d0*fermom)**2)then
	       sqrat=dsqrt(t)/fermom
	       supele=0.75*sqrat-sqrat**3/16.
	       supmag=supele
	     endif
        endif   
	geun=gep*dsqrt(supele*tarz)
	tarn=tara-tarz
	gmun=gep*dsqrt(supmag*(tarz*amm**2+tarn*amn**2))
      if(ire.eq.2)then
        gepo=geun
        gmpo=gmun
      else
	gepo=0.
	tarn=tara-tarz
	gmpo=gep*dsqrt(supmag*( 	  tarn*amn**2))
        endif
	end

********************** supst ************************************ 

      double precision function supst(t)
      implicit real*8(a-h,o-z)
      data chbar/.197328d0/
c
c     tf is t in fermi**(-2)
c
      tf=t/chbar**2
c
c   s.stein et al phys. rev. 12(1975)1884 (appendix 1)
c
               sqtf=dsqrt(tf)
               delff=(datan(sqtf/.93d0)-2.*datan(sqtf/3.19d0)+
     .    datan(sqtf/5.45d0))*1.580d0/sqtf
               delff=dmax1(0.d0,delff)
               supele=1.-delff**2
               supst=dmax1(0.d0,supele)
       return
      end


****************** ddexp **************************************

      double precision function ddexp(x)
      implicit real*8(a-h,o-z)
	ddexp=0.
	if(x.gt.-50.)ddexp=exp(x)
      return
      end


      subroutine dqn32(xl,xu,fct,y)
c
c
      external fct
      double precision xl,xu,y,a,b,c,fct
c
      a=.5d0*(xu+xl)
      b=xu-xl
      c=.49863193092474078d0*b
      y=.35093050047350483d-2*(fct(a+c)+fct(a-c))
      c=.49280575577263417d0*b
      y=y+.8137197365452835d-2*(fct(a+c)+fct(a-c))
      c=.48238112779375322d0*b
      y=y+.12696032654631030d-1*(fct(a+c)+fct(a-c))
      c=.46745303796886984d0*b
      y=y+.17136931456510717d-1*(fct(a+c)+fct(a-c))
      c=.44816057788302606d0*b
      y=y+.21417949011113340d-1*(fct(a+c)+fct(a-c))
      c=.42468380686628499d0*b
      y=y+.25499029631188088d-1*(fct(a+c)+fct(a-c))
      c=.39724189798397120d0*b
      y=y+.29342046739267774d-1*(fct(a+c)+fct(a-c))
      c=.36609105937014484d0*b
      y=y+.32911111388180923d-1*(fct(a+c)+fct(a-c))
      c=.33152213346510760d0*b
      y=y+.36172897054424253d-1*(fct(a+c)+fct(a-c))
      c=.29385787862038116d0*b
      y=y+.39096947893535153d-1*(fct(a+c)+fct(a-c))
      c=.25344995446611470d0*b
      y=y+.41655962113473378d-1*(fct(a+c)+fct(a-c))
      c=.21067563806531767d0*b
      y=y+.43826046502201906d-1*(fct(a+c)+fct(a-c))
      c=.16593430114106382d0*b
      y=y+.45586939347881942d-1*(fct(a+c)+fct(a-c))
      c=.11964368112606854d0*b
      y=y+.46922199540402283d-1*(fct(a+c)+fct(a-c))
      c=.7223598079139825d-1*b
      y=y+.47819360039637430d-1*(fct(a+c)+fct(a-c))
      c=.24153832843869158d-1*b
      y=b*(y+.48270044257363900d-1*(fct(a+c)+fct(a-c)))
      return
      end

       subroutine intbtk2(dbtk,nbtk,dtkmax)
      implicit real*8 (a-h,o-z)
c
c ***  choose 7 bins for integration over the theta k angle
       dimension  dbtk(8)
#include "radinclude/sxycom.inc"
#include "radinclude/cmpcom.inc"
#include "radinclude/ppicom.inc"
#include "radinclude/mcpolrad.inc"

       data dc /3.5d0/
c
c      bins are chosen according to the s peak and p peak
c      width of the peaks:

       dwids=dsqrt(ap*aml/s)/dc
       dwidp=dsqrt(ap*aml/x)/dc
       dts=acos( (s*sx+ap2*y)/sqls/sqly )
       dtp=acos( (x*sx-ap2*y)/sqlx/sqly )
       dpi=pi


c
c      define bin boundaries
       dbtk(1)=0.d0
       dbtk(2)=dmin1(4.d0*dc*dwids	,	dts   /3.0d0)
       dbtk(3)=dmax1(	 dts-dwids	,	dts   /1.5d0)
       dbtk(4)=dmin1(	 dts+dwids	,dts+(dtp-dts)/3.0d0)
       dbtk(5)=dmax1(	 dtp-dwidp	,dts+(dtp-dts)/1.5d0)
       dbtk(6)=dmin1(	 dtp+dwidp	,dtp+(dpi-dtp)/3.0d0)
       dbtk(7)=.98d0*dbtk(6)+.02d0*dpi
       dbtk(8)=dpi
c
c      check cut in theta which follows from infra red cut:
c
       nbtk=0
       do 10 ka=2,8
       nbtk=nbtk+1
       if(dbtk(nbtk+1).ge.dtkmax) goto 20
10     continue
       goto 30
20     dbtk(nbtk+1)=dmin1(dbtk(nbtk+1),dtkmax)
30     continue
c
c      nbtk=number of valid bins
c
       return
       end

       subroutine inttk2(isumtk,dbtk1,dbtk2,dsumtk,derrtk)
      implicit real*8 (a-h,o-z)
c      integrate over ntk bins of brems-gamma theta angle dtk
#include "radinclude/mcpolrad.inc"
#include "radinclude/sxycom.inc"
#include "radinclude/cmpcom.inc"
#include "radinclude/tailcom.inc"
#include "radinclude/ppicom.inc"

       ddtk=(dbtk2-dbtk1)/ntk
       dtk=dbtk1-ddtk*.5d0
       dsum=0.d0
       derr1=0.d0
c      loop over all bins
       do 10 itk=1,ntk
       isumtk=isumtk+1
       dtk=dtk+ddtk
	 itkcur=isumtk
c      calculate integrand
c	call intsai(dtk,dtsai)
c	dsum=dsum+dtsai
	ta=(sx-sqly*cos(dtk))/ap2
       dsigma=an*alfa/pi*rv2(ta) * sin(dtk)*sqly/ap2
       dsum=dsum+dsigma
c	write(9,*)dsigma,ta,dtk
       if(ita.le.3)then
	 dsitkm(isumtk,ita) = dsumtk + ddtk * dsum
	 dcvtkm(isumtk,ita) = dtk
	 ddetkm(isumtk,ita) = ddtk
       endif
c      add up errors:
c      errors of theta integration:
       if    (itk.gt.2.and.itk.lt.ntk) then
		 derr1=derr1+dabs(dsigma-dold)
       elseif(itk.eq.2.or.itk.eq.ntk) then
		 derr1=derr1+1.5d0*dabs(dsigma-dold)
       endif
c
       dold=dsigma
c
10     continue
c
c      integral on big bin is:
c	print *,dsumtk,ddtk,dsum
       dsumtk=dsumtk+ddtk*dsum
       derrtk=derrtk+ddtk*derr1
c
       return
       end

       subroutine grndm(rm,nm)
       dimension rm(nm)
       common/seed/iy

       do n=1,nm
c	   rm(n)=rndm(n)
	  rm(n)=urand(iy)
       enddo
       end

c$nodebug

c * * * * * * * * * * * * * * * * * * * * * * * * * * *
c *						      *
      function urand(iy)
c *						      *
c *   this is a standard pseudo-random generator      *
c *   that work on ibm-370 and ibm-pc. we don't       *
c *   know does it work on sun? 		      *
c *						      *
c * * * * * * * * * * * * * * * * * * * * * * * * * * *
      integer*4 iy,m2,ia,ic
      data s,m2,ia,ic/.46566128e-9,1073741824,843314861,453816693/
      iy=iy*ia+ic
      if(iy.lt.0)iy=(iy+m2)+m2
      urand=float(iy)*s
      end

c$debug

       subroutine itafun(ys,xs,pl,pn,ixytest,itagen)
#include "radinclude/mcpolrad.inc"

       parameter (nt=35)
       common/xytab/x(nt,nt),y(nt)
     .	 ,sig1g_u(nt,nt),sigrad_u(nt,nt),tbor_u(nt,nt)
     .	 ,tine_u(nt,nt),tnuc_u(nt,nt),tpro_u(nt,nt)
     .	 ,sig1g_p(nt,nt),sigrad_p(nt,nt),tbor_p(nt,nt)
     .	 ,tine_p(nt,nt),tnuc_p(nt,nt),tpro_p(nt,nt)
       common/density/denstk(nt,nt,245,3),width(nt,nt,7,3)

       dimension xym(2),na(2),a(2*nt),wrk(nt,nt)
       dimension xx(nt),work(nt),w(nt),d(nt),xwork1(nt),xwork2(nt)
       dimension rnd(1)


       do nny=nt,1,-1
	  if(ys.ge.y(nny))goto 10
       enddo
 5     print *,' ys=',ys,' out of region',nny,nt,xs,pl,pn,ixytest,itagen
       print *,' ys=',ys,' out of region',xs,ys,'*',y
        stop
10     if(nny.eq.nt)goto 5

       if(ixytest.eq.1)then

	  stop ' ixytest.eq.1 is not supported in the version'


       else

	xym(1)=xs
	xym(2)=ys
	na(1)=nt
	na(2)=nt
	do i=1,nt
	   a(i)=x(i,1)
	   a(i+nt)=y(i)
	enddo
	 tboru=fint(2,xym,na,a,tbor_u)
	 tborp=fint(2,xym,na,a,tbor_p)
	 sigradu=fint(2,xym,na,a,sigrad_u)
	 sigradp=fint(2,xym,na,a,sigrad_p)
	    sig1gu=fint(2,xym,na,a,sig1g_u)
	    sig1gp=fint(2,xym,na,a,sig1g_p)

       tbor=tboru+pl*pn*tborp
       sigrad=sigradu+pl*pn*sigradp
	    sig1g=sig1gu+pl*pn*sig1gp

      endif

      sigcor=sigrad/sig1g

c      write(*,*)xs,ys,sigcor,sigrad,sig1g

      call grndm(rnd,1)
      r1=rnd(1)
      rr1=r1*sigrad
      if(rr1.gt.sigrad-tbor)then
	itagen=0
	return
      else

	    tineu=fint(2,xym,na,a,tine_u)
	    tinep=fint(2,xym,na,a,tine_p)
	    tnucu=fint(2,xym,na,a,tnuc_u)
	    tnucp=fint(2,xym,na,a,tnuc_p)
	    tprou=fint(2,xym,na,a,tpro_u)
	    tprop=fint(2,xym,na,a,tnuc_p)

	    tine=tineu+pl*pn*tinep
	    tnuc=tnucu+pl*pn*tnucp
	    tpro=tprou+pl*pn*tprop



	 endif


c	 write(*,'(6g13.5)')sigrad,tbor,tine,tpro,tnuc,rr1
c	 pause

	 if(rr1.gt.(tpro+tnuc)) then
	       itagen=1
c	       scgen=rr1-tpro-tnuc

	 elseif(rr1.gt.tnuc)then
	       itagen=3
c	       scgen=rr1-tnuc
	 else
c	       scgen=rr1
	       itagen=2
	 endif
	    do ii=1,7*ntk
	       do ix=1,nt
	       do iy=1,nt
		   wrk(ix,iy)=denstk(ix,iy,ii,itagen)
	       enddo
	       enddo
	       dsitkm(ii,itagen)=fint(2,xym,na,a,wrk)
	    enddo

	    do ii=1,7
	       do ix=1,nt
	       do iy=1,nt
		   wrk(ix,iy)=width(ix,iy,ii,itagen)
	       enddo
	       enddo
               do iittkk=1,ntk
	        ddetkm((ii-1)*ntk+iittkk,itagen)=fint(2,xym,na,a,wrk)
	       enddo
	    enddo

	    dcvtkm(1,itagen)=ddetkm(1,itagen)/2d0
       	 do iittkk=2,7*ntk
	    dcvtkm(iittkk,itagen)=dcvtkm(iittkk-1,itagen)+
     6     (ddetkm(iittkk-1,itagen)+ddetkm(iittkk,itagen))/2d0
         enddo


      end


      subroutine xytabl(tname,e1,plrun,pnrun,ixytest,ire)
#include "radinclude/mcpolrad.inc"

       parameter (nt=35)
       common/xytab/x(nt,nt),y(nt)
     .	 ,sig1g_u(nt,nt),sigrad_u(nt,nt),tbor_u(nt,nt)
     .	 ,tine_u(nt,nt),tnuc_u(nt,nt),tpro_u(nt,nt)
     .	 ,sig1g_p(nt,nt),sigrad_p(nt,nt),tbor_p(nt,nt)
     .	 ,tine_p(nt,nt),tnuc_p(nt,nt),tpro_p(nt,nt)

       common/density/denstk(nt,nt,245,3),width(nt,nt,7,3)
       character*13 tname

c	data lun/66/ymax/0.951/ymin/0.2/
c     .	    xmax/0.661/xmin/0.12/
c	data lun/66/ymax/0.951/ymin/0.2/
c     .	    xmax/0.661/xmin/0.12/
	data lun/66/ymax/0.99/ymin/0.05/
     .	    xmax/0.99/xmin/0.01/

	open (unit=lun,file=tname)

	if(ixytest.eq.0)then

	  write(lun,'(1x,g11.4)')e1
	  do iy=1,nt

	  y(iy)=ymin+(ymax-ymin)/(nt-1)*(iy-1)
	  ys=y(iy)

	  step=(xmax-xmin)/(nt-1)
	  do ix=1,nt
	     x(ix,iy)=xmin+step*(ix-1)
	     bmp=0.938272
	     bmc2=1.151857d0
cc	       ylim=(bmc2-bmp**2)/(2.*bmp*e1*(1d0-x(ix,iy)))
	     dmu=bmp**2/(2.*bmp*e1)
	     dmuc=(bmc2-bmp**2)/(2.*bmp*e1)
	     if(ymin.lt.dmu)then
	      print *,'ymin = ',ymin
	      print *,'ymin must be larger then',dmu
	      stop 'ymin is too small'
	     endif
	     xmaxi=min(1.-dmuc/ys,(1.-ys)/dmu/ys)-1d-4
	     xixiy=min(xmaxi,x(ix,iy))
c	     write(*,'(7g11.4)')xmaxi,xixiy,x(ix,iy)
c     .	   ,1.-dmuc/ys,(1.-ys)/dmu/ys,dmuc,dmu
	     iphi=0
	     call mpolrad(e1,y(iy),xixiy,1.,plrun,pnrun,-1)
	     tbor_u(ix,iy)=tbor
	     sig1g_u(ix,iy)=sig1g
	     sigrad_u(ix,iy)=sigrad
	     tine_u(ix,iy)=tine
	     tnuc_u(ix,iy)=tnuc
	     tpro_u(ix,iy)=tpro

c	      write(*,'(a3,6g12.4)')' u ',
c     .	       y(iy),x(ix,iy),sig1g,tnuc,tpro,sigrad

	     write(lun,'(1x,14g11.4)')x(ix,iy),y(iy)
     .	 ,sig1g_u(ix,iy),sigrad_u(ix,iy),tbor_u(ix,iy)
     .	 ,tine_u(ix,iy),tnuc_u(ix,iy),tpro_u(ix,iy)

	do itkm=0,6
	 width(ix,iy,itkm+1,1)=ddetkm(itkm*ntk+1,1)
	 width(ix,iy,itkm+1,2)=ddetkm(itkm*ntk+1,2)
	if(ire.ne.1)
     .	 width(ix,iy,itkm+1,3)=ddetkm(itkm*ntk+1,3)
	enddo

	ittmax=3
	if(ire.eq.1)ittmax=2
	write(lun,*)((width(ix,iy,itkm+1,itt)
     .	,itkm=0,6),itt=1,ittmax),ndxtkm

       do itkm=1,ntk*7
	denstk(ix,iy,itkm,1)=dsitkm(itkm,1)/dsitkm(ntk*7,1)
	denstk(ix,iy,itkm,2)=dsitkm(itkm,2)/dsitkm(ntk*7,2)
	if(ire.ne.1)
     .	denstk(ix,iy,itkm,3)=dsitkm(itkm,3)/dsitkm(ntk*7,3)
       enddo
	   write(lun,'(35g12.4)')(denstk(ix,iy,itkm,1)
     .				 ,itkm=1,7*ntk)
	   write(lun,'(35g12.4)')(denstk(ix,iy,itkm,2)
     .				 ,itkm=1,7*ntk)
	if(ire.ne.1)
     .	   write(lun,'(35g12.4)')(denstk(ix,iy,itkm,3)
     .				 ,itkm=1,7*ntk)
	  enddo
	  enddo
         print *,'RADGEN close lun->',lun
	 close(lun)

	elseif(ixytest.gt.0)then
	  read(lun,'(1x,g11.4)')e1tab
	 if(abs(e1-e1tab).gt.1d-5)then
	  write(*,*)' Energies: ',e1,e1tab
	  stop 'The table is prepared for another value of energy'
	 endif
	do iy=1,nt
	do ix=1,nt
	  print *,'radgen grid: filling ',iy,ix
	     read(lun,*)x(ix,iy),y(iy)
     .	 ,sig1g_u(ix,iy),sigrad_u(ix,iy),tbor_u(ix,iy)
     .	 ,tine_u(ix,iy),tnuc_u(ix,iy),tpro_u(ix,iy)
	ittmax=3
	if(ire.eq.1)ittmax=2
	read(lun,*)((width(ix,iy,ite,itt),ite=1,7),itt=1,ittmax),ndxtkm
	   read(lun,*)(denstk(ix,iy,itkm,1)
     .				 ,itkm=1,7*ntk)
	   read(lun,*)(denstk(ix,iy,itkm,2)
     .		       ,itkm=1,7*ntk)
	 if(ire.ne.1)
     .	   read(lun,*)(denstk(ix,iy,itkm,3)
     .		       ,itkm=1,7*ntk)


	enddo
	enddo
      else
	stop 'xytabl -- ixytest < 0'
      endif
        print *,'RADGEN close lunend-->',lun
	close(lun)
      end






       SUBROUTINE R1990(DQ2,DX,DR)
       implicit real*8(a-h,o-z)

       dr=r1990f(dx,dq2)

       end

