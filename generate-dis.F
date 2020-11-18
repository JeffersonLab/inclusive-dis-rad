c
c   generates d\sigma/dE'd\Omega in nanobarn/Gev/steradian
c
       program gendis
       implicit none
#include "options.inc"
#include "const.inc"
#include "xfgrid.inc"
      integer iq,iw,ix,nskip,Nevents
      real q2,x,y,w2,rc,F1p,F2p
      real rndm,s2
      real E_beam,ztar,atar,nu,tmp
      data ztar/1.0/,atar/1.0/
      real eps,gam,xmin,xmax,ymin,ymax,q2min,q2max,thetar,fiel
      real px_el,py_el,pz_el,e_prime,p_prime,x_sec,pwidth,zwidth
      real xcalc,q2calc,rF1,rF2,rr,sf
      integer itype
      logical is_accepted,xsec_accept
c     radgen stuff
      integer ixytst,lst40
      character*2  CTARGET
c
       integer ifail,iaccept
      data ifail/0/,iaccept/0/
c
      call RNDMSEED()
      call cloptions()
       E_beam=cl_beam_energy
      if(cl_rad.gt.0) then    ! using table
       if(cl_target.eq.1) CTARGET='P'                ! can be also D,HI
       if(cl_target.eq.2) CTARGET='D'                ! can be also HI
       if(cl_target.eq.3) CTARGET='HI'                !helium
       LST40=0               ! unpolarized for now
c
        if(cl_rad.eq.2) then 
         ixytst=-1              ! calculate event by event
        else if(cl_rad.eq.3) then
         ixytst=0               ! create the table on flight
         if(cl_xmin.lt.0.01.or.cl_xmax.gt.0.99.or.
     6    cl_ymin.lt.0.055.or.cl_ymax.gt.0.99) then
          print *,'x,y-limits incompatible with radgen'
     6,   'x ',cl_xmin,cl_xmax,' y ',cl_ymin,cl_ymax
          stop
         endif
	 print *,'RADGEN WILL CREATE THE TABLE!',cl_rad
         call radgen_init(ctarget,E_beam,lst40,ixytst)
         print *,'RADGEN table created for Eb=',E_beam
         ixytst=2               ! create the table table
       else if(cl_rad.eq.1) then
         ixytst=2               ! reading from table
       call radgen_init(ctarget,E_beam,lst40,ixytst)
       endif
      endif
c
c
      if(cl_gridxsec.or.cl_gridsf) call reset_grids()  ! reset to 0 all grids      
c
c
      if(cl_gridxsec) call readgrid()   ! read x-sect-grid
      if(cl_gridsf) call readgrid()   ! read x-sect-grid
      if(cl_realxsec) then            ! will use grid to reject
       if(.not.cl_gridxsec) call fill_grid(itype,ztar,atar,E_beam) ! fill a new grid
       print *,'Filled xgrid for x-section'
      endif
       
c
      call get_datfilename()
c
      if (cl_iset.gt.0) then 
          print *,'Initialize-start CJ-set ',cl_iset
          call setCJ(cl_iset)
      endif
c
      if(cl_gridwrite) then   ! dump the grid and exit
        if (cl_iset.gt.0) then ! dump CJ
         itype=1
        else
         itype=0               ! dump F1F2
        endif
       call write_grid(itype,ztar,atar,E_beam)
      endif
c
      xmin=cl_xmin
      xmax=cl_xmax
      ymin=cl_ymin
      ymax=cl_ymax
      q2min=cl_qmin
      q2max=cl_qmax
      pwidth=cl_trad
      zwidth= cl_zwidth
      Nskip=cl_skip
      Nevents=cl_triggers
c
        iq=0
 100    continue
c
        x = xmin+(xmax-xmin)*rand(0)
        if(cl_xq2) then ! random Q^2    
         q2 = q2min+(q2max-q2min)*rand(0)
         y=q2/(2*am*x*E_beam)
        else            ! random y 
         y = ymin+(ymax-ymin)*rand(0)
         q2=(2*am*x*E_beam*y)
        endif
        fiel=2.0*pi*rand(0)
c
        nu=y*E_beam
        w2=am*am-Q2+2*am*nu 
        if(w2.lt.cl_wmin.or.q2.gt.cl_qmax
     6.or.q2.lt.cl_qmin.or.y.lt.ymin.or.y.gt.ymax) goto 100         ! not below proton 
c define the electron
        e_prime=E_beam-nu
        p_prime=sqrt(e_prime*e_prime-em*em)
        s2=q2/(4.*E_beam*e_prime)
        thetar=asin(sqrt(s2))*2.0
        if(s2.gt.1.0.or.s2.lt.0) goto 100         ! 0<sin^2<1
        px_el=p_prime*sin(thetar)*cos(fiel)
        py_el=p_prime*sin(thetar)*sin(fiel)
        pz_el=p_prime*cos(thetar)
c-------------------------------
        gam=q2/nu/nu
        q2calc=4*E_beam*e_prime*(sin(thetar/2.0))**2
        xcalc=q2calc/(2.0*am*nu)
        eps=(1-y-0.25*gam*gam*y*y)/(1-y+0.5*y*y+0.25*y*y*gam*gam) 
c
c       get F1,F2
c
        call get_xsec(ztar,atar,x,y,q2,w2,eps,gam,x_sec,F1p,F2p,rr,sf)

        if(x_sec.gt.1.E-12) then  ! write only x_sections above e-12
c
          if(cl_xq2) then 
            is_accepted=xsec_accept(x,q2)
          else
            is_accepted=xsec_accept(x,y)
          endif
c
          if(cl_realxsec.and..not.is_accepted) then
               ifail=ifail+1
c               print *,'-iq,x,y,accept',iq,x,y,ifail,iaccept
               goto 100   ! not good enough for real generation
          else
               iaccept=iaccept+1
c               print *,'+iq,x,y,accept',iq,x,y,ifail,iaccept
          endif
       call writeevent(x,q2,w2,px_el,py_el,pz_el,e_prime,x_sec
     6,sf,pwidth,zwidth)
        iq=iq+1
        endif !xsec
c
        if(iq.ge.Nevents) then
           goto 2  ! end the loop
        else 
           goto 100 ! continue looping
        endif !iq
c
c
 2       continue
               print *,'iq,ifail,iaccept',iq,ifail,iaccept
       stop
       end
c
       subroutine get_xsec(ztar,atar,x,y,q2,w2,eps,gam,x_sec,F1p,F2p
     6,rr,sf)
       implicit none
#include "options.inc"
#include "const.inc"
        real ztar,atar,x,y,q2,w2,eps,gam
        real x_sec,rF1,rF2,rFl,F1p,F2p,rr,sf
c
        if(cl_gridxsec) then ! get xsec from grid
          if(cl_xq2) then
            call calc_griddis(x,q2,x_sec,rF1,rF2,rFl,rr)
          else
            call calc_griddis(x,y,x_sec,rF1,rF2,rFl,rr)
          endif
          F1p=rF1
          F2p=rF2
          sf=2*y*y*x*F1p+2*(1.0-y-0.25*gam**2*y*y)*F2p
        else !calculate xsec from F1F2
          if (cl_gridsf) then    ! get F1f2 from grid
           call calc_griddis(x,y,x_sec,rF1,rF2,rFl,rr)
           F1p=rF1
           F2p=rF2
          else
           call calc_f1(ztar,atar, x, q2,w2, F1p, F2p,rr,gam)
          endif
          sf=2*y*y*x*F1p+2*(1.0-y-0.25*gam**2*y*y)*F2p
          x_sec=389379.0*alpha**2*sf/(q2*y*x)
        if(cl_xq2) then 
          x_sec=x_sec*y/q2                 ! add Jakobian if in x-vs-Q^2 instead y
        endif

c     6    (F1p+eps*(1+gam**2)/2/x/(1-eps)*F2p) 
        endif
c
        return
        end
c
c
c
      subroutine get_datfilename()
      implicit none
#include "options.inc"
c      character*4  cpol,cdis
      character*4  cpol
      character*5  cdis
      character*10  cxsec
      character*22  ccuts
c
c
ccc
      print *,'make the data file name',cl_beam_energy
c
      if(cl_beam_energy.lt.10.0) then
      write(ccuts,'(A,F5.3,A,F3.2,A,F4.2)') 'e',cl_beam_energy,
     6'.emin',cl_emin,'tmin',cl_tmin 
      else
      write(ccuts,'(A,F6.3,A,F4.2,A,F4.2)') 'e',cl_beam_energy,
     6'.emn',cl_emin,'tmn',cl_tmin 
      endif
c
      if(cl_rad.gt.0) then
      write(datfilename,'(A)') 'clas-dis-rad.'//ccuts//'.'
      else
      write(datfilename,'(A)') 'clas-dis.'//ccuts//'.'
      endif
      print *,'*********',datfilename
      return
      end
c
c  write to LUND file
c
      subroutine writeevent(x,q2,w2,px,py,pz,eprime,xsec
     6,sf,pwidth,zwidth)
      implicit none
#include "options.inc"
#include "const.inc"
      real x,q2,y,w2
      real px,py,pz,eprime,pprime,xsec,sf,pwidth,zwidth,vx,vy,vz
      integer   NumNuc,NumProt,Ntracks,beam_pid,targ_id,process_id
      integer irannom
      data beam_pid/11/,targ_id/2212/,process_id/0/ ! ep scattering
      data NumNuc/1/,NumProt/1/,Ntracks/1/                     ! proton target for now
      real e1,p1,vpgen(4),VPRAD(4),PhRAD(4),q2tr,utr,weight,vdotm
      real rand,rnd
c
          logical lFirst
          DATA lFirst /.TRUE./
          integer jc,m,itotal
          DATA itotal /0/
          character*8  c8name
          character*200  fname
          integer*4 iran,now(3),ifilenum
          DATA ifilenum /1/
c
          cl_pol=1.0    ! default values for weight
          cl_epol=0.0   ! default value for pol/sf
          if(cl_target.eq.2) NumNuc=2
          if(cl_target.eq.3) NumNuc=3
          if(cl_rad.gt.0) process_id=cl_rad  
c        
          if(itotal.ge.cl_nmax) then   ! each file correspond to number ~ xsection
           close (41)
           lFirst=.TRUE.
           ifilenum=ifilenum+1
           itotal=0
          endif
c
                itotal=itotal+1
          if (lFirst) then   ! open the file
           if(cl_docker) then
             write(fname,'(A)') 'inclusive-dis-rad.dat'
             cl_nmax=cl_triggers+100
           else
           print *,'filename ',datfilename
c
c    define the  file name
c
             jc = 0
             do m=1,len(datfilename)
                if ((datfilename(m:m).eq.' ').and.(jc.eq.0)) then
                   jc=m-1
                endif
             enddo
             if (jc.eq.0) jc=len(datfilename)
c
           if(ifilenum.le.9) then
             write(c8name,'(A,i1,A)') '000',ifilenum,'.dat'
           else if(ifilenum.ge.10.and.ifilenum.le.99) then
             write(c8name,'(A,i2,A)') '00',ifilenum,'.dat'
           else if(ifilenum.ge.100.and.ifilenum.le.999) then
             write(c8name,'(A,i3,A)') '0',ifilenum,'.dat'
           else
             write(c8name,'(i4,A)')  ifilenum,'.dat'
c stay with numbers below 9999 for now
           endif 
c
             write(fname,'(A,A)') datfilename(1:jc),c8name
             write(6,'(A,A)') 'Data file Name ',fname
          endif    ! docker check
             OPEN(unit=41,file=fname,status='unknown')
             lFirst=.FALSE.
          endif
c

           vx=(0.5-rand(0))*pwidth
           vy=(0.5-rand(0))*pwidth
           vz=cl_zpos+(0.5-rand(0))*zwidth
c
        if(cl_rad.gt.0) then
         e1=cl_beam_energy
         p1=sqrt(e1**2-aml2)
         pprime=sqrt(eprime**2-aml2)
c
         vpgen(1)=  -px   !pprime*sin(theta)*cos(phi)
         vpgen(2)=  -py   !pprime*sin(theta)*sin(phi)
         vpgen(3)=p1-pz   !pprime*cos(theta)
         vpgen(4)=e1-eprime

         y=q2/(2.0*am*x*e1)
c$$$       print *,'calling RADGEN',e1,eprime,x,y,vpgen(4),q2
c$$$     6,vpgen,vdotm(vpgen,vpgen,4),vprad,phrad
        
        if(y.gt.0.99) then
         print *,'not possible!!!!',y,q2,x,e1,am,cl_ymax 
         y=cl_ymax
        endif
        call RADGEN(e1,VPGEN,VPRAD,PhRAD,q2tr,utr,weight)
c
             if(0.lt.weight.and.100.gt.weight) then
              cl_pol=weight  ! temporary solution
             else
	      print *,'PROBLEM WITH WEIGHT?',itotal,weight,
     6        x,q2,w2,y,e1,eprime,' ** ',vpgen,'** ',PhRAD,'**',q2tr,utr
              itotal=itotal-1
              return
             endif
c
             if (Phrad(4).gt.0.001)  then
                Ntracks=2    ! only E_gamma>0.001 GeV
             else
               Ntracks=1
             endif
           xsec=xsec*weight
        endif                  ! end cl_rad>0
c
           cl_epol=sf  ! temporary solution
c
           write(41,101) Ntracks,NumNuc,NumProt,cl_pol,cl_epol
     6,beam_pid,cl_beam_energy,targ_id,process_id,xsec,sf,q2tr,utr
           write(41,102) 1,-1.0, 1, 11, 0, 0
     6 ,Px,Py,Pz,eprime,0.0005,vx,vy,vz
           if(Ntracks.eq.2) then   ! write also the photon
           write(41,102) 2,0.0, 1, 22, 1, 0
     6 ,Phrad(1),Phrad(2),Phrad(3),Phrad(4),0.0000,vx,vy,vz
           endif


            
c$$$           write(6,102) 1,-1.0, 1, 11, 0, 0
c$$$     6 ,Px,Py,Pz,e,0.0005,vx,vy,vz

 101      FORMAT(2x,I10,2I6,2F11.2,I5,F10.3,2I5,2E14.7,2F6.2)
 102      FORMAT(2x,I3,F4.0,I5,I7,2I3,5F10.4,2X,5F10.4)
      return
      end

       SUBROUTINE CLOPTIONS()
       implicit none
#include "options.inc"
       integer i
       character*20 cnumber
       INTEGER  numopts,ivalnum
       real valnum
c
c   grid selection
c
      cl_gridxsec=.FALSE.  ! don't read x-section grid
      cl_gridsf=.FALSE.  ! don't read struck function grid
      cl_gridwrite=.FALSE.  ! don't write the grid
      cl_realxsec=.FALSE.  ! do the weighted generation by default
      cl_xq2=.FALSE.       ! do x/y by default
      cl_docker=.FALSE.    ! set if 1 file is needed 
c
      cl_zpos=-3.0           ! z-position of RGA
      cl_zwidth=5.0          ! z-width ( zpos+/-zwidth/2)
      cl_trad=0.5             ! target radius
c
      cl_beam_energy=10.6   ! def eg1-2000 (e16-5.754,e1f-5.498)
      cl_emin=0.75           ! def e'min
      cl_tmin=0.05           ! e16->0.3           ! def e'tmin
      cl_tmax=1.57           ! def e'tmax
c
      cl_nxx=200
      cl_nyy=200
      cl_xmax=0.99
      cl_ymax=0.99
      cl_qmin=0.95
      cl_qmax=18.05
      cl_wmin=0.93827
      cl_wmax=400.0
c
      cl_triggers=10      ! all events
      cl_seed=0
      cl_nmax = 2000       ! max number of events in the file
      cl_skip = 0         ! number to skip
      cl_iset=0
      cl_rad=0            ! no radiative effects,1-std radgen
      cl_target=1           ! target-proton, 2-deuteron, 3-He
      cl_radsfs=1           ! RC model 1-radgen,2-Bosted
      numopts  = IARGC()
c      
      if(numopts.gt.0 ) then
        i=0
      else
        goto 2
      endif
 1    i=i+1
        CALL GETARG(i,cnumber)
        if(cnumber.eq.'--trig'.and.i.lt.numopts) then
          i=i+1
          CALL GETARG(i,cnumber)
          cl_triggers=valnum(cnumber)
        elseif(cnumber.eq.'--seed'.and.i.lt.numopts) then
          i=i+1
          CALL GETARG(i,cnumber)
          cl_seed=ivalnum(cnumber)
        elseif(cnumber.eq.'--maxev'.and.i.lt.numopts) then
          i=i+1
          CALL GETARG(i,cnumber)
          cl_maxev=valnum(cnumber)
        elseif(cnumber.eq.'--skip'.and.i.lt.numopts) then
          i=i+1
          CALL GETARG(i,cnumber)
          cl_skip=valnum(cnumber)
        elseif(cnumber.eq.'--nmax'.and.i.lt.numopts) then
          i=i+1
          CALL GETARG(i,cnumber)
          cl_nmax=valnum(cnumber)
        elseif(cnumber.eq.'--beam'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_beam_energy=valnum(cnumber)
        elseif(cnumber.eq.'--t'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_tmin=valnum(cnumber)*3.1415926/180.0
           i=i+1
           CALL GETARG(i,cnumber)
           cl_tmax=valnum(cnumber)*3.1415926/180.0
        elseif(cnumber.eq.'--e'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_emin=valnum(cnumber)
           i=i+1
           CALL GETARG(i,cnumber)
           cl_emax=valnum(cnumber)
        elseif(cnumber.eq.'--rad'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_rad=valnum(cnumber)
        elseif(cnumber.eq.'--targ'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_target=valnum(cnumber)
        elseif(cnumber.eq.'--radsfs'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_radsfs=valnum(cnumber)
        elseif(cnumber.eq.'--x'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_xmin=valnum(cnumber)
           i=i+1
           CALL GETARG(i,cnumber)
           cl_xmax=valnum(cnumber)
        elseif(cnumber.eq.'--y'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_ymin=valnum(cnumber)
           i=i+1
           CALL GETARG(i,cnumber)
           cl_ymax=valnum(cnumber)
        elseif(cnumber.eq.'--q'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_qmin=valnum(cnumber)
           i=i+1
           CALL GETARG(i,cnumber)
           cl_qmax=valnum(cnumber)
        elseif(cnumber.eq.'--nxy'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_nxx=valnum(cnumber)
           i=i+1
           CALL GETARG(i,cnumber)
           cl_nyy=valnum(cnumber)
        elseif(cnumber.eq.'--w'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_wmin=valnum(cnumber)
           i=i+1
           CALL GETARG(i,cnumber)
           cl_wmax=valnum(cnumber)
        elseif(cnumber.eq.'--zpos'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_zpos=valnum(cnumber)
        elseif(cnumber.eq.'--set'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_iset=valnum(cnumber)
        elseif(cnumber.eq.'--trad'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_trad=valnum(cnumber)
        elseif(cnumber.eq.'--zwidth'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_zwidth=valnum(cnumber)
        elseif(cnumber.eq.'--realxsec') then
           cl_realxsec=.TRUE.
        elseif(cnumber.eq.'--xgrid') then
           cl_gridxsec=.TRUE.
        elseif(cnumber.eq.'--sfgrid') then
           cl_gridsf=.TRUE.
        elseif(cnumber.eq.'--writegrid') then
           cl_gridwrite=.TRUE.
        elseif(cnumber.eq.'--xq2') then
           cl_xq2=.TRUE.
        elseif(cnumber.eq.'--docker') then
           cl_docker=.TRUE.
        else
          call printusage
          stop
        endif
        if(i.lt. numopts) goto 1
 2      continue
c
        if(cl_emax.le.0.0.or.cl_emax.gt.cl_beam_energy) then
         cl_emax=0.95*cl_beam_energy            ! put some limit
        endif
c      cl_xmin=1.0/2.0/cl_beam_energy  ! for 1 GeV^2 Q^2min,y_max=1
c      cl_ymin=1.0/2.0/cl_beam_energy  ! for 1 GeV^2 Q^2min
        if(cl_xmin.lt.0.055) cl_xmin=0.055   ! radgen couldn't hadle that
        if(cl_xmax.gt.0.99) cl_xmax=0.99   ! radgen couldn't hadle that
        if(cl_ymin.lt.0.055) cl_ymin=0.055  ! radgen couldn't hadle that
        if(cl_ymax.gt.0.99) cl_ymax=0.99  ! radgen couldn't hadle that
        call printvalues
c
c randomize
c
        call RNDMSEED()
        return 
        end
c
        subroutine printusage
        print *,'clasDIS [options]            '
        print *,'      option  value      default  comment  '
        print *,'  --v    verbos_level    0     additional printout '
        print *,'  --trig nevents         10     number of triggers '
        print *,'  --beam beam_energy     10.6   e- momentum in GeV '
c23456789012345678901234567890123456789012345678901234567890123456789012  
        print *,'  --e e-min e-max  min[750] max[0.95*e] Eprime in GeV '
        print *,'  --t tmin tmax  14.0    theta min/max for \theta(deg)'
        print *,'  --x xmin xmax        0.05 0.95    x_b min/max'
        print *,'  --y ymin ymax        0.05 0.99    x_b min/max'
        print *,'  --q Q2min Q2max      0.85 16.05 Q2 min max in GeV^2'
        print *,'  --w W2min W2max      4.0  50     W^2 min max GeV^2'
        print *,'  --nxy nxx nyy   200 200  number of bins'
        print *,'  --zpos z-position     -3    target z position in cm '
        print *,'  --zwidth z-width 5.0 width in z cm (zpos+/-zwidth/2)'
        print *,'  --set  iset  400    400-CJ15 LO PDFs'
        print *,'  --ski  0   number of rndm to skip'
        print *,'  --nmax   200   number of events in a file'
        print *,'  --xq2  .FALSE.   generate in x/y'
        print *,'  --xgrid  .FALSE.   use xsec grid'
        print *,'  --sfgrid  .FALSE.   use  F1F2 grid'
        print *,'  --realxsec .FALSE. use real xsection'
        print *,'  --xq2  .FALSE. use x/y generation by def.'
        print *,'  --writegrid  .FALSE.  dump the grid'
        print *,'  --rad  0   no radiation,1-read-grid,2-calc,3-create'
        print *,'  --targ  1  proton,2-deut,3-He'
        print *,'  --radsfs  1  1-radgen,2-bosted'
        return
        end
c
        subroutine printvalues
        implicit none
        include 'options.inc'
        print *,'clasDIS input values:'    
        print *,'NEVENTS:             ',cl_triggers    
        print *,'Maxevents in file:   ',cl_maxev    
        print *,'nmax-events/per-file   : ',cl_nmax
        print *,'BEAM ENERGY:         ',cl_beam_energy  
        print *,'Eprime   min/max:         ',cl_emin,cl_emax  
        print *,'Eprime tmin/tmax:         ',cl_tmin,cl_tmax  
        print *,'xmin/xmax:         ',cl_xmin,cl_xmax  
        print *,'q2min/q2max:         ',cl_qmin,cl_qmax  
        print *,'ymin/ymax:         ',cl_ymin,cl_ymax  
        print *,'rc-mode:         ',cl_rad 
        print *,'xgrid:  use xsecion grid',cl_gridxsec
        print *,'sfgrid  use F1/F2 grid',cl_gridsf
        print *,'realxsec .TRUE. use real xsection',cl_realxsec
        print *,'xq2 .TRUE. use x/Q^2 generation',cl_xq2
c
        print *,'TARGET type:   ', cl_target 
        print *,'RC-Model 1-radgen/2-bosted:   ', cl_radsfs 
        print *,'TARGET z-position:   ', cl_zpos 
        print *,'TARGET 1/2 width :   ', cl_zwidth
        print *,'TARGET radius :      ', cl_trad
        print *,'__________________________________'
        return
        end
c
c
c
      integer function length(str)
c     return the string length without the blanks characters

      implicit integer (k-l)
      character *(*) str    
      
      lmax=len(str)      
       
c     search the last non blank character
      do i=lmax,1,-1
      if(str(i:i).ne.' ')then
      length=i
      return
      end if
      end do

      length=lmax

      return
      end     

      real function valnum(str)
c     return the real value contained into a string

      implicit integer (k-l)
      character *(*) str
      logical segno
      

      segno=.false.
      valnum=0.00
      lu=length(str)
      
c     check the number sign
      if(str(1:1).eq.'-')then
      segno=.true.
      str=str(2:lu)
      lu=lu-1
      end if

c     check if number is float or integer
      if(index(str,'.').ne.0)then
      iin=index(str,'.')-1
      else
      iin=lu
      end if

      ifr=lu-(iin+1)

c     translate the integer portion
      doi=1,iin  
      k=ichar(str(i:i))-48
      valnum=valnum+float(k)*10.00**float(iin-i)
      end do         

      if(iin.eq.lu)goto 10
      str=str(iin+2:lu)

c     translate the decimal portion
      doi=1,ifr  
      k=ichar(str(i:i))-48
      valnum=valnum+float(k)/10.00**float(i)
      end do

10    if(segno)valnum=-valnum

      return
      end     
c 
      integer function ivalnum(str)
c     return the real value contained into a string

      implicit integer (k-l)
      character *(*) str
      logical segno
      double precison valnum
      

      segno=.false.
      valnum=0.00
      lu=length(str)
      
c     check the number sign
      if(str(1:1).eq.'-')then
      segno=.true.
      str=str(2:lu)
      lu=lu-1
      end if

c     check if number is float or integer
      if(index(str,'.').ne.0)then
      iin=index(str,'.')-1
      else
      iin=lu
      end if

      ifr=lu-(iin+1)

c     translate the integer portion
      doi=1,iin  
      k=ichar(str(i:i))-48
      valnum=valnum+float(k)*10.00**float(iin-i)
      end do         

      if(iin.eq.lu)goto 10
      str=str(iin+2:lu)

c     translate the decimal portion
      doi=1,ifr  
      k=ichar(str(i:i))-48
      valnum=valnum+float(k)/10.00**float(i)
      end do

10    if(segno)valnum=-valnum
       ivalnum=int(valnum)
      return
      end     
c
      real function vdotm(a,b,n)
       implicit none
       real a(4),b(4),s
       integer i,n
       s=0.0
       do i=1,3
         s=s+a(i)*b(i)
       enddo
       if(n.eq.4) s=s-a(n)*b(n)
       vdotm=s
       return
       end
c
c
       subroutine RNDMSEED()
       implicit none
#include "options.inc"
       integer*4 now(3),ifilenum
       integer iran,irannom
       real rand,rnd
c
c
          if(cl_seed.eq.0) then
          call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
	  irannom=now(3)+now(1)*3600+now(2)*60
	  else 
	  irannom=cl_seed
	  endif
           CALL SRAND(irannom)
          print *,'Rndm numb->',irannom,now(1),now(2),now(3),rand()
c           DO 100 iran=1,irannom
c 100       RND=rand()
c
	   return
	   end
