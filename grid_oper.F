C--------------------------------------- 
c attention  grids limited by 300x300 !
c
      subroutine reset_grids()
#include "xfgrid.inc"
         call VZERO(xgrid,90000)
         call VZERO(F1g,90000)
         call VZERO(F2g,90000)
         call VZERO(F3g,90000)
         call VZERO(Flg,90000)
         return
         end

C--------------------------------------- 
      subroutine write_json_file(itype)
C--------------------------------------- 
      implicit none
#include "options.inc"
      character*40  c0
      character*38  c1
      character*55 c5
      character*13  c2,c3
      character*14  c4
      character*6  c6
      character*150 xline,yline,zline,ptline
      integer itype,nxx,nyy,nzz,npt
      real xmin,xmax,zmin,zmax,ymin,ymax,ptmin,ptmax
c
         OPEN(unit=25,file='mydis.json',status='unknown')
c
      write(c0,'(A,F6.3)') '       "Beam Energy": ',cl_beam_energy
       nxx=cl_nxx
       nyy=cl_nyy
       xmin=cl_xmin
       xmax=cl_xmax
       ymin=cl_ymin
       ymax=cl_ymax
c
       write(25,'(A)') '{'
       if(itype.eq.0) then
       write(25,'(A)') '    "model": "F1F2009"'
       write(25,'(A)') '    "description": "x-section in 1/GeV^2"'
       write(25,'(A)') '    "reference": "Christi et al"'
       elseif (itype.eq.1) then
       write(25,'(A)') '    "model": "CJ15"'
       write(25,'(A)') '    "description": "x-section in 1/GeV^2"'
       write(25,'(A)') '    "reference": "N. Sato"'
       endif
c       
       write(25,'(A)') c0
       write(25,'(A)') '       "lepton-polarization": "0"'
       write(25,'(A)') '       "nucleon-polarization": "0"'
       write(25,'(A)') '       "variables":["x-sec,F1,F2,R"'
       write(25,'(A)') '       "axis":['
c
       write(c6,'(A)')  '  '
c
      if(nxx.lt.10) then
        write(c1,'(A,I1,A)')  '          {"name":"xb","bins": ',nxx,','
      else if(nxx.lt.100) then
        write(c1,'(A,I2,A)')  '          {"name":"xb","bins": ',nxx,','
      else
        write(c1,'(A,I3,A)')  '          {"name":"xb","bins": ',nxx,','
      endif
c
      write(c2,'(A,F5.3,A)') '"min": ',xmin,','
      write(c3,'(A,F5.3,A)') '"max": ',xmax,','
      write(c4,'(A)') '"scale":"lin",'
      write(c5,'(A)') '"description":"x_i=min+(2i-1)*(max-min)/2/Nx"}'

      write(xline,'(A,A,A,A,A,A)')  c6,c1,c2,c3,c4,c5
c
        write(c1,'(A,I3,A)')  '          {"name":" y","bins": ',nyy,','
c
      write(c2,'(A,F5.3,A)') '"min": ',ymin,','
      write(c3,'(A,F5.3,A)') '"max": ',ymax,','
      write(c4,'(A)') '"scale":"lin",'
      write(c5,'(A)') '"description":"y_j=min+(2j-1)*(max-min)/2/Nq"}'
      write(yline,'(A,A,A,A,A,A)') c6,c1,c2,c3,c4,c5

       write(25,'(A)') xline
       write(25,'(A)') yline
       write(25,'(A)') '    ],'
       write(25,'(A)') '       "parameters":['
       write(25,'(A)') '    ]'
       write(25,'(A)') '}'
c
      return
      end
C-----------------------------------------------------
      subroutine readgrid()
C-----------------------------------------------------
      implicit none
#include "options.inc"
#include "xfgrid.inc"
#include "const.inc"
      logical ex
      integer iostat,i,j,k,l,i1,j1,ix1,itotal
      real sxbmax,symax,ebeam
      real gam2,xa,ya,q2,q2a,xcalc,xcalc2,x_sec
      character*200  pdffile,clasdispdf
c
      inquire(file='xdata/xdata.json',exist=ex)
c
       if(.not.ex) then
        print *,'Error: missing file xdata/xdata.json'
        iostat=-100
        stop
       else ! read SFs
       call getenv('DISRAD_PDF',clasdispdf)
       print *,'READING the xdata/xdata.json for beam_E=',ebeam,
     6' From directory $DISRAD_PDF=',trim(clasdispdf)
       itotal=0
       write(pdffile,'(A,A)') trim(clasdispdf),'/xdata.json'
       open(unit=41,file=pdffile,status='OLD')
        read(41,*,IOSTAT=iostat,ERR=2) ebeam
         if(abs(ebeam-cl_beam_energy).gt.0.0001) then
          write(6,*),'Inconsistent set, stop!! ', ebeam,cl_beam_energy
          stop
         endif
          read(41,*,IOSTAT=iostat,ERR=2) nxx,sxbmax,sxbmin
          read(41,*,IOSTAT=iostat,ERR=2) nyy,symax,symin
          dxx=(sxbmax-sxbmin)/nxx
          dyy=(symax-symin)/nyy
          xsec_max=0.0
c
c         make the min values according to grid
c
       print *,'READING the xdata/xdata.json',nxx,nyy
 3     CONTINUE
       read(41,*,IOSTAT=iostat,ERR=2) i,j,xa,ya,q2a,
     6 F2g(i+1,j+1),Flg(i+1,j+1),F3g(i+1,j+1),xgrid(i+1,j+1)
       q2=2.0*am*ya*cl_beam_energy*xa
       gam2=4.0*am2*xa*xa/q2
       xcalc=6.2829*alpha**2/xa/ya/q2a*(2*F2g(i+1,j+1)
     6*(1.0-ya+ya*ya/2.0+ya*ya*gam2/4.0)-ya*ya*Flg(i+1,j+1)
     6-xa*F3g(i+1,j+1))
       xcalc2=6.2829*alpha**2/xa/ya/q2a*(2*F2g(i+1,j+1)
     6*(1.0-ya+ya*ya/2.0+ya*ya*gam2/4.0)-ya*ya*Flg(i+1,j+1))
       F1g(i+1,j+1)=((1.0+gam2)*F2g(i+1,j+1)-Flg(i+1,j+1))/2.0/xa
         write(6,*) i,j,xa,ya,q2a, 
     6F2g(i+1,j+1),Flg(i+1,j+1),F3g(i+1,j+1),xgrid(i+1,j+1),xcalc,
     6xgrid(i+1,j+1)/xcalc,xcalc2/xcalc,F1g(i+1,j+1)
       itotal=itotal+1
       x_sec=xgrid(i+1,j+1)
       if(xsec_max.lt.x_sec) xsec_max=x_sec         ! define the maximum just in case
       goto 3   ! read until end of file
 2     CONTINUE
        write(6,*),'Total #ov lines = ',itotal
        close(41)
       endif
      return
      end
c
c
c
      logical function  xsec_accept(x,y)
#include "xfgrid.inc"
      real rlim,vallim,xsec_grid
      real x,y
      vallim=xsec_grid(x,y)/xsec_max
      rlim=rand(0)
c      print *,'in  xsec_accept',vallim,rlim,xsec_max
      if(vallim.gt.rlim) then
      xsec_accept=.TRUE.
      else
      xsec_accept=.FALSE.
      endif
      return
      end
c
c
      real function  xsec_grid(x,y)
      implicit none
#include "xfgrid.inc"
#include "options.inc"
      integer ix1,jy1,ix,jy
      real x,y,ddx,ddy,frx,frq2,xsec,rF1,rF2,rr,rFt,rFl
       xsec=0.0
       ddx=x-sxbmin
       ddy=y-symin
c
c
       ix=ddx/dxx+1
       jy=ddy/dyy+1
       if(ix.lt.1.or.ix.gt.nxx.or.jy.lt.1.or.jy.gt.nyy) return
       ix1=ix+1
       if(ix1.gt.nxx) ix1=nxx
       jy1=jy+1
       if(jy1.gt.nyy) jy1=nyy
       xsec_grid=xgrid(ix,jy)
       return
       end

      subroutine  calc_griddis(x,y,xsec,rF1,rF2,rFl,rr)
      implicit none
#include "xfgrid.inc"
#include "options.inc"
      integer ix1,jy1,ix,jy
      real x,y,ddx,ddy,frx,frq2,xsec,rF1,rF2,rr,rFt,rFl
       xsec=0.0
       rF1=0.0
       rF2=0.0
       
       ddx=x-sxbmin
       ddy=y-symin
c
c
       ix=ddx/dxx+1
       jy=ddy/dyy+1
       if(ix.lt.1.or.ix.gt.nxx.or.jy.lt.1.or.jy.gt.nyy) return
       ix1=ix+1
       if(ix1.gt.nxx) ix1=nxx
       jy1=jy+1
       if(jy1.gt.nyy) jy1=nyy
c
       rF1=F1g(ix,jy)
       rF2=F2g(ix,jy)
       rFt=2.0*x*rF1
       rFl=Flg(ix,jy)
       xsec=xgrid(ix,jy)
       rr=0.0
       if(rFt.gt.0) rr=rFl/rFt ! R=L/T ratio
       return
       end

c
c      write the grid
c
       subroutine write_grid(itype,ztar,atar,E_beam)
       implicit none
#include "options.inc"
#include "const.inc"
       real ztar,atar,E_Beam,q2,w2,nu,eps,gam
       real x_sec,rF1,rF2,rr,F1p,F2p,sf,Fl
       real dxx,dyy,x,y
       integer itype,iix,iiy
c
c     test
c
        x_sec=0.0
         call write_json_file(itype)
         print *,'Write the grid for \sigma>1.E-10 and exit!!!'
         dxx=(cl_xmax-cl_xmin)/cl_nxx
         dyy=(cl_ymax-cl_ymin)/cl_nyy
        do iix=1,cl_nxx
        x=cl_xmin+(2*iix-1)*dxx/2.0
c
        do iiy=1,cl_nyy
        y=cl_ymin+(2*iiy-1)*dyy/2.0
c
c       define kin variables
c
        q2=(2*am*x*E_beam*y)
        nu=y*E_beam
        w2=am*am-Q2+2*am*nu 
        gam=q2/nu/nu
c
        eps=(1-y-0.25*gam*gam*y*y)/(1-y+0.5*y*y+0.25*y*y*gam*gam) 
c 
        call get_xsec(ztar,atar,x,y,q2,w2,eps,gam,x_sec,F1p,F2p,rr,sf)
c
        if(x_sec.gt.1.E-10)  then
          Fl=(1+gam*gam)*F2p-2.0*x*F1p
          print *,'checkk', iix-1,iiy-1,x,y,F2p,Fl,x_sec,F1p
     6,'**',q2,nu,gam,eps
          write(25,101) iix-1,iiy-1,x,y,q2,F2p,Fl,rr,x_sec
        endif
 101    FORMAT(2I10,10E11.4)
c
        enddo
        enddo
        close(25)
        stop
        return
        end

c
c      fill the grid and store the maximum of x-section
c
       subroutine fill_grid(itype,ztar,atar,E_beam)
       implicit none
#include "xfgrid.inc"
#include "options.inc"
#include "const.inc"
       real ztar,atar,E_Beam,q2,w2,nu,eps,gam
       real x_sec,rF1,rF2,rr,F1p,F2p,sf,Fl
       real x,y
       integer itype,iix,iiy
c
        x_sec=0.0
        xsec_max=0.0

         print *,'Fill the grid for x-section!!!'
         nxx=cl_nxx
         dxx=(cl_xmax-cl_xmin)/nxx
         sxbmin=cl_xmin
c
         if(cl_xq2) then
          nyy=cl_nyy
          dyy=(cl_qmax-cl_qmin)/nyy
          symin=cl_qmin
         else
          nyy=cl_nyy
          dyy=(cl_ymax-cl_ymin)/nyy
          symin=cl_ymin
         endif
c
        do iix=1,nxx
        x=sxbmin+(2*iix-1)*dxx/2.0
c
        do iiy=1,nyy
        if(cl_xq2) then ! loop on Q^2    
         q2=symin+(2*iiy-1)*dyy/2.0
         y=q2/(2*am*x*E_beam)
        else
         y=symin+(2*iiy-1)*dyy/2.0
         q2=(2*am*x*E_beam*y)
        endif
c
c       define kin variables
c
        q2=(2*am*x*E_beam*y)
        nu=y*E_beam
        w2=am*am-Q2+2*am*nu 
        gam=q2/nu/nu
c
        eps=(1-y-0.25*gam*gam*y*y)/(1-y+0.5*y*y+0.25*y*y*gam*gam) 
c 
        call get_xsec(ztar,atar,x,y,q2,w2,eps,gam,x_sec,F1p,F2p,rr,sf)
        if(xsec_max.lt.x_sec) xsec_max=x_sec   !  define the maximum x_section
c
        if(x_sec.gt.1.E-12)  then
         xgrid(iix,iiy)=x_sec
        else
         xgrid(iix,iiy)=0.0
        endif
        print *,'grid-fill',cl_xq2,iix,iiy,x,y,q2,dxx,dyy,xgrid(iix,iiy)
c
        enddo
        enddo
c        stop
        return
        end

