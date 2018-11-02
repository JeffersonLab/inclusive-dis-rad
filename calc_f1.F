         subroutine calc_f1(Z, A, x, Q2, w2, F1, F2, rc,gam)
         implicit none
#include "options.inc"
#include "const.inc"
        real Z,A,x,q2,w2,F1, F2, rc,gam
        real*8 ZD,AD,qsq,wsq,F1p, F2p, r
        double precision xd,qd,CJpdf
        real f1u,f1d,f1ub,f1db,F1calc,F2calc,slst,rsl
        real cj_h0,cj_h1,cj_h2,F2wht
        data cj_h0/-3.2874/,cj_h1/1.9274/,cj_h2/-2.07/
c
         ZD=Z
         AD=A
         QSQ=Q2
         WSQ=w2 
         call F1F2IN09(ZD, AD, QSQ, Wsq, F1p, F2p,r)
         F1=F1p
         F2=F2p
         rc=r
c
c        compare with CJ15
c
       if(cl_iset.gt.0) then
       xd=x
       Qd=sqrt(Q2)
       f1u=CJpdf(1,xd,Qd)
       f1d=CJpdf(2,xd,Qd)
       f1ub=CJpdf(-1,xd,Qd)
       f1db=CJpdf(-2,xd,Qd)
       F1calc=0.5*(4*(f1u+f1ub)+(f1d+f1db))/9.0
       F2calc=xd*(4*(f1u+f1ub)+(f1d+f1db))/9.0
       F2wht=1.0+cj_h0*x**cj_h1*(1.0+cj_h2*x)/Q2
       rsl=(1.0+gam*gam)*F2/(2.0*x*F1)-1
       print *,'x,q2-- ',x,q2,w2,f1,f1calc,f2,rsl,f1u
     6,f1d,f1ub,f1db,F2calc,F2wht
c
c      replacing F_2 with model
c
       F2=F2calc
       endif

         

c        call F1F2IN09(0.D0, 1.D0, q2, w2, F1n, F2n,rc)
c        call F1F2IN09(1.D0, 2.D0, q2, w2, F1d, F2d,rc)
c        call F1F2IN09(4.D0, 9.D0, q2, w2, F1be, F2be,rbe)
         return
         end
