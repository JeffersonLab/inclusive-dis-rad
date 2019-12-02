	real	polrad_xmin, polrad_xmax, polrad_ymin, polrad_ymax
	parameter	(polrad_xmin=0.01)
	parameter	(polrad_xmax=0.99)
	parameter	(polrad_ymin=0.02)
	parameter	(polrad_ymax=0.99)

      double precision sigrad,tine,tnuc,tpro,tbor,demin,phipoi
     . ,sig1g,sigcor,vac,vertex,small,redfac
     . ,dsts,dcts
     . ,taout,rrout,dsitkm,dcvtkm,ddetkm,dsigmr,drcurr,ddeler
      integer ntk,nrr,itkcur,iphi,ndxtkm
      COMMON /mcpolrad/sigrad,tine,tnuc,tpro,tbor,demin
     . ,sig1g,sigcor,vac,vertex,small,redfac
     . ,dsts,dcts
     . ,phipoi,taout,rrout
     . ,dsitkm(400,3),dcvtkm(400,3),ddetkm(400,3)
     . ,dsigmr(200,400),drcurr(200,400),ddeler(200,400)
     . ,ntk,nrr ,itkcur,iphi,ndxtkm(3)
