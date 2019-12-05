      real cl_pol,cl_epol
      common /dispol/cl_pol,cl_epol
      logical cl_gridxsec,cl_gridsf,cl_gridwrite
     6,cl_realxsec,cl_xq2,cl_docker
      common /disgrids/cl_gridxsec,cl_gridsf,cl_gridwrite
     6,cl_realxsec,cl_xq2,cl_docker
c
      INTEGER cl_nxx,cl_nyy
      common /disxdump/cl_nxx,cl_nyy
      REAL cl_beam_energy,cl_zpos,cl_zwidth,cl_trad
      REAL  cl_emin,cl_emax,cl_tmin,cl_tmax,cl_ymin,cl_ymax
      REAL  cl_xmin,cl_xmax,cl_qmin,cl_qmax,cl_wmin,cl_wmax
      integer cl_triggers,cl_skip,cl_nmax,cl_iset
      common /disrealb/cl_beam_energy,cl_zpos,cl_zwidth,cl_trad
      common /disreal/  cl_xmin,cl_xmax,cl_qmin,cl_qmax,
     6cl_wmin,cl_wmax,cl_ymin,cl_ymax,cl_tmin,cl_tmax,cl_emin,cl_emax
      common /disint/ cl_triggers,cl_skip,cl_nmax,cl_iset 
c
      integer cl_maxev,total_events
      character*180 datfilename
      common /elasfile/ cl_maxev,total_events,datfilename
c
      integer cl_rad
      common/disrad/cl_rad
      integer cl_target,cl_radsfs
      common/distarg/cl_target,cl_radsfs
