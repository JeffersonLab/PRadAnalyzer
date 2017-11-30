      integer nv,nt1,nz
      parameter(nt1=30)
      parameter(nz=60)
      common/merad_grid/grt1(nt1),grz(nz)
      common/merad_dist/distsit1(4*nt1),distart1(4*nt1),
     .                  distsiz(nz),distarz(nz)
      real*8 grt1,grz,distsit1,distart1,distsiz,distarz

