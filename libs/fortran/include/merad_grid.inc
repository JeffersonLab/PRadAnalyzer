      integer nv,nt1,nz
      parameter(nv=300)
      parameter(nt1=30)
      parameter(nz=60)
      common/merad_grid/grv(nv),grt1(nt1),grz(nz)
      common/merad_dist/distsiv(nv),distarv(nv),
     .                  distsit1(4*nt1),distart1(4*nt1),
     .                  distsiz(nz),distarz(nz)
      real*8 grv,grt1,grz,distsiv,distarv,
     .       distsit1,distart1,distsiz,distarz

