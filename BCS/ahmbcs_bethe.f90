!Set variables in common for the BCS solution
module BCS_COMMON
  implicit none
  integer             :: L
  real(8)             :: beta,u
  real(8)             :: n0
  real(8),allocatable :: epsi(:),dos(:),csi(:),Ep(:)
end module BCS_COMMON

program AHM_MF
  USE BCS_COMMON
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter      :: ndim=2,M=5000
  real(8)                :: x(ndim),fvec(ndim)
  real(8)                :: xmu
  real(8)                :: D,de,tol
  real(8)                :: n,delta,wmax,eps
  integer                :: i,ik,iflag,info
  real(8),allocatable    :: wr(:)
  complex(8)             :: zdet,zeta1,zeta2,x1,x2
  complex(8),allocatable :: fg(:,:),zeta(:)
  logical                :: printf
  external bcs_funcs
  namelist/bcsvars/L,beta,u,D,wmax,eps,tol,n0,printf
  call parse_cmd_variable(L,"L",default=5000)
  call parse_cmd_variable(beta,"BETA",default=1000.d0)
  call parse_cmd_variable(u,"U",default=1.d0)
  call parse_cmd_variable(D,"D",default=1.d0)
  call parse_cmd_variable(wmax,"WMAX",default=5.d0)
  call parse_cmd_variable(eps,"EPS",default=1.d-4)
  call parse_cmd_variable(tol,"TOL",default=1.d-15)
  call parse_cmd_variable(n0,"N0",default=1.d0)
  call parse_cmd_variable(printf,"PRINTF",default=.false.)
  write(*,nml=bcsvars)

  allocate(epsi(L),dos(L),csi(L),Ep(L))
  call bethe_lattice(dos,epsi,L,D)

  x=(/0.d0,0.001d0/)
  write(*,"(A5,4A16)")"Iter","mu","delta","f(mu)","f(delta)"
  call fsolve(bcs_funcs,x,tol,info)

  !Post-processing:
  csi  = epsi - x(1)
  Ep   = sqrt(csi**2+x(2)**2)
  n    = 1.d0-sum(dos*csi/Ep)
  xmu  = x(1)
  delta= x(2)

  allocate(wr(M),zeta(M))
  allocate(fg(2,M))
  wr = linspace(-wmax,wmax,M)
  fg=zero
  zeta = cmplx(wr,eps,8) + xmu
  do i=1,M
     zeta1 = zeta(i)
     zeta2 = conjg(zeta(M+1-i))
     x1 = 0.5d0*((zeta1+zeta2) + sqrt((zeta1-zeta2)**2 - 4.d0*delta**2 ))
     x2 = 0.5d0*((zeta1+zeta2) - sqrt((zeta1-zeta2)**2 - 4.d0*delta**2 ))
     fg(1,i) = zeta2/(x2-x1)*(gfbether(wr(i),x1,D)-gfbether(wr(i),x2,D))
     fg(2,i) =-delta/(x2-x1)*(gfbether(wr(i),x1,D)-gfbether(wr(i),x2,D))
     ! do ik=1,L
     !    zdet = (zeta1-epsi(ik))*(zeta2-epsi(ik)) + delta**2
     !    fg(1,i)=fg(1,i) + dos(ik)*(zeta2-epsi(ik))/zdet
     !    fg(2,i)=fg(2,i) - dos(ik)*delta/zdet
     ! enddo
  enddo
  call splot("DOS.bcs",wr,-dimag(fg(1,:))/pi,append=printf)
  call splot("G_realw.bcs",wr,fg(1,:),append=printf)
  call splot("F_realw.bcs",wr,fg(2,:),append=printf)
  call splot("observables.bcs",u,x(1),beta,n,abs(delta),append=printf)
end program AHM_MF


subroutine bcs_funcs(n,x,fvec,iflag)
  USE BCS_COMMON
  implicit none
  integer :: n
  real(8),intent(in) :: x(n)
  real(8) :: fvec(n)
  integer :: iflag,iter=0
  save iter
  iter=iter+1
  csi = epsi - x(1)
  Ep  = sqrt(csi**2+x(2)**2)
  fvec(1) = 1.d0-sum(dos*csi/Ep)-n0
  fvec(2) = u/2.d0*sum(dos*tanh(beta/2.d0*Ep)/Ep)-1.d0
  write(*,"(I5,4F16.9)")iter,x(1),x(2),fvec(1),fvec(2)
end subroutine bcs_funcs
