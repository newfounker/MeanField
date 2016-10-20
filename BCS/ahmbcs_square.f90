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
  real(8)                :: ts,de,tol,delta0
  real(8)                :: n,delta,wmax,eps,kx,ky
  integer                :: i,ix,iy,ik,iflag,info,Nkx,Lk
  real(8),allocatable    :: wr(:)
  complex(8)             :: zdet,zeta1,zeta2,x1,x2
  complex(8),allocatable :: fg(:,:),zeta(:)
  logical                :: printf,bool

  external bcs_funcs

  call parse_input_variable(L,"L","inputBCS.conf",default=5000)
  call parse_input_variable(Nkx,"NKX","inputBCS.conf",default=100)
  call parse_input_variable(beta,"BETA","inputBCS.conf",default=1000.d0)
  call parse_input_variable(u,"U","inputBCS.conf",default=1.d0)
  call parse_input_variable(ts,"TS","inputBCS.conf",default=0.5d0)
  call parse_input_variable(wmax,"WMAX","inputBCS.conf",default=5.d0)
  call parse_input_variable(eps,"EPS","inputBCS.conf",default=0.01d0)
  call parse_input_variable(tol,"TOL","inputBCS.conf",default=1.d-15)
  call parse_input_variable(n0,"N0","inputBCS.conf",default=1.d0)
  call parse_input_variable(delta0,"DELTA0","inputBCS.conf",default=1.d-2)
  call parse_input_variable(printf,"PRINTF","inputBCS.conf",default=.false.)
  call save_input_file("inputBCS.conf")

  Lk=Nkx**2
  allocate(epsi(Lk),dos(Lk),csi(Lk),Ep(Lk))
  ik=0
  do ix=1,Nkx
     kx = -pi + (ix-1)*pi2/Nkx
     do iy=1,Nkx
        ky = -pi + (iy-1)*pi2/Nkx
        ik=ik+1
        epsi(ik) = -2d0*ts*(cos(kx)+cos(ky))
     enddo
  enddo
  dos = 1d0/Lk
  call get_free_dos(epsi,dos)

  xmu=0d0
  delta=delta0
  x=[xmu,delta]
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
     do ik=1,Lk
        zdet = (zeta1-epsi(ik))*(zeta2-epsi(ik)) + delta**2
        fg(1,i)=fg(1,i) + dos(ik)*(zeta2-epsi(ik))/zdet
        fg(2,i)=fg(2,i) - dos(ik)*delta/zdet
     enddo
  enddo
  call splot("DOS.bcs",wr,-dimag(fg(1,:))/pi,append=printf)
  call splot("G_realw.bcs",wr,fg(1,:),append=printf)
  call splot("F_realw.bcs",wr,fg(2,:),append=printf)
  call splot("observables.bcs",u,abs(delta),xmu,beta,n,append=printf)
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
