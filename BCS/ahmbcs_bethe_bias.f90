program AHM_MF
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter :: ndim=1
  !Model parameters:
  integer           :: L,M
  real(8)           :: beta,U,Wband,wmax,eps,tol,vbias
  logical           :: printf  
  real(8),allocatable :: epsi(:),dos(:),csi(:),Ep(:)
  real(8)                :: x(ndim),fvec(ndim)
  real(8)                :: xmu
  real(8)                :: de
  real(8)                :: n,delta,delta0,Egap
  integer                :: i,ik,iflag,info,Wgap(1)
  real(8),allocatable    :: wr(:)
  complex(8)             :: det,zeta1,zeta2,x1,x2
  complex(8),allocatable :: fg(:,:),zeta(:)
  integer :: iter=0
  real(8), parameter :: epsabs = 0.d00
  real(8), parameter :: epsrel = 1.d-8
  integer ier
  integer neval
  real(8) result
  real(8) abserr

  character(len=60) :: finput

  call parse_cmd_variable(finput,"FINPUT",default="inputBCS_VBIAS.conf")
  call parse_input_variable(L,"L",finput,default=20000)
  call parse_input_variable(beta,"BETA",finput,default=1000.d0)
  call parse_input_variable(delta0,"DELTA0",finput,default=0.01d0)
  call parse_input_variable(u,"U",finput,default=0.25d0)
  call parse_input_variable(Wband,"Wband",finput,default=1.d0)
  call parse_input_variable(wmax,"WMAX",finput,default=5.d0)
  call parse_input_variable(eps,"EPS",finput,default=1.d-4)
  call parse_input_variable(tol,"TOL",finput,default=1.d-15)
  call parse_input_variable(vbias,"vbias",finput,default=0d0)
  call parse_input_variable(printf,"PRINTF",finput,default=.false.)
  call save_input_file(finput)
  M=L
  allocate(epsi(L),dos(L),csi(L),Ep(L))
  call bethe_lattice(dos,epsi,L,Wband)
  de=epsi(2)-epsi(1)

  x=[delta0]
  write(*,"(A5,2A16)")"Iter","delta","f(delta)"
  call fsolve(bcs_funcs,x,tol,info)

  !Post-processing:
  delta= x(1)
  csi  = epsi
  Ep   = sqrt(csi**2+delta**2)
  n    = 1.d0-simps(-wband,wband,dos/de*csi/Ep)

  allocate(wr(M),zeta(M))
  allocate(fg(2,M))
  wmax=delta+1d0
  wr = linspace(-wmax,wmax,M,mesh=de)
  fg=zero
  zeta = cmplx(wr,eps,8) + xmu
  do i=1,M
     zeta1 = zeta(i)
     zeta2 = conjg(zeta(M+1-i))
     x1 = 0.5d0*((zeta1+zeta2) + sqrt((zeta1-zeta2)**2 - 4.d0*delta**2 ))
     x2 = 0.5d0*((zeta1+zeta2) - sqrt((zeta1-zeta2)**2 - 4.d0*delta**2 ))
     fg(1,i) = zeta2/(x2-x1)*(gfbether(wr(i),x1,Wband)-gfbether(wr(i),x2,Wband))
     fg(2,i) =-delta/(x2-x1)*(gfbether(wr(i),x1,Wband)-gfbether(wr(i),x2,Wband))
  enddo
  call splot("G_realw.bcs",wr,-dimag(fg(1,:))/pi,dreal(fg(1,:)),append=printf)
  call splot("F_realw.bcs",wr,fg(2,:),append=printf)

  Wgap = maxloc(-dimag(fg(1,M/2+1:))/pi)
  Egap = wr(M/2+1+Wgap(1))
  call splot("observables.bcs",u,vbias,n,abs(delta),abs(delta)/u/2,Egap,append=printf)

contains

  subroutine bcs_funcs(n,x,fvec,iflag)
    integer :: n
    real(8) :: x(n)
    real(8) :: fvec(n)
    integer :: iflag,iter=0
    real(8) :: func(L)
    save iter
    iter=iter+1
    delta=x(1)
    ! csi = epsi
    ! Ep  = sqrt(csi**2+x(1)**2)
    !fvec(1) = u/2.d0*sum(dos*tanh(beta/2.d0*Ep)/Ep)-1.d0
    !func = (1d0 - fermi(Ep+vbias/2d0,beta)-fermi(Ep-vbias/2d0,beta))/Ep
    ! fvec(1) = u/2.d0*sum(dos*func)-1.d0
    call qags (fintegrand, -wband, wband, epsabs, epsrel, result, abserr, neval, ier )
    fvec(1) = u/2*result-1d0
    write(*,"(I5,4F16.9)")iter,abs(x(1)),fvec(1)
  end subroutine bcs_funcs


  function fintegrand(e) result(f)
    real(8) :: e
    real(8) :: f
    real(8) :: Ep
    Ep  = sqrt(e**2+delta**2)
    f = (1d0 - fermi(Ep+vbias/2d0,beta)-fermi(Ep-vbias/2d0,beta))/Ep
    f = dens_bethe(e,wband)*f
  end function fintegrand

end program AHM_MF


