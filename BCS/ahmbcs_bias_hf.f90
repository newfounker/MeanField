!Set variables in common for the BCS solution
module BCS_COMMON
  USE COMMON_VARS
  USE FUNCTIONS
  USE INTEGRATE
  USE MATRIX
  implicit none
  integer             :: L
  real(8)             :: beta,u,phi
  real(8)             :: n0
  real(8),allocatable :: epsi(:),dos(:),Ep(:),num1(:),num2(:),gamma1(:)
end module BCS_COMMON

program AHM_MF
  USE BCS_COMMON
  USE PARSE_CMD
  USE TOOLS
  USE IOTOOLS
  USE ZEROS
  implicit none
  integer,parameter      :: ndim=1,M=5000
  real(8)                :: x(ndim),fvec(ndim)
  real(8)                :: xmu
  real(8)                :: D,de,tol
  real(8)                :: n,delta,wmax,eps
  integer                :: i,ik,iflag,info
  real(8),allocatable    :: wr(:)
  complex(8)             :: det,zeta1,zeta2,x1,x2
  complex(8),allocatable :: fg(:,:),zeta(:)
  logical                :: printf
  external bcs_funcs
  namelist/bcsvars/L,beta,u,D,wmax,eps,tol,phi,printf
  call parse_cmd_variable(L,"L",default=5000)
  call parse_cmd_variable(beta,"BETA",default=1000.d0)
  call parse_cmd_variable(u,"U",default=1.d0)
  call parse_cmd_variable(D,"D",default=1.d0)
  call parse_cmd_variable(wmax,"WMAX",default=5.d0)
  call parse_cmd_variable(eps,"EPS",default=1.d-4)
  call parse_cmd_variable(tol,"TOL",default=1.d-15)
  call parse_cmd_variable(phi,"PHI",default=0.d0)!;phi=phi/2.d0
  call parse_cmd_variable(printf,"PRINTF",default=.false.)
  write(*,nml=bcsvars)

  
  allocate(epsi(L),dos(L),Ep(L),num1(L),num2(L),gamma1(L))
  call bethe_lattice(dos,epsi,L,D)

  x=(/0.001d0/)
  !write(*,"(A5,2A16)")"Iter","delta","f(delta)"
  call fsolve(bcs_funcs,x,tol,info)

  !Post-processing:
  delta= x(1)
  Ep   = sqrt(epsi**2+delta**2)
  n    = 1.d0-sum(dos*epsi/Ep)

  allocate(wr(M),zeta(M))
  allocate(fg(2,M))
  wr = linspace(-wmax,wmax,M)
  fg=zero
  zeta = cmplx(wr,eps,8)
  do i=1,M
     zeta1 = zeta(i)
     zeta2 = conjg(zeta(M+1-i))
     x1 = 0.5d0*((zeta1+zeta2) + sqrt((zeta1-zeta2)**2 - 4.d0*delta**2 ))
     x2 = 0.5d0*((zeta1+zeta2) - sqrt((zeta1-zeta2)**2 - 4.d0*delta**2 ))
     fg(1,i) = zeta2/(x2-x1)*(gfbether(wr(i),x1,D)-gfbether(wr(i),x2,D))
     fg(2,i) =-delta/(x2-x1)*(gfbether(wr(i),x1,D)-gfbether(wr(i),x2,D))
     ! do ik=1,L
     !    det = (zeta1-epsi(ik))*(zeta2-epsi(ik)) + delta**2
     !    fg(1,i)=fg(1,i) + dos(ik)*(zeta2-epsi(ik))/det
     !    fg(2,i)=fg(2,i) - dos(ik)*delta/det
     ! enddo
  enddo
  call splot("DOS.bcs",wr,-dimag(fg(1,:))/pi,append=printf)
  call splot("G_realw.bcs",wr,fg(1,:),append=printf)
  call splot("F_realw.bcs",wr,fg(2,:),append=printf)
  call splot("observables.bcs",phi,delta,u,beta,n,append=printf)
end program AHM_MF


subroutine bcs_funcs(n,x,fvec,iflag)
  USE BCS_COMMON
  implicit none
  integer :: n,i,ik,j
  real(8) :: x(n),delta,sumj
  real(8) :: fvec(n),int,e,eg,n1,n2,gp,gm,g1,nf,uk,vk
  integer :: iflag,iter=0  
  save iter
  iter=iter+1
  delta=x(1)
  int=0.d0
  do i=1,L
     e  = epsi(i)
     eg = sqrt(delta**2/2.d0+phi**2+e**2 + 0.5d0*sqrt(delta**4+(2.d0*delta*phi)**2+(4.d0*phi*e)**2))
     nf = fermi(Eg,beta)
     n1 = (phi**2+(e-eg)**2) * (phi**2-(e+Eg)**2)**2 + delta**2*(e-Eg)**2*(phi**2+(e+Eg)**2)
     n2 = (phi**2+(e+eg)**2) * (phi**2-(e-Eg)**2)**2 + delta**2*(e+Eg)**2*(phi**2+(e-Eg)**2)
     gp = phi**2-(e+Eg)**2
     gm = phi**2-(e-Eg)**2
     g1 = gp*gm
     int=int-dos(i)/sqrt(n1*n2)*( &
          (e**2-Eg**2)*(g1+phi**2*delta**2)*(1.d0-nf) + &
          (delta**2*(e**2-Eg**2)+phi**2*g1)*nf&
          )
  enddo
  fvec(1)=int*u/delta-1.d0
  !write(*,"(I5,4F16.9)")iter,x(1),fvec(1)
end subroutine bcs_funcs
