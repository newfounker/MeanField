program HM_MF
  USE CONSTANTS
  USE PARSE_INPUT
  USE IOTOOLS
  USE OPTIMIZE, only: fsolve
  USE FUNCTIONS, only: fermi,bethe_lattice
  USE INTEGRATE, only: trapz,simps
  USE ARRAYS, only: linspace
  implicit none
  integer,parameter      :: n=1
  real(8)                :: x(1)
  real(8)                :: xmu,dens
  real(8)                :: D,de,tol
  real(8)                :: wmax,eps
  integer                :: i,ik,iflag,info
  integer                :: Lw,Le
  real(8)                :: temp,u
  real(8)                :: n0
  real(8),allocatable    :: epsi(:),dos(:)
  real(8),allocatable    :: wr(:)
  complex(8),allocatable :: fg(:),zeta(:)
  character(len=32)      :: finput

  call parse_cmd_variable(finput,"FINPUT",default="inputHFHM.in")
  call parse_input_variable(Le,"Le",finput,default=5000,comment="Size of the energy mesh.")
  call parse_input_variable(Lw,"Lw",finput,default=5000,comment="Size of the frequency mesh.")
  call parse_input_variable(temp,"TEMP",finput,default=0.001d0,comment="Temperature, 0.0 is not admitted.")
  call parse_input_variable(u,"U",finput,default=0.d0,comment="Interaction strenght.")
  call parse_input_variable(D,"D",finput,default=1.d0,comment="Half-bandwidth, default=1.")
  call parse_input_variable(wmax,"WMAX",finput,default=5.d0,comment="Largest frequency value for GF calculation.")
  call parse_input_variable(eps,"EPS",finput,default=1.d-4,comment="Spectral density broadening.")
  call parse_input_variable(tol,"TOL",finput,default=1.d-15,comment="Tolerance for non-linear equation solution.")
  call parse_input_variable(n0,"N0",finput,default=1.d0,comment="Target density.")

  allocate(epsi(Le),dos(Le))
  call bethe_lattice(dos,epsi,Le,D)
  de = epsi(2)-epsi(1)

  x=(/0.d0/)
  dens=1.d0
  write(*,"(A5,4A16)")"Iter","mu","f(mu)"
  call fsolve(hfhm_func,x,tol,info)

  !Post-processing:
  allocate(wr(Lw),zeta(Lw),fg(Lw))
  xmu  = x(1)
  dens = 2.d0*simps(de,dos(:)/de*fermi(epsi(:)-xmu+U*n0,1.d0/temp))
  wr   = linspace(-wmax,wmax,Lw)
  fg   = zero
  zeta = cmplx(wr,eps,8) + xmu - U*n0
  do i=1,Lw
     fg(i) = trapz(de,dos(:)/de/(zeta(i)-epsi(:)))
  enddo
  call splot("GF.dat",wr,-dimag(fg(:))/pi,dreal(fg(:)))
  call splot("observables.dat",u,x(1),1.d0/temp,dens,append=.true.)

contains

  function hfhm_func(x) result(fvec)
    real(8) :: x(:)
    real(8) :: fvec(size(x))
    integer :: iter=0
    save iter
    iter=iter+1
    fvec(1) =  n0/2.d0-simps(de,dos(:)/de*fermi(epsi(:)-x(1)+U*n0,1.d0/temp))
    write(*,"(I5,4F16.9)")iter,x(1),fvec(1)
  end function hfhm_func

end program HM_MF


