program HM_MF
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter      :: n=1
  integer                :: p,q
  real(8)                :: x(1),e
  real(8)                :: xmu,dens
  real(8)                :: D,de,tol
  real(8)                :: wmax,eps
  integer                :: i,ik,iflag,info
  integer                :: Lw,Le
  real(8)                :: temp,u,sum
  real(8)                :: n0
  real(8),allocatable    :: epsi(:),dos(:),mesh(:)
  real(8),allocatable    :: wr(:)
  complex(8),allocatable :: fg(:),zeta(:)
  character(len=32)      :: finput
  integer                :: imesh
  integer :: iter=0
  real(8), parameter :: epsabs = 0.d00
  real(8), parameter :: epsrel = 1.d-8
  integer ier
  integer neval
  real(8) result
  real(8) abserr


  call parse_cmd_variable(finput,"FINPUT",default="inputHFHM.in")
  call parse_input_variable(Le,"Le",finput,default=10000,comment="Size of the energy mesh.")
  call parse_input_variable(P,"P",finput,default=20,comment=".")
  call parse_input_variable(Q,"Q",finput,default=500,comment=".")
  call parse_input_variable(Lw,"Lw",finput,default=5000,comment="Size of the frequency mesh.")
  call parse_input_variable(temp,"TEMP",finput,default=0.01d0,comment="Temperature, 0.0 is not admitted.")
  call parse_input_variable(u,"U",finput,default=0.d0,comment="Interaction strenght.")
  call parse_input_variable(D,"D",finput,default=1.d0,comment="Half-bandwidth, default=1.")
  call parse_input_variable(wmax,"WMAX",finput,default=3.d0,comment="Largest frequency value for GF calculation.")
  call parse_input_variable(eps,"EPS",finput,default=2.d-3,comment="Spectral density broadening.")
  call parse_input_variable(tol,"TOL",finput,default=1.d-9,comment="Tolerance for non-linear equation solution.")
  call parse_input_variable(n0,"N0",finput,default=1.d0,comment="Target density.")
  call parse_input_variable(imesh,"imesh",finput,default=0,comment="Select mesh type: 0=linear, 1=power-law (see P,Q), 2=log")
  call save_input_file(finput)

  select case(imesh)
  case default
     allocate(epsi(Le),dos(Le),mesh(Le))     
     epsi = linspace(-D,D,Le,istart=.false.,iend=.false.,mesh=de)
  case (1)
     Le=2*P*Q+1
     allocate(epsi(Le),dos(Le),mesh(Le))
     epsi = upminterval(-D+1.d-4,D-1.d-4,0.d0,P,Q)
  case(2)
     allocate(epsi(Le),dos(Le),mesh(Le))
     epsi(Le/2+1:Le)=(D-1.d-3)-logspace(D-1.d-3,0.d0,Le/2)
     epsi(1:Le/2)=-epsi(Le/2+1:Le)
     call sort_array(epsi)
  end select

  print*,"Using Le=",Le


  !SETUP THE 1-DIMENSIONAL DOS
  do i=1,Le
     dos(i)=1.d0/pi/sqrt(D**2-epsi(i)**2)*heaviside(D-abs(epsi(i)))
  enddo
  dos = dos/trapz(epsi,dos)
  open(100,file="dos1d.dat")
  do i=1,Le
     write(100,*)epsi(i),dos(i)
  enddo
  close(100)


  !SOLVE THE SELF-CONSISTENCY CONDITION FINDING THE ZERO OF THE USER DEFINED FUNCTION BELOW
  x=0.d0
  dens=1.d0
  write(*,"(A5,4A16)")"Iter","mu","f(mu)"
  call fsolve(hfhm_func,x,tol,info)
  if(info/=1)print*,"INFO from hybrd=",info

  !POST-PROCESSING: PRINT THE GREEN'S FUNCTION
  allocate(wr(Lw),zeta(Lw),fg(Lw))
  xmu  = x(1)
  call qags (fintegrand, -D, D, epsabs, epsrel, dens, abserr, neval, ier )
  !dens = 2.d0*trapz(epsi,dos(:)*fermi(epsi(:)-xmu+U*n0,1.d0/temp))
  wr   = linspace(-wmax,wmax,Lw)
  fg   = zero
  zeta = cmplx(wr,eps,8) + xmu - U*n0/2d0
  do i=1,Lw
     fg(i) = trapz(epsi,dos(:)/(zeta(i)-epsi(:)))
  enddo
  open(100,file="GF.dat")
  open(101,file="observables.dat")
  do i=1,Lw
     write(100,*)wr(i),-dimag(fg(i))/pi,dreal(fg(i))
  enddo
  write(101,"(4F15.8)")u,x(1),1.d0/temp,2d0*dens

contains

  function hfhm_func(x) result(fvec)
    real(8) :: x(:)
    real(8) :: fvec(size(x))
    integer :: iter=0
    real(8), parameter :: epsabs = 0.d00
    real(8), parameter :: epsrel = 1.d-8
    integer ier
    integer neval
    real(8) result
    real(8) abserr
    save iter
    iter=iter+1
    xmu = x(1)
    call qags (fintegrand, -D, D, epsabs, epsrel, result, abserr, neval, ier )
    fvec(1) =  n0-2d0*result !trapz(epsi,dos(:)*fermi(epsi(:)-xmu+U*n0/2d0,1.d0/temp))
    write(*,"(I5,4F16.9)")iter,x(1),fvec(1),2*result
  end function hfhm_func

  function fintegrand(e) result(f)
    real(8) :: e
    real(8) :: f
    f = 1.d0/pi/sqrt(D**2-e**2)*heaviside(D-abs(epsi(i)))*fermi(e-xmu + U*n0/2d0,1.d0/temp)
  end function fintegrand
end program HM_MF


