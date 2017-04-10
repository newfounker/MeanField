program bdg_haldane
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter                           :: Norb=1,Nspin=1,Nlat=2
  integer                                     :: Nso,Nlso
  integer                                     :: Nloop
  integer                                     :: Lmats,Lreal
  integer                                     :: Lk
  !variables for the model:
  integer                                     :: Nk,Nkpath
  real(8)                                     :: Uloc
  real(8)                                     :: ts,tsp,phi,delta,Mh,wmixing,xmu,eps
  real(8)                                     :: beta
  real(8)                                     :: nsc
  real(8)                                     :: deltasc
  real(8)                                     :: wini,wfin
  real(8)                                     :: nread,tol
  logical                                     :: igetgf  
  !
  complex(8),allocatable,dimension(:,:,:,:)   :: Hk
  complex(8),allocatable,dimension(:,:,:)     :: HkNambu
  real(8)                                     :: a,a0,bklen,chern,kx,ky
  real(8),allocatable,dimension(:,:)          :: kvec
  real(8),dimension(:),allocatable        :: Wtk
  !
  real(8),dimension(2)                        :: d1,d2,d3
  real(8),dimension(2)                        :: a1,a2,a3
  real(8),dimension(2)                        :: bk1,bk2,pointK,pointKp
  real(8),dimension(:),allocatable            :: nii,pii,xarray   ![Nlat]
  integer                                     :: ik,i,is,ix,iy,info
  character(len=50)                           :: finput
  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats,Greal,Smats,Sreal ![2][Nlat][Norb][Norb][L]


  ! READ INPUT FILES !
  call parse_cmd_variable(finput,"FINPUT",default="inputBDG.conf")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS",finput,default=1d0)
  call parse_input_variable(tsp,"TSP",finput,default=1.d0/3/sqrt(3d0))
  call parse_input_variable(mh,"MH",finput,default=0d0)
  call parse_input_variable(phi,"PHI",finput,default=0d0)
  call parse_input_variable(uloc,"ULOC",finput,default=1d0,comment="Values of the local interaction")
  call parse_input_variable(nread,"NREAD",finput,default=0d0,comment="Target density")
  call parse_input_variable(xmu,"XMU",finput,default=0d0)
  call parse_input_variable(eps,"EPS",finput,default=1d-3)
  call parse_input_variable(beta,"BETA",finput,default=1000d0)
  call parse_input_variable(tol,"TOL",finput,default=1.d-5)
  call parse_input_variable(nsc,"NSC",finput,default=0.5d0,comment="Value of the initial density per spin (1/2=half-filling).")
  call parse_input_variable(deltasc,"DELTASC",finput,default=0.02d0,comment="Value of the SC symmetry breaking term.")
  call parse_input_variable(nloop,"NLOOP",finput,default=500,comment="Max number of iterations.")
  call parse_input_variable(Lmats,"LMATS",finput,default=2000,comment="Number of Matsubara frequencies.")
  call parse_input_variable(Lreal,"LREAL",finput,default=2000,comment="Number of real-axis frequencies.")
  call parse_input_variable(wini,"WINI",finput,default=-15.d0,comment="Smallest real-axis frequency")
  call parse_input_variable(wfin,"WFIN",finput,default=15.d0,comment="Largest real-axis frequency")
  call parse_input_variable(igetgf,"IGETGF",finput,default=.true.,comment="Flag to evaluate or not GF")


  call save_input_file(finput)


  ! SET THE COMPRESSION THRESHOLD TO 1Mb (1024Kb)
  call set_store_size(1024)

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")
  phi=phi*pi

  Nso=Nspin*Norb
  Nlso=Nlat*Nso

  a=1d0
  a0=a*sqrt(3d0)

  !Lattice basis (a=1; a0=sqrt3*a) is:
  !\a_1 = a0 [ sqrt3/2 , 1/2 ]
  !\a_2 = a0 [ sqrt3/2 ,-1/2 ]
  !
  !nearest neighbor: A-->B, B-->A
  d1= a*[  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= a*[  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= a*[ -1d0     , 0d0           ]
  !
  !
  !next nearest-neighbor displacements: A-->A, B-->B
  a1 = d2-d3                    !a*sqrt(3)[sqrt(3)/2,-1/2]
  a2 = d3-d1                    !a*sqrt(3)[-sqrt(3)/2,-1/2]
  a3 = d1-d2                    !a*sqrt(3)[0, 1]
  !
  !
  !RECIPROCAL LATTICE VECTORS:
  bklen=4d0*pi/sqrt(3d0)
  bk1=bklen*[ sqrt(3d0)/2d0 ,  1d0/2d0 ]
  bk2=bklen*[ sqrt(3d0)/2d0 , -1d0/2d0 ]
  !
  !
  pointK = [2*pi/3, 2*pi/3/sqrt(3d0)]
  pointKp= [2*pi/3,-2*pi/3/sqrt(3d0)]


  ! BUILD THE RECIPROCAL VECTOR GRID !
  Lk=Nk*Nk
  write(*,*)"# of k-points     :",Lk
  allocate(kvec(2,Lk))
  ik=0
  do iy=1,Nk
     ky = dble(iy-1)/Nk
     do ix=1,Nk
        ik=ik+1
        kx=dble(ix-1)/Nk
        kvec(:,ik) = kx*bk1 + ky*bk2
     enddo
  enddo



  ! SOLVE BDG PROBLEM !
  allocate(nii(Nlso))
  allocate(pii(Nlso))

  nii=nsc
  pii=deltasc

  allocate(xarray(Nlso+1))
  xarray(1)=xmu
  xarray(2:)=pii
  call read_data("xmu_pii.restart",xarray)
  xmu=xarray(1)
  pii=xarray(2:)
  deallocate(xarray)

  if(nread==0d0)then     
     call solve_BdG(bool=.true.)
  else
     allocate(xarray(1))
     write(*,"(A5,4A16)")"Iter","xmu","f(xmu)"
     xarray(1) = xmu
     call fsolve(search_mu_BdG,xarray,tol,info)
     xmu = xarray(1)
     call solve_BdG(bool=.true.)
  endif


  call get_LocalObservables()

  allocate(Hknambu(2*Nlso,2*Nlso,Lk))
  do ik=1,Lk
     Hknambu(:,:,ik) =  hk_NambuHaldane_model(kvec(:,ik),2*Nlso)
  enddo
  call get_Chern_number(Hknambu,[Nk,Nk],2,Nk/pi2*Nk/pi2,Chern)


  call Build_Eigenbands()


  call store_data("xmu_pii.restart",[xmu,pii])


  if(igetgf)then
     allocate(Gmats(2,Nlat,Norb,Norb,Lmats),Smats(2,Nlat,Norb,Norb,Lmats))
     allocate(Greal(2,Nlat,Norb,Norb,Lreal),Sreal(2,Nlat,Norb,Norb,Lreal))
     allocate(Hk(2,Nlso,Nlso,Lk))
     allocate(Wtk(Lk))
     Smats=zero
     Sreal=zero
     Gmats=zero
     Greal=zero
     forall(i=1:Nlso)
        Smats(1,:,i,i,:)=-Uloc*(nii(i)-0.5d0)
        Sreal(1,:,i,i,:)=-Uloc*(nii(i)-0.5d0)
        Smats(2,:,i,i,:)=-Uloc*pii(i)
        Sreal(2,:,i,i,:)=-Uloc*pii(i)
     end forall
     !
     do ik=1,Lk
        Hk(1,:,:,ik)    =     hk_haldane_model(kvec(:,ik),Nlso)
        Hk(2,:,:,ik)    =    -transpose(hk_haldane_model(-kvec(:,ik),Nlso))
     enddo
     Wtk=1d0/Lk
     !
     call dmft_gloc_matsubara_superc(Hk,Wtk,Gmats,Smats,iprint=4)
     call dmft_gloc_realaxis_superc(Hk,Wtk,Greal,Sreal,iprint=4)
  endif



contains





  subroutine solve_BdG(bool)
    real(8),dimension(Nlso) :: array
    integer                 :: niter
    real(8)                 :: fresult
    logical                 :: bool
    !
    array(1:Nlso)        = pii
    if(bool)print*,"Guess:",array
    call fmin_cg(array,get_FreeEnergy,niter,fresult,ftol=tol)
    pii = array(1:Nlso)
    if(bool)print*,"Result Phi :",pii
    !
  end subroutine solve_BdG



  subroutine search_mu_BdG(n,x,fvec,iflag)
    integer                  :: n
    real(8)                  :: x(n)
    real(8)                  :: fvec(n)
    integer                  :: iflag,iter=0
    real(8),dimension(Nlso)  :: dens
    save iter
    iter=iter+1
    xmu = x(1)
    call solve_BdG(bool=.false.)
    call get_Density(dens)
    fvec(1) = nread - 2*sum(dens)
    write(*,"(I5,4F16.9)")iter,x(1),fvec(1)
  end subroutine search_mu_BdG



  function get_FreeEnergy(a) result(Omega)
    real(8),dimension(:)                  :: a
    real(8)                               :: Omega
    complex(8),dimension(2*Nlso,2*Nlso)   :: Hknambu
    real(8),dimension(2*Nlso)             :: Eknambu
    complex(8)                            :: epsk
    real(8)                               :: Epk,Emk
    real(8)                               :: Delta(Nlso),dens(Nlso)
    integer                               :: i,ik
    !
    pii = a(1:Nlso)
    !
    Delta = -Uloc*pii
    !
    Omega=0d0
    do ik=1,Lk
       !
       Hknambu = hk_NambuHaldane_model(kvec(:,ik),2*Nlso)
       !
       epsk    = HkNambu(1,1) - Mh !get the diagonal dispersion and take out the Mass term!
       !
       call eigh(HkNambu,EkNambu)
       !
       Epk = EkNambu(1)
       Emk = EkNambu(2)
       !

       Omega = Omega - (abs(Epk) + abs(Emk) - 2*epsk)/dble(Lk)
    enddo
    !
    Omega = Omega + (Delta(1)**2 + Delta(2)**2)/Uloc !+ Uloc*(Dens(1)**2 + Dens(2)**2)
    !
    print*,Omega,Pii!,Dens
    !
  end function get_FreeEnergy









  subroutine get_LocalObservables()
    real(8),dimension(Nlso)             :: dens
    real(8),dimension(Nlso)             :: phii
    complex(8),dimension(2*Nlso,2*Nlso) :: Hknambu
    real(8),dimension(2*Nlso)           :: Eknambu
    real(8),dimension(2*Nlso)           :: RhoDiag  !diagonal density matrix 
    complex(8),dimension(2*Nlso,2*Nlso) :: RhoNambu !Nambu density matrix
    real(8),dimension(Nlso,Lk)          :: nk,phik
    integer                             :: i,ik,iorb
    !
    do ik=1,Lk
       !
       Hknambu = hk_NambuHaldane_model(kvec(:,ik),2*Nlso)
       !
       call eigh(HkNambu,EkNambu)
       !
       RhoDiag  = fermi(EkNambu,beta)
       RhoNambu = matmul(HkNambu, matmul(diag(RhoDiag),transpose(conjg(HkNambu)) ) )
       !
       forall(i=1:Nlso)
          nk(i,ik)  = dreal(RhoNambu(i,i))
          phik(i,ik) = dreal(RhoNambu(i,i+Nlso))
       end forall
       !
    enddo
    !
    dens = sum(nk,dim=2)/Lk
    phii = sum(phik,dim=2)/Lk
    !
    print*,"n_a  =",2*dens(1)
    print*,"n_b  =",2*dens(2)
    print*,"n_tot=",2*sum(dens)
    !
    open(10,file="density.bdg")
    write(10,"(10F20.12)")(2*dens(iorb),iorb=1,Nlso)
    close(10)
    open(10,file="phi.bdg")
    write(10,"(10F20.12)")(phii(iorb),iorb=1,Nlso)
    close(10)
    !
  end subroutine Get_LocalObservables


  subroutine get_Density(dens)
    real(8),dimension(Nlso)             :: dens
    complex(8),dimension(2*Nlso,2*Nlso) :: Hknambu
    real(8),dimension(2*Nlso)           :: Eknambu
    real(8),dimension(2*Nlso)           :: RhoDiag  !diagonal density matrix 
    complex(8),dimension(2*Nlso,2*Nlso) :: RhoNambu !Nambu density matrix
    real(8),dimension(Nlso,Lk)          :: nk,phik
    integer                             :: i,ik,iorb
    !
    do ik=1,Lk
       !
       Hknambu = hk_NambuHaldane_model(kvec(:,ik),2*Nlso)
       !
       call eigh(HkNambu,EkNambu)
       !
       RhoDiag  = fermi(EkNambu,beta)
       RhoNambu = matmul(HkNambu, matmul(diag(RhoDiag),transpose(conjg(HkNambu)) ) )
       !
       forall(i=1:Nlso)nk(i,ik)  = dreal(RhoNambu(i,i))
       !
    enddo
    !
    dens = sum(nk,dim=2)/Lk
    !
  end subroutine Get_Density





  !--------------------------------------------------------------------!
  !Haldane HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_haldane_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: h0,hx,hy,hz
    complex(8)                      :: h11,h22,h12,h21
    real(8)                         :: kdotd(3),kdota(3)
    !(k.d_j)
    kdotd(1) = dot_product(kpoint,d1)
    kdotd(2) = dot_product(kpoint,d2)
    kdotd(3) = dot_product(kpoint,d3)
    !(k.a_j)
    kdota(1) = dot_product(kpoint,a1)
    kdota(2) = dot_product(kpoint,a2)
    kdota(3) = dot_product(kpoint,a3)
    ! !
    !
    h11 = -2*tsp*sum( cos(kdota(:) + phi) ) + Mh - xmu
    h12 =    -ts*sum( exp(xi*kdotd(:)) )
    h21 =    -ts*sum( exp(-xi*kdotd(:)) )
    h22 = -2*tsp*sum( cos(kdota(:) - phi) ) - Mh - xmu
    !
    hk  = reshape([h11,h21,h12,h22],[2,2])
  end function hk_haldane_model

  function hk_NambuHaldane_model(kpoint,N) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: N
    complex(8),dimension(N,N)       :: hk
    complex(8),dimension(Nlso,Nlso) :: hk11,hk22
    real(8),dimension(Nlso)         :: sigma,self
    real(8)                         :: h0,hx,hy,hz
    if(2*Nlso/=N)stop "hk_NambuHaldane_model ERROR: 2*Nlso != N"
    !
    hk11  = hk_haldane_model(kpoint,Nlso)
    hk22  = -transpose(hk_haldane_model(-kpoint,Nlso))
    !
    sigma =  -Uloc*(nii(:)-0.5d0)
    self  =  -Uloc*pii(:)
    !
    Hk(1:Nlso,1:Nlso)               =  hk11 !+ diag(sigma)
    Hk(1:Nlso,Nlso+1:2*Nlso)        =       + diag(self)
    Hk(Nlso+1:2*Nlso,1:Nlso)        =       + diag(self)
    Hk(Nlso+1:2*Nlso,Nlso+1:2*Nlso) =  hk22 !- diag(sigma) !Hk22=-H^t(-k)
  end function hk_NambuHaldane_model




  subroutine build_EigenBands()
    real(8),dimension(4,2)                 :: kpath
    !
    KPath(1,:)=[0d0,0d0]
    KPath(2,:)=pointK
    Kpath(3,:)=pointKp
    KPath(4,:)=[0d0,0d0]
    call TB_Solve_model(hk_NambuHaldane_model,2*Nlso,KPath,Nkpath,&
         colors_name=[red1,blue1,green1,orange1],&
         points_name=[character(len=10) :: "G","K","K`","G"],&
         file="Eigenbands.bcs")
  end subroutine build_EigenBands






  ! !******************************************************************
  ! !******************************************************************





  ! calcola il numero di chern di un generico stato dipendente da k con il metodo di Resta
  subroutine Get_Chern_number(Hk,Nkvec,Noccupied,one_over_area,Chern)
    complex(8),intent(in),dimension(:,:,:)    :: Hk    ![Nlso][Nlso][Nktot]
    integer,intent(in),dimension(2)           :: Nkvec ![Nk1][Nk2]: prod(Nkvec)=Nktot
    integer,intent(in)                        :: Noccupied
    real(8),intent(in)                        :: one_over_area
    real(8),intent(out)                       :: Chern
    !
    integer                                   :: Nlso
    integer                                   :: Nktot
    integer                                   :: Nkx,Nky
    integer                                   :: ikx,iky
    integer                                   :: ikxP,ikyP
    integer                                   :: ik,iocc,i,j
    complex(8),dimension(:,:),allocatable     :: Eigvec ![Nlso][Nlso]
    real(8),dimension(:),allocatable          :: Eigval ![Nlso]
    complex(8),dimension(:,:),allocatable     :: Gmat
    complex(8),dimension(:,:,:,:),allocatable :: BlochStates ![Nkx][Nky][Noccupied][Nlso]
    complex(8),dimension(4)                   :: Ulink
    real(8),dimension(:,:),allocatable        :: BerryCurvature
    real(8)                                   :: berry_phase
    integer                                   :: unit
    !
    Nlso  = size(Hk,1)
    Nktot = size(Hk,3)
    Nkx   = Nkvec(1)
    Nky   = Nkvec(2)
    call assert_shape(Hk,[Nlso,Nlso,Nktot],"Get_Chern_NUmber_NEW","Hk")
    if(Nkx*Nky/=Nktot)stop "ERROR Get_Chern_Number_NEW: Nktot = prod(Nkvec)"
    !
    !
    !1. Get the Bloch states from H(:,:,k)
    allocate(Eigvec(Nlso,Nlso))
    allocate(Eigval(Nlso))
    allocate(BlochStates(Nkx,Nky,Noccupied,Nlso))
    allocate(BerryCurvature(Nkx,Nky))
    allocate(Gmat(Noccupied,Noccupied))
    ik=0
    do ikx=1,Nkx
       do iky=1,Nky
          ik=ik+1
          Eigvec = Hk(:,:,ik)
          call eigh(Eigvec,Eigval)
          do iocc=1,Noccupied
             BlochStates(ikx,iky,iocc,:) = Eigvec(:,iocc)
          enddo
       enddo
    enddo
    deallocate(Eigvec,Eigval)
    !
    !
    !2. Evaluate the Berry Curvature
    chern=0d0
    do ikx= 1, Nkx
       ikxP = modulo(ikx,Nkx) + 1
       !ikxM = modulo(ikx-2,Nkx) + 1
       do iky= 1, Nky
          ikyP = modulo(iky,Nky) + 1
          !ikyM = modulo(iky-2,Nky) + 1
          !
          if(Noccupied==1)then
             Ulink(1) = dot_product(BlochStates(ikx,iky,1,:)  , BlochStates(ikx,ikyP,1,:))
             Ulink(2) = dot_product(BlochStates(ikx,ikyP,1,:) , BlochStates(ikxP,ikyP,1,:))
             Ulink(3) = dot_product(BlochStates(ikxP,ikyP,1,:), BlochStates(ikxP,iky,1,:))
             Ulink(4) = dot_product(BlochStates(ikxP,iky,1,:) , BlochStates(ikx,iky,1,:))
             !
          else
             !
             forall(i=1:Noccupied,j=1:Noccupied)&
                  gmat(i,j)=dot_product(BlochStates(ikx,iky,i,:)  , BlochStates(ikx,ikyP,j,:))
             Ulink(1) = det(gmat)
             !
             forall(i=1:Noccupied,j=1:Noccupied)&
                  gmat(i,j) = dot_product(BlochStates(ikx,ikyP,i,:) , BlochStates(ikxP,ikyP,j,:))
             Ulink(2) = det(gmat)
             !
             forall(i=1:Noccupied,j=1:Noccupied)&
                  gmat(i,j) = dot_product(BlochStates(ikxP,ikyP,i,:), BlochStates(ikxP,iky,j,:))
             Ulink(3) = det(gmat)
             !
             forall(i=1:Noccupied,j=1:Noccupied)&
                  gmat(i,j) = dot_product(BlochStates(ikxP,iky,i,:) , BlochStates(ikx,iky,j,:))
             Ulink(4) = det(gmat)
             !
          endif
          !
          berry_phase = -dimag(zlog( product(Ulink(:))  ))
          chern = chern + berry_phase
          BerryCurvature(ikx,iky) = berry_phase*one_over_area
          !
       enddo
    enddo
    !
    chern=chern/pi2/3
    !
    open(unit=free_unit(unit),file="Chern_Number.bdg")
    write(unit,*)chern
    write(*,*)"C    =",chern
    close(unit)
    !
    call splot3d("Berry_Curvature.bdg",&
         linspace(0d0,pi2,Nkx,iend=.false.),&
         linspace(0d0,pi2,Nky,iend=.false.),&
         BerryCurvature(:Nkx,:Nky))
    !
  end subroutine Get_Chern_number










end program bdg_haldane
