program bdg_haldane
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter                           :: Norb=1,Nspin=1,Nlat=2
  integer                                     :: Nso,Nlso
  integer                                     :: Nloop
  integer                                     :: Lmats,Lreal
  integer                                     :: Nsuccess
  integer                                     :: Lk

  !variables for the model:
  integer                                     :: Nk,Nkpath
  real(8)                                     :: Uloc
  real(8)                                     :: ts,tsp,phi,delta,Mh,wmixing,xmu,eps
  real(8)                                     :: beta
  real(8)                                     :: Wdis
  real(8)                                     :: nsc
  real(8)                                     :: deltasc
  real(8)                                     :: wini,wfin
  real(8)                                     :: bdg_error


  complex(8),allocatable,dimension(:,:,:)     :: Hk
  complex(8),allocatable,dimension(:,:)       :: halHloc
  complex(8),allocatable,dimension(:,:,:,:,:) :: Hloc
  real(8),allocatable,dimension(:)            :: Wtk
  real(8)                                     :: a,a0,bklen

  real(8),dimension(2)                        :: d1,d2,d3
  real(8),dimension(2)                        :: a1,a2,a3
  real(8),dimension(2)                        :: bk1,bk2,pointK,pointKp


  real(8),dimension(:),allocatable            :: nii,pii ![Nlat]
  logical                                     :: converged,bool
  integer                                     :: i,is,iloop,nrandom,idum
  real(8),dimension(:),allocatable            :: wm,wr
  character(len=50)                           :: finput

  ! READ INPUT FILES !
  call parse_cmd_variable(finput,"FINPUT",default="inputBDG.conf")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(ts,"TS",finput,default=1d0)
  call parse_input_variable(tsp,"TSP",finput,default=1.d0/3/sqrt(3d0))
  call parse_input_variable(mh,"MH",finput,default=0d0)
  call parse_input_variable(phi,"PHI",finput,default=0d0)
  call parse_input_variable(uloc,"ULOC",finput,default=-1d0,comment="Values of the local interaction")
  call parse_input_variable(xmu,"XMU",finput,default=0d0)
  call parse_input_variable(eps,"EPS",finput,default=1d-3)
  call parse_input_variable(beta,"BETA",finput,default=1000d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(nsc,"NSC",finput,default=0.5d0,comment="Value of the initial density per spin (1/2=half-filling).")
  call parse_input_variable(deltasc,"DELTASC",finput,default=0.02d0,comment="Value of the SC symmetry breaking term.")
  call parse_input_variable(nloop,"NLOOP",finput,default=500,comment="Max number of iterations.")
  call parse_input_variable(Lmats,"LMATS",finput,default=2000,comment="Number of Matsubara frequencies.")
  call parse_input_variable(Lreal,"LREAL",finput,default=2000,comment="Number of real-axis frequencies.")
  call parse_input_variable(wini,"WINI",finput,default=-15.d0,comment="Smallest real-axis frequency")
  call parse_input_variable(wfin,"WFIN",finput,default=15.d0,comment="Largest real-axis frequency")
  call parse_input_variable(bdg_error,"BDG_ERROR",finput,default=0.00001d0,comment="Error threshold for the convergence")
  call parse_input_variable(nsuccess,"NSUCCESS",finput,default=1,comment="Number of successive iterations below threshold for convergence")
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
  !
  !nearest neighbor: A-->B, B-->A
  d1= a*[  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= a*[  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= a*[ -1d0     , 0d0           ]
  !
  !
  !next nearest-neighbor displacements: A-->A, B-->B \== \nu_1,\nu_2, \nu_3=\nu_1-\nu_2
  a1=a0*[ sqrt(3d0)/2d0, 1d0/2d0]
  a2=a0*[ sqrt(3d0)/2d0,-1d0/2d0]
  a3=a2-a1
  !
  !
  !RECIPROCAL LATTICE VECTORS:
  bklen=4d0*pi/3d0
  bk1=bklen*[ 1d0/2d0 ,  sqrt(3d0)/2d0 ]
  bk2=bklen*[ 1d0/2d0 , -sqrt(3d0)/2d0 ]
  !
  !
  pointK = [2*pi/3, 2*pi/3/sqrt(3d0)]
  pointKp= [2*pi/3,-2*pi/3/sqrt(3d0)]


  ! ALLOCATE ALL THE REQUIRED VARIABLES !
  allocate(nii(Nlso))
  allocate(pii(Nlso))

  ! ALLOCATE MATSUBARA AND REAL FREQUENCIES !
  allocate(wm(Lmats),wr(Lreal))
  wr = linspace(wini,wfin,Lreal)
  wm = pi/beta*(2*arange(1,Lmats)-1)



  ! GET THE TIGHT BINDING HAMILTONIAN FOR THE SQUARE LATTICE 
  call build_hk("Hkfile_Haldane_sc.bdg")


  nii=nsc
  pii=deltasc
  call read_data("nVSisite.bdg",nii)
  call read_data("phiVSisite.bdg",pii)
  iloop=0;
  converged=.false.;
  !+-------------------------------------+!
  do while(.not.converged.AND.iloop<nloop) 
     iloop=iloop+1
     call start_loop(iloop,nloop,"BdG-loop")

     call BdG_Solve(nii,pii,Hk)
     converged = check_convergence_local(pii,bdg_error,nsuccess,nloop,file="error.err")
     call print_sc_out(converged)
     call end_loop()
  enddo
  !+-------------------------------------+!

  call ChernNumber()

contains



  subroutine BdG_Solve(nii,pii,H)
    real(8),dimension(Nlso),intent(inout)         :: nii
    real(8),dimension(Nlso),intent(inout)         :: pii
    complex(8),dimension(Nlso,Nlso,Lk),intent(in) :: H
    real(8),dimension(Nlso)                       :: Sigma_HFB
    real(8),dimension(Nlso)                       :: Self_HFB
    real(8),dimension(2*Nlso,2*Nlso)              :: Hknambu
    real(8),dimension(2*Nlso)                     :: Eknambu
    real(8),dimension(2*Nlso)                     :: RhoDiag  !diagonal density matrix 
    real(8),dimension(2*Nlso,2*Nlso)              :: RhoNambu !Nambu density matrix
    real(8),dimension(Nlso,Lk)                    :: nk,phik
    integer                                       :: i,ik
    
    !
    Sigma_HFB(:) =  Uloc*(nii-0.5d0)
    Self_HFB(:)  =  Uloc*pii(:)
    do ik=1,Lk
       !
       Hknambu(1:Nlso,1:Nlso)               =  H(:,:,ik) + diag(Sigma_HFB)
       Hknambu(1:Nlso,Nlso+1:2*Nlso)        =    + diag(Self_HFB)
       Hknambu(Nlso+1:2*Nlso,1:Nlso)        =    + diag(Self_HFB)
       Hknambu(Nlso+1:2*Nlso,Nlso+1:2*Nlso) = -conjg(H(:,:,ik) + diag(Sigma_HFB))
       !
       call eigh(HkNambu,EkNambu)
       !
       RhoDiag  = fermi(EkNambu,beta)
       RhoNambu = matmul(HkNambu, matmul(diag(RhoDiag),transpose(HkNambu)) )
       !
       forall(i=1:Nlso)
          nk(i,ik)  = RhoNambu(i,i)
          phik(i,ik) = RhoNambu(i,i+Nlso)
       end forall
       !
    enddo
    !
    nii = sum(nk,dim=2)/Lk
    pii = sum(phik,dim=2)/Lk
    !
  end subroutine BdG_Solve




  !--------------------------------------------------------------------!
  !Haldane HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_haldane_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: h0,hx,hy,hz
    real(8)                         :: kdotd(3),kdota(3)
    !(k.d_j)
    kdotd(1) = dot_product(kpoint,d1)
    kdotd(2) = dot_product(kpoint,d2)
    kdotd(3) = dot_product(kpoint,d3)
    !(k.a_j)
    kdota(1) = dot_product(kpoint,a1)
    kdota(2) = dot_product(kpoint,a2)
    kdota(3) = dot_product(kpoint,a3)
    !
    h0 = -2*tsp*cos(-phi)*sum( cos(kdota(:)) )
    hx =-ts*sum( cos(kdotd(:)) )
    hy =-ts*sum( sin(kdotd(:)) )
    hz = -2*tsp*sin(-phi)*sum( sin(kdota(:)) ) + Mh 
    hk = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z
  end function hk_haldane_model






  !---------------------------------------------------------------------
  !PURPOSE: Get Haldane Model Hamiltonian
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional             :: file
    integer                               :: i,j,ik
    integer                               :: ix,iy
    real(8)                               :: kx,ky    
    integer                               :: iorb,jorb
    integer                               :: isporb,jsporb
    integer                               :: ispin,jspin
    real(8),dimension(:),allocatable      :: kxgrid,kygrid
    real(8),dimension(2)                  :: kvec
    real(8),dimension(:,:),allocatable    :: kpath
    !
    write(*,*)"Build H(k) for Haldane model (now also Nobel Laureate):"
    Lk=Nk*Nk
    write(*,*)"# of k-points     :",Lk
    write(*,*)"# of SO-bands     :",Nlso
    !
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    !
    allocate(Hk(Nlso,Nlso,Lk))
    allocate(Wtk(Lk))
    print*,"Build Hk(Nlso,Nlso) for the Haldane model"
    ik=0
    do iy=1,Nk
       ky = dble(iy-1)/Nk
       do ix=1,Nk
          ik=ik+1
          kx=dble(ix-1)/Nk
          kvec = kx*bk1 + ky*bk2
          Hk(:,:,ik) = hk_haldane_model(kvec,Nlso)
       enddo
    enddo
    Wtk = 1d0/Lk
    !
    !
  end subroutine build_hk




  subroutine build_EigenBands()
    integer                               :: i,j,ik
    integer                               :: ix,iy
    real(8)                               :: kx,ky    
    integer                               :: iorb,jorb
    integer                               :: isporb,jsporb
    integer                               :: ispin,jspin
    real(8)                               :: foo,n0(Nlso)
    integer                               :: unit
    real(8),dimension(:),allocatable      :: kxgrid,kygrid
    real(8),dimension(2)                  :: kvec
    real(8),dimension(:,:),allocatable    :: kpath
    !
    allocate(Kpath(4,2))
    KPath(1,:)=[0d0,0d0]
    KPath(2,:)=pointK
    Kpath(3,:)=pointKp
    KPath(4,:)=[0d0,0d0]
    call TB_Solve_path(hk_haldane_model,Nlso,KPath,Nkpath,&
         colors_name=[red1,blue1],&
         points_name=[character(len=10) :: "G","K","K`","G"],&
         file="Eigenbands.bdg")
  end subroutine build_EigenBands













  ! !******************************************************************
  ! !******************************************************************






  subroutine ChernNumber
    complex(8),dimension(:,:,:),allocatable :: BlochStates
    real(8),dimension(:,:),allocatable      :: Berry_curvature
    integer                                 :: ik,ix,iy
    complex(8)                              :: Eigvec(Nlso,Nlso)
    real(8)                                 :: Chern,BZ_area,eigval(Nlso)
    real(8),dimension(:),allocatable        :: kxgrid
    !CHERN NUMBERS:
    allocate(BlochStates(Nlso,Nk,Nk))
    allocate(Berry_curvature(Nk,Nk))
    chern=0d0
    ik=0
    do ix=1,Nk
       do iy=1,Nk
          ik=ik+1
          Eigvec = Hk(:,:,ik) !-  nnn2lso_reshape(S0,Nlat,Nspin,Norb)
          call eigh(Eigvec,Eigval)
          BlochStates(:,ix,iy) = Eigvec(:,1)
       enddo
    enddo
    call get_Chern_number(BlochStates,Chern,Berry_curvature,Nk/pi2*Nk/pi2)
    write(*,*)"Chern =",chern
    allocate(kxgrid(Nk))
    kxgrid=linspace(0d0,pi2,Nk)
    call splot3d("Berry_Curvature.nint",kxgrid,kxgrid,Berry_Curvature)
    open(10,file="ChernNumber.dat")
    write(10,"(I3,F16.12)")nint(chern),chern
    close(10)
  end subroutine ChernNumber
  !
  subroutine Get_Chern_number(State,Chern_number,Berry_curvatures,one_over_area)
    complex(8),intent(in),dimension(:,:,:)                     :: state !(nhilb, nk1, nk2)
    real(8),intent(out)                                        :: chern_number
    real(8),intent(out),dimension(size(state,2),size(state,3)) :: berry_curvatures
    real(8),intent(in)                                         :: one_over_area
    integer                                                    :: nhilb,nk1,nk2
    integer                                                    :: i1,i2,i3,ix,i1p,i1m,i2m,i2p,it
    complex(8)                                                 :: path_scalar_products(4)
    real(8)                                                    :: berry_phase
    !
    Nhilb= size(state,1)
    Nk1  = size(state,2)
    Nk2  = size(state,3)
    chern_number = zero
    do i1= 1, nk1
       i1p = modulo(i1,nk1) + 1
       i1m = modulo(i1-2,nk1) + 1
       do i2= 1, nk2           !faccio l'integrale sulla bz
          i2p = modulo(i2,nk2) + 1
          i2m = modulo(i2-2,nk2) + 1
          path_scalar_products(1) = dot_product(state(:,i1,i2),state(:,i1, i2p))
          path_scalar_products(2) = dot_product(state(:,i1,i2p),state(:,i1p, i2p))
          path_scalar_products(3) = dot_product(state(:,i1p,i2p),state(:,i1p, i2))
          path_scalar_products(4) = dot_product(state(:,i1p,i2),state(:,i1,i2))
          berry_phase = -dimag(zlog( product(path_scalar_products)  ))
          berry_curvatures(i1,i2) = berry_phase*one_over_area
          chern_number = chern_number + berry_phase
       enddo
    enddo
    chern_number = chern_number/pi2
  end subroutine Get_Chern_number





  subroutine print_sc_out(converged)
    integer                               :: iorb,i,j,is,row,col,unit
    logical                               :: converged
    real(8)                               :: foo,n0(Nlso),phi0(Nlso)
    complex(8),dimension(2,Nlso,Nlso,Lmats) :: Gmats,fooSmats
    complex(8),dimension(2,Nlso,Nlso,Lreal) :: Greal,fooSreal
    !


    print*,"<nimp>  =",2*nii
    print*,"<phi>   =",pii

    call splot("nVSiloop.bdg",iloop,2*nii(1),2*nii(2),append=.true.)
    call splot("phiVSiloop.bdg",iloop,pii(1),pii(2),append=.true.)


    !Save the guess:
    call store_data("nVSisite.bdg",nii(:))
    call store_data("phiVSisite.bdg",pii(:))

    if(converged)then
       !Build the local GF:

       open(10,file="density.bdg")
       write(10,"(10F20.12)")(2*nii(iorb),iorb=1,Nlso)
       close(10)
       open(10,file="phi.bdg")
       write(10,"(10F20.12)")(pii(iorb),iorb=1,Nlso)
       close(10)

       do i=1,Lmats
          fooSmats(1,:,:,i)=diag(Uloc*(nii-0.5d0))
          fooSmats(2,:,:,i)=diag(Uloc*pii(:))
       enddo
       call dmft_gloc_matsubara_superc(Hk,Wtk,Gmats,fooSmats,iprint=4)

       do i=1,Lreal
          fooSreal(1,:,:,i)=diag(Uloc*(nii-0.5d0))
          fooSreal(2,:,:,i)=diag(Uloc*pii(:))
       enddo
       call dmft_gloc_realaxis_superc(Hk,Wtk,Greal,fooSreal,iprint=4)

       ! do iorb=1,Nlso
       !    n0(iorb) = fft_get_density(Gmats(1,iorb,iorb,:),beta)
       !    phi0(iorb) = fft_get_density(Gmats(2,iorb,iorb,:),beta,notail=.true.)
       ! enddo
       ! open(10,file="density_from_GF.bdg")
       ! write(10,"(10F20.12)")(n0(iorb),iorb=1,Nlso)
       ! close(10)
       ! open(10,file="phi_from_GF.bdg")
       ! write(10,"(10F20.12)")(phi0(iorb),iorb=1,Nlso)
       ! close(10)
    endif

  end subroutine print_sc_out



end program bdg_haldane
