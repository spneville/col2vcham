module global
  
  use constants

  implicit none

  save

  ! I/O
  integer                                       :: freqtyp,qctyp,&
                                                   nlambdafiles,idrt
  character(len=120)                            :: freqfile
  character(len=120), dimension(50)             :: qcfile
  character(len=250), dimension(200)            :: lambdafile
  logical                                       :: outau,inpfile

  ! System information
  integer                                       :: natm,ncoo,nmodes,&
                                                   nsta,smax
  integer, dimension(:), allocatable            :: atnum
  real(d), dimension(:), allocatable            :: mass,xcoo0
  character(len=2), dimension(:), allocatable   :: atlbl
  logical                                       :: ldeuterate

  ! Energies, gradients, NACTs and 
  ! dipole matrix elements
  real(d), dimension(:), allocatable            :: ener
  real(d), dimension(:,:), allocatable          :: grad,kappa
  real(d), dimension(:,:,:), allocatable        :: nact,lambda
  real(d), dimension(:,:,:), allocatable        :: dipole
  logical                                       :: lgrad,lnact,&
                                                   laprxlambda

  ! Normal modes
  real(d), dimension(:,:), allocatable          :: nmcoo,coonm
  real(d), dimension(:), allocatable            :: freq
  character(len=3), dimension(:), allocatable   :: nmlab

  ! H_ML (Interaction with an external field)
  real(d)                                       :: omega,t0,sigma,I0
  logical                                       :: hml

end module global
