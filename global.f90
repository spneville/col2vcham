module global
  
  use constants

  implicit none

  save

  ! I/O
  integer                                     :: freqtyp,qctyp
  character(len=120)                          :: freqfile
  character(len=120), dimension(50)           :: qcfile
  logical                                     :: outau

  ! System information
  integer                                     :: natm,ncoo,nmodes,nsta
  integer, dimension(:), allocatable          :: atnum
  real(d), dimension(:), allocatable          :: mass,xcoo0
  character(len=2), dimension(:), allocatable :: atlbl

  ! Energies, gradients, NACTs and 
  ! dipole matrix elements
  real(d), dimension(:), allocatable          :: ener
  real(d), dimension(:,:), allocatable        :: grad,kappa
  real(d), dimension(:,:,:), allocatable      :: nact,lambda
  real(d), dimension(:,:,:), allocatable      :: dipole
  logical                                     :: lgrad,lnact

  ! Normal modes
  real(d), dimension(:,:), allocatable        :: nmcoo,coonm
  real(d), dimension(:), allocatable          :: freq
  character(len=3), dimension(:), allocatable :: nmlab

  ! H_ML (Interaction with an external field)
  real(d)                                     :: omega,t0,sigma,I0
  logical                                     :: hml

end module global
