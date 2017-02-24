module global
  
  use constants

  implicit none

  save

  ! I/O
  integer                                     :: ityp
  character(len=120)                          :: freqfile,coldir
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

  ! Normal modes
  real(d), dimension(:,:), allocatable        :: nmcoo,coonm
  real(d), dimension(:), allocatable          :: freq
  character(len=3), dimension(:), allocatable :: nmlab

  ! H_ML (Interaction with an external field)
  real(d)                                     :: omega,t0,sigma,I0
  logical                                     :: hml

end module global
