module global
  
  use constants

  implicit none

  save

  ! I/O
  integer                                     :: ityp
  character(len=120)                          :: freqfile,coldir

  ! System information
  integer                                     :: natm,ncoo,nmodes,nsta
  integer, dimension(:), allocatable          :: atnum
  real(d), dimension(:), allocatable          :: mass,xcoo0
  character(len=2), dimension(:), allocatable :: atlbl

  ! Energies, gradients and NACTs
  real(d), dimension(:), allocatable          :: ener
  real(d), dimension(:,:), allocatable        :: grad,kappa
  real(d), dimension(:,:,:), allocatable      :: nact,lambda

  ! Normal modes
  real(d), dimension(:,:), allocatable        :: nmcoo,coonm
  real(d), dimension(:), allocatable          :: freq
  character(len=3), dimension(:), allocatable :: nmlab

end module global
