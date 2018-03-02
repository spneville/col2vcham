module transmod

use constants

  implicit none

  save

  ! I/O
  integer                                     :: freqtypA,freqtypB
  character(len=120)                          :: freqfileA,freqfileB

  ! System information
  real(d), dimension(:), allocatable          :: xcoo0A,xcoo0B

  ! Normal modes
  real(d), dimension(:,:), allocatable        :: nmcooA,coonmA,&
                                                 nmcooB,coonmB
  real(d), dimension(:), allocatable          :: freqA,freqB
  character(len=3), dimension(:), allocatable :: nmlabA,nmlabB

  ! Normal mode shift
  real(d), dimension(:), allocatable          :: DeltaQ
  
  ! Rotation martix
  real(d), dimension(:,:), allocatable        :: Smatrix

  ! Transformed zeroth-order potential
  real(d), dimension(:), allocatable          :: kappaB
  real(d), dimension(:,:), allocatable        :: gammaB
  
end module transmod
