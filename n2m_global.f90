module n2m_global

  use constants

  implicit none

  save

  ! System dimensions
  integer              :: nmodes,nsta

  ! Expansion order
  integer              :: order
  
  ! Frequencies
  real(d), allocatable :: freq(:)

  ! Expansion coefficients
  real(d), allocatable :: const(:,:)
  real(d), allocatable :: kappa(:,:)
  real(d), allocatable :: lambda(:,:,:)
  real(d), allocatable :: gamma(:,:,:)
  real(d), allocatable :: mu(:,:,:,:)

  ! NADVIBS file
  character(len=120)   :: nadvibsfile
  
end module n2m_global
