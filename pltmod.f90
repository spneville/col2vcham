module pltmod

  use constants

  implicit none

  save

  integer                              :: mplt,npnts
  real(d)                              :: qi,qf,ei,ef
  real(d), dimension(:,:), allocatable :: surf

end module pltmod
