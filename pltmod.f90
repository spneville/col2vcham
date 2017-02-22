module pltmod

  use constants

  implicit none

  save

  integer                              :: mplt,npnts,si,sf,surftyp
  real(d)                              :: qi,qf,ei,ef
  real(d), dimension(:,:), allocatable :: surf
  logical                              :: leps

end module pltmod
