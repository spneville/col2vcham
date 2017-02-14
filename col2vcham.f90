program col2vcham

  use constants
  use global
  use ioqc

  implicit none

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
  call initialise

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdinp

!----------------------------------------------------------------------
! Determine the system dimensions
!----------------------------------------------------------------------
  call getdim

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  call alloc

!----------------------------------------------------------------------
! Read the Cartesian gradients and NACTs
!----------------------------------------------------------------------
  call rdgrad
  call rdnact

!----------------------------------------------------------------------
! Read the frequency file: normal modes and reference geometry
!----------------------------------------------------------------------
  call rdfreqfile

!----------------------------------------------------------------------
! Create the transformation matrices
!----------------------------------------------------------------------
  call nm2xmat

!----------------------------------------------------------------------
! Finalisation and deallocataion
!----------------------------------------------------------------------
  call finalise

contains

!######################################################################

  subroutine initialise
    
    use channels
    use iomod

    implicit none

!----------------------------------------------------------------------
! Open the logfile
!----------------------------------------------------------------------
    call freeunit(ilog)
    open(ilog,file='col2vcham.log',form='formatted',status='unknown')

    return

  end subroutine initialise

!######################################################################

  subroutine rdinp

    use constants
    use global
    use iomod
    use channels

    implicit none

    integer            :: i,n
    character(len=120) :: string1,string2

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    freqfile=''
    coldir=''

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
    n=0
5   continue
    
    n=n+1
    call getarg(n,string1)
       
    if (string1.eq.'-d') then
       n=n+1
       call getarg(n,coldir)
    else if (string1.eq.'-f') then
       n=n+1
       call getarg(n,freqfile)
    else
       write(6,'(/,2(2x,a),/)') 'Unknown keyword:',trim(string1)
       stop
    endif

    if (n.lt.iargc()) goto 5

!----------------------------------------------------------------------
! Make sure that all the required information has been given
!----------------------------------------------------------------------
    if (freqfile.eq.'') then
       errmsg="The frequency file has not been given"
       call error_control
    endif

    if (coldir.eq.'') then
       errmsg="The Columbus directory has not been given"
       call error_control
    endif

!----------------------------------------------------------------------
! Output some information to the log file
!----------------------------------------------------------------------
    write(ilog,'(/,2(2x,a))') "Frequency file:",trim(freqfile)
    write(ilog,'(/,2(2x,a))') "Columbus directory:",trim(coldir)

    return
 
  end subroutine rdinp

!######################################################################

  subroutine alloc

    use global

    implicit none

    ! Energy gradients
    allocate(grad(ncoo,nsta))
    grad=0.0d0

    ! NACTs
    allocate(nact(ncoo,nsta,nsta))
    grad=0.0d0

    ! Atomic masses
    allocate(mass(ncoo))
    mass=0.0d0

    ! Reference coordinates
    allocate(xcoo0(ncoo))
    xcoo0=0.0d0

    ! Atomic numbers
    allocate(atnum(natm))
    atnum=0

    ! Atom labels
    allocate(atlbl(natm))
    atlbl=''

    ! Normal mode labels
    allocate(nmlab(nmodes))
    nmlab=''

    ! Normal mode frequencies
    allocate(freq(nmodes))
    freq=0.0d0

    ! Transformation matrices
    allocate(nmcoo(ncoo,nmodes))
    allocate(coonm(nmodes,ncoo))
    nmcoo=0.0d0
    coonm=0.0d0

    return

  end subroutine alloc

!######################################################################

  subroutine finalise

    use channels
    use iomod

    implicit none

!----------------------------------------------------------------------
! Close the logfile
!----------------------------------------------------------------------
    close(ilog)

!----------------------------------------------------------------------
! Deallocation
!----------------------------------------------------------------------
    deallocate(grad)
    deallocate(nact)
    deallocate(mass)
    deallocate(xcoo0)
    deallocate(atnum)
    deallocate(atlbl)
    deallocate(nmlab)
    deallocate(freq)

    return
    
  end subroutine finalise

!######################################################################

end program col2vcham
