program pltlvc

  use constants

  implicit none

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
  call initialise

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdpltinp

!----------------------------------------------------------------------
! Get the system dimensions
!----------------------------------------------------------------------
  call getdim

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  call alloc

!----------------------------------------------------------------------
! Read the LVC parmaters from the lvc.dat file
!----------------------------------------------------------------------
  call rddatfile

!----------------------------------------------------------------------
! Calculate the surfaces
!----------------------------------------------------------------------
  call calcsurf

!----------------------------------------------------------------------
! Write the gnuplot file and plot the surfaces to the screen
!----------------------------------------------------------------------
  call wrgnuplot

!----------------------------------------------------------------------
! Finalisation
!----------------------------------------------------------------------
  call finalise

contains

!######################################################################

  subroutine initialise

    use channels
    use iomod

!----------------------------------------------------------------------
! Open the logfile
!----------------------------------------------------------------------
    call freeunit(ilog)
    open(ilog,file='pltlvc.log',form='formatted',status='unknown')

    return

  end subroutine initialise

!######################################################################

  subroutine rdpltinp

    use pltmod
    use iomod

    implicit none

    integer :: n
    character(len=120) :: string1,string2

!----------------------------------------------------------------------
! Set default values
!----------------------------------------------------------------------
    ! Plotting mode
    mplt=-1

    ! Initial and final states
    si=-1
    sf=-1

    ! Coordinate interval
    qi=-7.0d0
    qf=+7.0d0

    ! Energy interval
    ei=-999.0d0
    ef=-999.0d0

    ! No. points
    npnts=1000

    ! eps output
    leps=.false.

    ! Surface type: 1 <-> adiabatic potentials
    !               2 <-> diabatic potentials
    surftyp=1

!----------------------------------------------------------------------
! Exit if no arguments have been given
!----------------------------------------------------------------------
    if (iargc().eq.0) then
       errmsg='No arguments have been given!'
       call error_control
    endif

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
    n=0
5   continue
    
    n=n+1
    call getarg(n,string1)

    if (string1.eq.'-h') then
       call wrhelp
    else if (string1.eq.'-m') then
       n=n+1
       call getarg(n,string2)
       read(string2,*) mplt
    else if (string1.eq.'-xrange') then
       n=n+1
       call getarg(n,string2)
       read(string2,*) qi
       n=n+1
       call getarg(n,string2)
       read(string2,*) qf
    else if (string1.eq.'-yrange') then
       n=n+1
       call getarg(n,string2)
       read(string2,*) ei
       n=n+1
       call getarg(n,string2)
       read(string2,*) ef
    else if (string1.eq.'-npnts') then
       n=n+1
       call getarg(n,string2)
       read(string2,*) npnts
    else if (string1.eq.'-si') then
       n=n+1
       call getarg(n,string2)
       read(string2,*) si
    else if (string1.eq.'-sf') then
       n=n+1
       call getarg(n,string2)
       read(string2,*) sf
    else if (string1.eq.'-eps') then
       leps=.true.
    else if (string1.eq.'-adiab') then
       surftyp=1
    else if (string1.eq.'-diab') then
       surftyp=2
    else
       errmsg='Unknown keyword: '//trim(string1)
       call error_control
    endif

    if (n.lt.iargc()) goto 5

!----------------------------------------------------------------------
! Make sure that all the required information has been given
!----------------------------------------------------------------------
    if (mplt.eq.-1) then
       errmsg='The plotting mode has not been given'
       call error_control
    endif

    return

  end subroutine rdpltinp

!######################################################################

  subroutine wrhelp

    implicit none

    integer :: i

!----------------------------------------------------------------------
! Write the input options to screen, then quit
!----------------------------------------------------------------------
    ! Purpose
    write(6,'(/,25a)') ('-',i=1,25)
    write(6,'(a)') 'Purpose'
    write(6,'(25a)') ('-',i=1,25)
    write(6,'(a)') 'Plots the model potentials of the LVC &
         Hamiltonian constructed using the qc2vcham program'

    ! Usage
    write(6,'(/,25a)') ('-',i=1,25)
    write(6,'(a)') 'Usage'
    write(6,'(25a)') ('-',i=1,25)
    write(6,'(a)') 'pltlvc -m (-xrange -yrange -npnts -si -sf -eps &
         -diab -adiab)'

    ! Options
    write(6,'(/,25a)') ('-',i=1,25)
    write(6,'(a)')   'Options'
    write(6,'(25a)') ('-',i=1,25)
    write(6,'(a)') '-m M              : &
         The plot is along the Mth mode'
    write(6,'(a)') '-xrange QI QF     : &
         The plot is over the coordinate interval [QI,QF]'
    write(6,'(a)') '-yrange EI EF     : &
         The plot is over the energy interval [EI,EF]'
    write(6,'(a)') '-npnts N          : &
         The plot will be made using N points for each state'
    write(6,'(a)') '-si I             : &
         I is the index of the lowest state potential to be plotted'
    write(6,'(a)') '-sf F             : &
         F is the index of the highest state potential to be plotted'
    write(6,'(a)') '-eps              : &
         Also write the model potentials to an eps file'
    write(6,'(a)') '-diab             : &
         Plot the diabatic potentials (default)'
    write(6,'(a)') '-adiab            : &
         Plot the adiabatic potentials'
    
    write(6,'(/)')
    STOP
    
    return
    
  end subroutine wrhelp
    
!######################################################################
  
  subroutine rddatfile
    
    use constants
    use iomod
    use global

    implicit none

    integer :: unit,m,s,s1,s2

!----------------------------------------------------------------------
! Open the lvc.dat file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file='lvc.dat',form='unformatted',status='old')

!----------------------------------------------------------------------
! Write the lvc.dat file
!----------------------------------------------------------------------
    ! Dimensions
    read(unit) nmodes
    read(unit) nsta

    ! Frequencies
    do m=1,nmodes
       read(unit) freq(m)
    enddo

    ! Energies
    do s=1,nsta
       read(unit) ener(s)
    enddo

    ! 1st-order intrastate coupling terms (kappa)
    do m=1,nmodes
       do s=1,nsta
          read(unit) kappa(m,s)
       enddo
    enddo
    
    ! 1st-order interstate coupling terms (lambda)
    do m=1,nmodes
       do s1=1,nsta-1
          do s2=s1+1,nsta
             read(unit) lambda(m,s1,s2)
             lambda(m,s2,s1)=lambda(m,s1,s2)
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Close the lvc.dat file
!----------------------------------------------------------------------
    close(unit)    

    return

  end subroutine rddatfile

!######################################################################

  subroutine getdim

    use constants
    use iomod
    use global

    implicit none
    
    integer :: unit

!----------------------------------------------------------------------
! Open the lvc.dat file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file='lvc.dat',form='unformatted',status='old')

!----------------------------------------------------------------------
! Read the system dimensions
!----------------------------------------------------------------------
    read(unit) nmodes
    read(unit) nsta

!----------------------------------------------------------------------
! Close the lvc.dat file
!----------------------------------------------------------------------
    close(unit)

    return

  end subroutine getdim

!######################################################################

  subroutine alloc

    use global
    use pltmod

    implicit none

    ! Surface to be plotted
    allocate(surf(npnts,nsta))
    surf=0.0d0

    ! Frequencies
    allocate(freq(nmodes))
    freq=0.0d0

    ! Energies
    allocate(ener(nsta))
    ener=0.0d0
    
    ! kappa
    allocate(kappa(nmodes,nsta))
    kappa=0.0d0
    
    ! lambda
    allocate(lambda(nmodes,nsta,nsta))
    lambda=0.0d0

    return

  end subroutine alloc

!######################################################################

  subroutine calcsurf

    use constants
    use iomod
    use global
    use pltmod

    implicit none

    integer                    :: unit,i,j,s
    real(d)                    :: dq
    real(d), dimension(nmodes) :: q
    character(len=2)           :: am,as
    character(len=80)          :: datfile,filename,string

!----------------------------------------------------------------------
! Calculate the LVC adiabatic potentials
!----------------------------------------------------------------------
    ! Step size
    dq=(qf-qi)/(npnts-1) 

    q=0.0d0

    ! Loop over points
    do i=1,npnts
       
       ! Current geometry
       q(mplt)=qi+(i-1)*dq

       ! Current set of energies
       surf(i,:)=surface(q)

    enddo

!----------------------------------------------------------------------
! Open the data file
!----------------------------------------------------------------------
    write(am,'(i2)') mplt
    datfile='plot_q'//trim(adjustl(am))//'.dat'
    
    call freeunit(unit)
    open(unit,file=datfile,form='formatted',status='unknown')

!----------------------------------------------------------------------
! Write the data file
!----------------------------------------------------------------------
    do i=1,npnts
       write(unit,*) qi+(i-1)*dq,(surf(i,j),j=1,nsta)
    enddo

!----------------------------------------------------------------------
! Close the data file
!----------------------------------------------------------------------
    close(unit)

    return

  end subroutine calcsurf

!######################################################################

  function surface(q) result(func)

    use constants
    use global
    use pltmod

    implicit none
    
    real(d), dimension(nmodes) :: q
    real(d), dimension(nsta)   :: func

    select case(surftyp)

    case(1) ! adiabatic potentials
       func=adiabpot(q)

    case(2) ! diabatic potentials
       func=diabpot(q)

    end select

    return

  end function surface

!######################################################################

  function adiabpot(q) result(v)

    use constants
    use iomod
    use global

    implicit none

    integer                       :: e2,error
    real(d), dimension(nmodes)    :: q
    real(d), dimension(nsta)      :: v
    real(d), dimension(nsta,nsta) :: w
    real(d), dimension(3*nsta)    :: work

!----------------------------------------------------------------------
! Construct the LVC potential
!----------------------------------------------------------------------
    w=lvcpot(q)

!----------------------------------------------------------------------
! Diagonalise the LVC potential to yield the model adiabatic
! potentials
!----------------------------------------------------------------------
    e2=3*nsta
    call dsyev('V','U',nsta,w,nsta,v,work,e2,error)

    if (error.ne.0) then
       errmsg='Diagonalisation of the LVC potential failed'
       call error_control
    endif

    return

  end function adiabpot

!######################################################################

  function diabpot(q) result(wii)

    use constants
    use global

    implicit none

    integer                       :: i
    real(d), dimension(nmodes)    :: q
    real(d), dimension(nsta)      :: wii
    real(d), dimension(nsta,nsta) :: w

!----------------------------------------------------------------------
! Construct the LVC potential
!----------------------------------------------------------------------
    w=lvcpot(q)

!----------------------------------------------------------------------
! Return the on-diagonal elements
!----------------------------------------------------------------------
    do i=1,nsta
       wii(i)=w(i,i)
    enddo

    return

  end function diabpot

!######################################################################

  function lvcpot(q) result(w)

    use constants
    use global

    implicit none

    integer                       :: m,s,s1,s2
    real(d), dimension(nmodes)    :: q
    real(d), dimension(nsta,nsta) :: w

!----------------------------------------------------------------------
! Initialisation of the LVC potential
!----------------------------------------------------------------------
    w=0.0d0

!----------------------------------------------------------------------
! Zeroth-order contributions
!----------------------------------------------------------------------
    ! Energies
    do s=1,nsta
       w(s,s)=w(s,s)+ener(s)
    enddo

    ! Harmonic potentials
    do s=1,nsta
       do m=1,nmodes
          w(s,s)=w(s,s)+0.5d0*freq(m)*q(m)**2
       enddo
    enddo

!----------------------------------------------------------------------
! First-order contributions
!----------------------------------------------------------------------
    ! kappa
    do s=1,nsta
       do m=1,nmodes
          w(s,s)=w(s,s)+kappa(m,s)*q(m)
       enddo
    enddo

    ! lambda
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m=1,nmodes
             w(s1,s2)=w(s1,s2)+lambda(m,s1,s2)*q(m)
             w(s2,s1)=w(s2,s1)+lambda(m,s1,s2)*q(m)
          enddo
       enddo
    enddo

    return

  end function lvcpot

!######################################################################

  subroutine wrgnuplot

    use constants
    use iomod
    use global
    use pltmod

    implicit none

    integer           :: unit,s
    character(len=2)  :: am,as
    character(len=80) :: filename,datfile,string

!----------------------------------------------------------------------
! Filenames
!----------------------------------------------------------------------
    write(am,'(i2)') mplt

    ! data file
    datfile='plot_q'//trim(adjustl(am))//'.dat'

    ! gnuplot file
    filename='plot_q'//trim(adjustl(am))//'.gnu'

!----------------------------------------------------------------------
! Energy ranges and states
!----------------------------------------------------------------------
    if (si.eq.-1) si=1
    if (sf.eq.-1) sf=nsta
    if (ei.eq.-999.0d0) ei=0.98d0*minval(surf(:,si:sf))
    if (ef.eq.-999.0d0) ef=1.02d0*maxval(surf(:,si:sf))

!----------------------------------------------------------------------
! Open the gnuplot file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='unknown')

!----------------------------------------------------------------------
! Write the gnuplot file
!----------------------------------------------------------------------
    ! Set up
    write(unit,'(a)') 'set size square'
    write(unit,'(a)') 'unset key'
    !write(unit,'(a)') 'monitorSize=system("xrandr | awk &
    !     ''/\*/{sub(/x/,\",\");print $1;exit}''")'
    !write(unit,'(a,/)') 'set terminal x11 size @monitorSize'

    ! Axis labels
    write(unit,'(a)') 'set ylabel ''Energy (eV)'''
    write(am,'(i2)') mplt
    string='set xlabel ''Q_{'//trim(adjustl(am))//'}'''
    write(unit,'(a,/)') trim(string)

    ! Ranges
    write(unit,'(2(a,F6.2),a)') 'set xrange [',qi,':',qf,']'
    write(unit,'(2(a,F6.2),a,/)') 'set yrange [',ei,':',ef,']'

    ! State si
    string='plot '''//trim(datfile)//''' u 1:'
    write(as,'(i2)') si+1
    string=trim(string)//trim(adjustl(as))//' w l lw 4'
    write(unit,'(a)') trim(string)

    ! States si+1 to sf
    do s=si+1,sf
       string='replot '''//trim(datfile)//''' u 1:'
       write(as,'(i2)') s+1
       string=trim(string)//trim(adjustl(as))//' w l lw 4'
       write(unit,'(a)') trim(string)
    enddo

    ! eps output
    if (leps) then
       write(unit,'(/,a)') 'set terminal postscript eps &
            enhanced color solid "Helvetica" 20'
       write(unit,'(a)') 'set output '''&
            //filename(1:index(filename,'.gnu'))//'eps'''
       write(unit,'(a)') 'replot'
    endif

    ! pause -1
    write(unit,'(/,a)') 'pause -1'

!----------------------------------------------------------------------
! Close the gnuplot file
!----------------------------------------------------------------------
    close(unit)

!----------------------------------------------------------------------
! Plot the surfaces to the screen
!----------------------------------------------------------------------
    call system('gnuplot '//trim(filename))

    return

  end subroutine wrgnuplot

!######################################################################

  subroutine finalise

    use channels
    use global

    implicit none

    ! Close the logfile
    close(ilog)

    ! Deallocate arrays
    deallocate(freq)
    deallocate(ener)
    deallocate(kappa)
    deallocate(lambda)

    return

  end subroutine finalise

!######################################################################

end program pltlvc
