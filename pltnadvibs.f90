program pltnadvibs

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
! Allocate arrays
!----------------------------------------------------------------------
  call alloc

!----------------------------------------------------------------------
! Read the NADVIBS input
!----------------------------------------------------------------------
  call rdnadvibsinp

!----------------------------------------------------------------------
! Re-order the parameters s.t. the mode indices are ordered in terms
! of increasing frequency
!----------------------------------------------------------------------
  call reorderpar

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
    open(ilog,file='pltnadvibs.log',form='formatted',status='unknown')

    return

  end subroutine initialise
  
!######################################################################

    subroutine rdpltinp

    use pltmod
    use n2m_global
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

    ! Number of modes
    nmodes=0

    ! Number of states
    nsta=0

    ! NADVIBS input file name
    nadvibsfile=''

    ! Expansion order
    order=0

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
    if (iargc().ne.0) then

       n=0

5      continue
       n=n+1
       call getarg(n,string1)
       
       if (string1.eq.'-ns') then
          n=n+1
          call getarg(n,string2)
          read(string2,*) nsta
          
       else if (string1.eq.'-nm') then
          n=n+1
          call getarg(n,string2)
          read(string2,*) nmodes
          
       else if (string1.eq.'-f') then
          n=n+1
          call getarg(n,nadvibsfile)
          
       else if (string1.eq.'-o') then
          n=n+1
          call getarg(n,string2)
          read(string2,*) order
          
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

    endif
       
!----------------------------------------------------------------------
! If arguments are missing, get them from the user
!----------------------------------------------------------------------
    ! Number of states
    if (nsta.eq.0) then
       write(6,'(a)') 'Please enter the number of states'
       read(5,*) nsta
    endif

    ! Number of modes
    if (nmodes.eq.0) then
       write(6,'(a)') 'Please enter the number of modes'
       read(5,*) nmodes
    endif

    ! NADVIBS input file name
    if (nadvibsfile.eq.'') then
       write(6,'(a)') 'Please enter the name of the NADVIBS input file'
       read(5,'(a)') nadvibsfile
    endif

    ! Expansion order
    if (order.eq.0) then
       write(6,'(a)') 'Please enter the expansion order'
       read(5,*) order
    endif

!----------------------------------------------------------------------
! Check that all the required information has been given
!----------------------------------------------------------------------
    ! Number of states
    if (nsta.eq.0) then
       errmsg='The number of states has not been given'
       call error_control
    endif

    ! Number of modes
    if (nmodes.eq.0) then
       errmsg='The number of modes has not been given'
       call error_control
    endif

    ! NADVIBS input file name
    if (nadvibsfile.eq.'') then
       errmsg='The name of the NADVIBS input file has not been given'
       call error_control
    endif

    ! Expansion order
    if (order.eq.0) then
       errmsg='The expansion order has not been given'
       call error_control
    endif

!----------------------------------------------------------------------
! Exit if the expansion order is not supported
!----------------------------------------------------------------------
    if (order.gt.2) then
       errmsg='Expansion orders > 2 are not yet supported'
       call error_control
    endif
       
!----------------------------------------------------------------------
! If arguments are missing, get them from the user
!----------------------------------------------------------------------
    ! Number of states
    if (nsta.eq.0) then
       write(6,'(a)') 'Please enter the number of states'
       read(5,*) nsta
    endif

    ! Number of modes
    if (nmodes.eq.0) then
       write(6,'(a)') 'Please enter the number of modes'
       read(5,*) nmodes
    endif

    ! NADVIBS input file name
    if (nadvibsfile.eq.'') then
       write(6,'(a)') 'Please enter the name of the NADVIBS input file'
       read(5,'(a)') nadvibsfile
    endif

    ! Expansion order
    if (order.eq.0) then
       write(6,'(a)') 'Please enter the expansion order'
       read(5,*) order
    endif

!----------------------------------------------------------------------
! Check that all the required information has been given
!----------------------------------------------------------------------
    ! Number of states
    if (nsta.eq.0) then
       errmsg='The number of states has not been given'
       call error_control
    endif

    ! Number of modes
    if (nmodes.eq.0) then
       errmsg='The number of modes has not been given'
       call error_control
    endif

    ! NADVIBS input file name
    if (nadvibsfile.eq.'') then
       errmsg='The name of the NADVIBS input file has not been given'
       call error_control
    endif

    ! Expansion order
    if (order.eq.0) then
       errmsg='The expansion order has not been given'
       call error_control
    endif

    ! Plotting mode
    if (mplt.eq.-1) then
       errmsg='The plotting mode has not been given'
       call error_control
    endif
    
!----------------------------------------------------------------------
! Exit if the expansion order is not supported
!----------------------------------------------------------------------
    if (order.gt.2) then
       errmsg='Expansion orders > 2 are not yet supported'
       call error_control
    endif

    return
    
  end subroutine rdpltinp
    
!######################################################################

  subroutine alloc

    use constants
    use iomod
    use n2m_global
    use pltmod
    
    implicit none

    ! Frequencies
    allocate(freq(nmodes))
    freq=0.0d0
    
    ! Expansion coefficients
    allocate(const(nsta,nsta))
    allocate(kappa(nmodes,nsta))
    allocate(lambda(nmodes,nsta,nsta))
    allocate(gamma(nmodes,nmodes,nsta))
    allocate(mu(nmodes,nmodes,nsta,nsta))
    const=0.0d0
    kappa=0.0d0
    lambda=0.0d0
    gamma=0.0d0
    mu=0.0d0

    ! Surface to be plotted
    allocate(surf(npnts,nsta))
    surf=0.0d0
    
    return
    
  end subroutine alloc
  
!######################################################################

    subroutine rdnadvibsinp

    use constants
    use iomod
    use n2m_global
    
    implicit none

    integer                                 :: unit,i,j,m,m1,m2,n
    real(d), dimension(nmodes)              :: tmp1
    real(d), dimension(nmodes*(nmodes+1)/2) :: tmp2
    
!----------------------------------------------------------------------
! Open the NADVIBS input file
!----------------------------------------------------------------------
    call freeunit(unit)
    open (unit,file=nadvibsfile,form='formatted',status='old')
    
!----------------------------------------------------------------------
! Read the frequencies
!----------------------------------------------------------------------
    read(unit,*)
    read(unit,*) (freq(m),m=1,nmodes)

!-----------------------------------------------------------------------
! Read the expansion coefficients
!-----------------------------------------------------------------------
    ! Loop over columns
    do j=1,nsta

       ! Loop over rows
       do i=1,j

          ! Constants
          read(unit,*)
          read(unit,*) const(i,j)
          const(j,i)=const(i,j)

          ! First-order coefficients
          if (order.ge.1) then
             read(unit,*)
             read(unit,*) (tmp1(m),m=1,nmodes)
             if (i.eq.j) then
                kappa(:,i)=tmp1(:)
             else
                lambda(:,i,j)=tmp1(:)
                lambda(:,j,i)=tmp1(:)
             endif
          endif

          ! Second-order coefficients
          if (order.ge.2) then
             read(unit,*)
             read(unit,*) (tmp2(m),m=1,nmodes*(nmodes+1)/2)
             n=0
             do m2=1,nmodes
                do m1=1,m2
                   n=n+1
                   if (i.eq.j) then
                      gamma(m1,m2,i)=tmp2(n)
                      gamma(m2,m1,i)=tmp2(n)
                   else
                      mu(m1,m2,i,j)=tmp2(n)
                      mu(m1,m2,j,i)=tmp2(n)
                      mu(m2,m1,i,j)=tmp2(n)
                      mu(m2,m1,j,i)=tmp2(n)
                   endif
                enddo
             enddo
          endif
          
       enddo

    enddo
    
!----------------------------------------------------------------------
! Close the NADVIBS input file
!----------------------------------------------------------------------
    close(unit)

!-----------------------------------------------------------------------
! Unit conversions: frequency-weight and convert to eV
!-----------------------------------------------------------------------
    ! Constants
    const=const*eh2ev

    ! First-order coefficients
    do m=1,nmodes
       kappa(m,:)=kappa(m,:)/sqrt(freq(m))
       lambda(m,:,:)=lambda(m,:,:)/sqrt(freq(m))
    enddo
    kappa=kappa*eh2ev
    lambda=lambda*eh2ev

    ! Second-order coefficients
    if (order.ge.2) then
       do m1=1,nmodes
          do m2=1,nmodes
             gamma(m1,m2,:)=gamma(m1,m2,:)&
                  /sqrt(freq(m1))/sqrt(freq(m2))
             mu(m1,m2,:,:)=mu(m1,m2,:,:)&
                  /sqrt(freq(m1))/sqrt(freq(m2))
          enddo
       enddo
    endif
    gamma=gamma*eh2ev
    mu=mu*eh2ev

    ! Frequencies: this need to be done last, else the frequencies
    ! used in the frequency-weighting have to be scaled by 1/ev2eh
    freq=freq*eh2ev

    return
    
  end subroutine rdnadvibsinp
  
!######################################################################

    subroutine reorderpar

    use constants
    use utils
    use n2m_global
    
    implicit none

    integer                                     :: m,m1,m2,s,s1,s2
    integer, dimension(nmodes)                  :: indx
    real(d), dimension(nmodes)                  :: freq1
    real(d), dimension(nmodes,nsta)             :: kappa1
    real(d), dimension(nmodes,nsta,nsta)        :: lambda1
    real(d), dimension(nmodes,nmodes,nsta)      :: gamma1
    real(d), dimension(nmodes,nmodes,nsta,nsta) :: mu1
    
!----------------------------------------------------------------------
! Get the list of mode indices in order of increasing frequency
!----------------------------------------------------------------------
    call dsortindxa1('A',nmodes,freq,indx)

!----------------------------------------------------------------------
! Reorder the parameters so that the mode indices are in order of
! increasing frequency
!----------------------------------------------------------------------
    ! Frequencies
    do m=1,nmodes
       freq1(m)=freq(indx(m))
    enddo
    freq=freq1

    ! 1st-order intrastate coupling constants
    do m=1,nmodes
       kappa1(m,:)=kappa(indx(m),:)
    enddo
    kappa=kappa1

    ! 1st-order interstate coupling constants
    do m=1,nmodes
       lambda1(m,:,:)=lambda(indx(m),:,:)
    enddo
    lambda=lambda1

    ! 2nd-order intrastate coupling constants
    do m1=1,nmodes
       do m2=1,nmodes
          gamma1(m1,m2,:)=gamma(indx(m1),indx(m2),:)
       enddo
    enddo
    gamma=gamma1

    ! 2nd-order interstate coupling constants
    do m1=1,nmodes
       do m2=1,nmodes
          mu1(m1,m2,:,:)=mu(indx(m1),indx(m2),:,:)
       enddo
    enddo
    mu=mu1
    
    return
    
  end subroutine reorderpar

!######################################################################

    subroutine calcsurf

    use constants
    use iomod
    use n2m_global
    use pltmod

    implicit none

    integer                    :: unit,i,j,s
    real(d)                    :: dq
    real(d), dimension(nmodes) :: q
    character(len=2)           :: am,as
    character(len=80)          :: datfile,filename,string

!----------------------------------------------------------------------
! Calculate the model potentials
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
    use n2m_global
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
    use n2m_global

    implicit none

    integer                       :: e2,error
    real(d), dimension(nmodes)    :: q
    real(d), dimension(nsta)      :: v
    real(d), dimension(nsta,nsta) :: w
    real(d), dimension(3*nsta)    :: work

!----------------------------------------------------------------------
! Construct the LVC potential
!----------------------------------------------------------------------
    w=potfunc(q)

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
    use n2m_global

    implicit none

    integer                       :: i
    real(d), dimension(nmodes)    :: q
    real(d), dimension(nsta)      :: wii
    real(d), dimension(nsta,nsta) :: w

!----------------------------------------------------------------------
! Construct the LVC potential
!----------------------------------------------------------------------
    w=potfunc(q)

!----------------------------------------------------------------------
! Return the on-diagonal elements
!----------------------------------------------------------------------
    do i=1,nsta
       wii(i)=w(i,i)
    enddo
    
    return

  end function diabpot
  
!######################################################################

  function potfunc(q) result(w)

    use constants
    use n2m_global

    implicit none

    integer                       :: m,m1,m2,s,s1,s2
    real(d), dimension(nmodes)    :: q
    real(d), dimension(nsta,nsta) :: w

!----------------------------------------------------------------------
! Initialisation of the diabatic potential
!----------------------------------------------------------------------
    w=0.0d0

!----------------------------------------------------------------------
! Zeroth-order contributions
!----------------------------------------------------------------------
    ! Energies
    do s1=1,nsta
       do s2=1,nsta
          w(s1,s2)=w(s1,s2)+const(s1,s2)
       enddo
    enddo

    ! Harmonic potentials
    do s=1,nsta
       do m=1,nmodes
          w(s,s)=w(s,s)+0.5d0*freq(m)*q(m)**2
       enddo
    enddo

!----------------------------------------------------------------------
! First-order terms
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

!----------------------------------------------------------------------
! Second-order terms: note that the nadvibs coefficient already
! account for the 1/n! prefactor and also assume a restricted summation
! over unique n-tuples of mode indices
!----------------------------------------------------------------------
    ! gamma
    do s=1,nsta
       do m1=1,nmodes
          do m2=m1,nmodes
             w(s,s)=w(s,s)+gamma(m1,m2,s)*q(m1)*q(m2)
          enddo
       enddo
    enddo

    ! mu
    do s1=1,nsta-1
       do s2=s1+1,nsta
          do m1=1,nmodes
             do m2=m1,nmodes
                w(s1,s2)=w(s1,s2)+mu(m1,m2,s1,s2)*q(m1)*q(m2)
                w(s2,s1)=w(s2,s1)+mu(m1,m2,s1,s2)*q(m1)*q(m2)
             enddo
          enddo
       enddo
    enddo
       
    return

  end function potfunc

!######################################################################

  subroutine wrgnuplot

    use constants
    use iomod
    use n2m_global
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
    write(unit,'(a)') 'monitorSize=system("xrandr | awk &
         ''/\*/{sub(/x/,\",\");print $1;exit}''")'
    write(unit,'(a,/)') 'set terminal x11 size @monitorSize'

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

    use constants
    use iomod
    use channels
    use n2m_global
    use pltmod
    
    implicit none

!----------------------------------------------------------------------
! Close the log file
!----------------------------------------------------------------------
    close(ilog)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(freq)
    deallocate(const)
    deallocate(kappa)
    deallocate(lambda)
    deallocate(gamma)
    deallocate(mu)
    deallocate(surf)

    return

  end subroutine finalise

!######################################################################
    
end program pltnadvibs
