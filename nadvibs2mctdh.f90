!#######################################################################
! nadvibs2mctdh: a program to convert a NADVIBS input file to an MCTDH
!                operator file
!#######################################################################

program nadvibs2mctdh

  implicit none

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
  call initialise

!----------------------------------------------------------------------
! Read the input
!----------------------------------------------------------------------
  call rdinput

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
! Write the MCTDH operator file
!----------------------------------------------------------------------
  call wroper
  
!----------------------------------------------------------------------
! Finalisation
!----------------------------------------------------------------------
  call finalise
  
contains

!######################################################################

  subroutine initialise

    use channels
    use iomod
    
    implicit none

!----------------------------------------------------------------------
! Open the log file
!----------------------------------------------------------------------
    call freeunit(ilog)
    open(ilog,file='nadvibs2mctdh.log',form='formatted',&
         status='unknown')
    
    return
    
  end subroutine initialise

!######################################################################

  subroutine rdinput

    use constants
    use iomod
    use n2m_global
    
    implicit none

    integer            :: n
    character(len=120) :: string1,string2
    
!----------------------------------------------------------------------
! Set default values
!----------------------------------------------------------------------
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
    
    return
    
  end subroutine rdinput
  
!######################################################################

  subroutine alloc

    use constants
    use iomod
    use n2m_global
    
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
    
    return
    
  end subroutine alloc

!######################################################################

  subroutine finalise

    use channels
    use n2m_global
    
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
    
    return
    
  end subroutine finalise

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
                do m=1,nmodes
                   kappa(m,i)=tmp1(m)
                enddo
             else
                do m=1,nmodes
                   lambda(m,i,j)=tmp1(m)
                   lambda(m,j,i)=tmp1(m)
                enddo
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

  subroutine wroper

    use constants
    use iomod
    use n2m_global
    
    implicit none

    integer            :: unit,fel,m,m1,m2,s,s1,s2,i,k,nl,ncurr
    real(d), parameter :: thrsh=5e-4
    character(len=2)   :: am,am1,am2,as,as1,as2,afel
    character(len=90)  :: string
    character(len=3)   :: amode
    character(len=8)   :: atmp
!----------------------------------------------------------------------
! Open the MCTDH operator file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file='mctdh.op',form='formatted',status='unknown')

!----------------------------------------------------------------------
! Index of the electronic DOF
!----------------------------------------------------------------------
    fel=nmodes+1
    write(afel,'(i2)') fel
    
!----------------------------------------------------------------------
! Write the op_define section
!----------------------------------------------------------------------
    write(unit,'(a)') 'op_define-section'
    write(unit,'(a)') 'title'
    write(unit,'(a)') 'MCTDH operator file created by nadvibs2mctdh'
    write(unit,'(a)') 'Note that NADVIBS assumes a restricted &
         summation over mode indices'
    write(unit,'(a)') 'and that the 1/n! prefactors are folded into &
         the parameter values'
    write(unit,'(a)') 'end-title'
    write(unit,'(a)') 'end-op_define-section'

!----------------------------------------------------------------------
! Write the parameter section
!----------------------------------------------------------------------
    ! Starting line
    write(unit,'(/,a)') 'parameter-section'

    ! Frequencies
    write(unit,'(/,a)') '# Frequencies'
    do m=1,nmodes
       write(am,'(i2)') m
       write(unit,'(a,F9.6,a)') &
            'omega_'//adjustl(am)//' =',freq(m),' , ev'
    enddo

    ! Energies
    write(unit,'(/,a)') '# Energies'
    do s=1,nsta
       write(as,'(i2)') s
       write(unit,'(a,F9.6,a)') &
            'E'//adjustl(as)//' = ',const(s,s),' , ev'
    enddo
    
    ! 1st-order intrastate coupling constants (kappa)
    write(unit,'(/,a)') '# 1st-order intrastate coupling &
         constants (kappa)'
    do s=1,nsta
       write(as,'(i2)') s
       do m=1,nmodes
          if (abs(kappa(m,s)).lt.thrsh) cycle
          write(am,'(i2)') m
          write(unit,'(a,F9.6,a)') 'kappa'//trim(adjustl(as))&
               //'_'//adjustl(am)//' = ',kappa(m,s),' , ev'
       enddo
    enddo

    ! 1st-order interstate coupling constants (lambda)
    write(unit,'(/,a)') '# 1st-order interstate coupling &
         constants (lambda)'
    do s1=1,nsta-1
       write(as1,'(i2)') s1
       do s2=s1+1,nsta
          write(as2,'(i2)') s2
          do m=1,nmodes
             if (abs(lambda(m,s1,s2)).lt.thrsh) cycle
             write(am,'(i2)') m
             write(unit,'(a,F9.6,a)') 'lambda'//trim(adjustl(as1))&
               //'_'//trim(adjustl(as2))//'_'//adjustl(am)//&
               ' = ',lambda(m,s1,s2),' , ev'
          enddo
       enddo
    enddo

    ! 2nd-order intrastate coupling constants (gamma)
    if (order.ge.2) then
       write(unit,'(/,a)') '# 2nd-order intrastate coupling &
         constants (gamma)'
       do s=1,nsta
          write(as,'(i2)') s
          do m1=1,nmodes
             write(am1,'(i2)') m1
             do m2=m1,nmodes
                if (abs(gamma(m1,m2,s)).lt.thrsh) cycle
                write(am2,'(i2)') m2
                write(unit,'(a,F9.6,a)') 'gamma'//trim(adjustl(as))&
               //'_'//trim(adjustl(am1))//'_'//adjustl(am2)//&
               ' = ',gamma(m1,m2,s),' , ev'
             enddo
          enddo
       enddo
    endif

    ! 2nd-order interstate coupling constants (mu)
    if (order.ge.2) then
       write(unit,'(/,a)') '# 2nd-order interstate coupling &
         constants (mu)'
       do s1=1,nsta-1
          write(as1,'(i2)') s1
          do s2=s1+1,nsta
             write(as2,'(i2)') s2
             do m1=1,nmodes
                write(am1,'(i2)') m1
                do m2=m1,nmodes
                   if (abs(mu(m1,m2,s1,s2)).lt.thrsh) cycle
                   write(am2,'(i2)') m2
                   write(unit,'(a,F9.6,a)') 'mu'//trim(adjustl(as1))&
                        //'_'//trim(adjustl(as2))&
                        //'_'//trim(adjustl(am1))//'_'//adjustl(am2)//&
                        ' = ',mu(m1,m2,s1,s2),' , ev'
                enddo
             enddo
          enddo
       enddo
    endif
    
    ! Finishing line
    write(unit,'(/,a)') 'end-parameter-section'

!----------------------------------------------------------------------
! Write the Hamiltonian section
!----------------------------------------------------------------------
    ! Starting line
    write(unit,'(/,a)') 'hamiltonian-section'

    ! Modes section
    write(unit,'(/,38a)') ('-',i=1,38)
    m=0
    nl=ceiling((real(nmodes+1))/10.0d0)
    do i=1,nl
       ncurr=min(10,nmodes+1-10*(i-1))
       string='modes|'
       do k=1,ncurr
          m=m+1
          if (m.lt.nmodes+1) then
             write(amode,'(i3)') m
             write(atmp,'(a)') ' v'//adjustl(amode)//' '
             if (m.lt.10) then
                string=trim(string)//trim(atmp)//'  |'
             else
                string=trim(string)//trim(atmp)//' |'
             endif
          else
             if (m.eq.nmodes+1) then
                string=trim(string)//' el'
             else
                string=trim(string)//'  | Time'
             endif
          endif
       enddo
       write(unit,'(a)') trim(string)
    enddo
    write(unit,'(38a)') ('-',i=1,38)

    ! Kinetic energy operator
    write(unit,'(/,a)') '# Kinetic energy'
    do m=1,nmodes
       write(am,'(i2)') m
       write(unit,'(a)') &
            'omega_'//adjustl(am)//'  |'//adjustl(am)//'  KE'
    enddo

    ! Zeroth-order potential: VEEs
    write(unit,'(/,a)') '# Zeroth-order potential: VEEs'
    do s=1,nsta
       write(as,'(i2)') s
       write(unit,'(a)') &
            'E'//adjustl(as)//'  |'//adjustl(afel)&
            //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
    enddo

    ! Zeroth-order potential: Harmonic oscillators
    write(unit,'(/,a)') '# Zeroth-order potential: &
         Harmonic oscillators'
    do m=1,nmodes
       write(am,'(i2)') m
       write(unit,'(a)') &
            '0.5*omega_'//adjustl(am)//'  |'//adjustl(am)//'  q^2'
    enddo

    ! 1st-order intrastate coupling terms (kappa)
    write(unit,'(/,a)') '# 1st-order intrastate coupling &
         constants (kappa)'
    do m=1,nmodes
       write(am,'(i2)') m
       do s=1,nsta
          if (abs(kappa(m,s)).lt.thrsh) cycle
          write(as,'(i2)') s
          write(unit,'(a)') 'kappa'//trim(adjustl(as))&
               //'_'//adjustl(am)&
               //'  |'//adjustl(am)//'  q'//'  |'//adjustl(afel)&
               //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
       enddo
    enddo

    ! 1st-order interstate coupling terms (lambda)
    write(unit,'(/,a)') '# 1st-order interstate coupling &
         constants (lambda)'
    do s1=1,nsta-1
       write(as1,'(i2)') s1
       do s2=s1+1,nsta
          write(as2,'(i2)') s2
          do m=1,nmodes
             write(am,'(i2)') m
             if (abs(lambda(m,s1,s2)).lt.thrsh) cycle
             write(unit,'(a)') 'lambda'//trim(adjustl(as1)) &
                  //'_'//trim(adjustl(as2))//'_'//adjustl(am) &
                  //'  |'//adjustl(am)//'  q'//'  |'//adjustl(afel)&
                  //'  S'//trim(adjustl(as1))//'&'//trim(adjustl(as2))
          enddo
       enddo
    enddo

    ! Quadratic 2nd-order intrastate coupling terms (gamma)
    write(unit,'(/,a)') '# Quadratic 2nd-order intrastate &
         coupling constants (gamma)'
    do s=1,nsta
       write(as,'(i2)') s
       do m=1,nmodes
          write(am,'(i2)') m
          if (abs(gamma(m,m,s)).lt.thrsh) cycle
          write(unit,'(a)') 'gamma'//trim(adjustl(as))&
               //'_'//trim(adjustl(am))//'_'//adjustl(am)&
               //'  |'//adjustl(am)//'  q^2'//'  |'//adjustl(afel)&
               //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
       enddo
    enddo

    ! Bi-linear 2nd-order intrastate coupling terms (gamma)
    write(unit,'(/,a)') '# Bi-linear 2nd-order intrastate &
         coupling constants (gamma)'
    do s=1,nsta
       write(as,'(i2)') s
       do m1=1,nmodes-1
          write(am1,'(i2)') m1
          do m2=m1+1,nmodes
             write(am2,'(i2)') m2             
             if (abs(gamma(m1,m2,s)).lt.thrsh) cycle
             write(unit,'(a)') 'gamma'//trim(adjustl(as))&
                  //'_'//trim(adjustl(am1))//'_'//adjustl(am2)&
                  //'  |'//adjustl(am1)//'  q'&
                  //'  |'//adjustl(am2)//'  q'//'  |'//adjustl(afel)&
                  //'  S'//trim(adjustl(as))//'&'//trim(adjustl(as))
          enddo
       enddo
    enddo

    ! Quadratic 2nd-order interstate coupling terms (mu)
    write(unit,'(/,a)') '# Quadratic 2nd-order interstate &
         coupling constants (mu)'
    do s1=1,nsta-1
       write(as1,'(i2)') s1
       do s2=s1+1,nsta
          write(as2,'(i2)') s2
          do m=1,nmodes
             write(am,'(i2)') m
             if (abs(mu(m,m,s1,s2)).lt.thrsh) cycle
             write(unit,'(a)') 'mu'//trim(adjustl(as1))&
                  //'_'//trim(adjustl(as2))&
                  //'_'//trim(adjustl(am))//'_'//adjustl(am)&
                  //'  |'//adjustl(am)//'  q^2'//'  |'//adjustl(afel)&
                  //'  S'//trim(adjustl(as1))//'&'//trim(adjustl(as2))
          enddo
       enddo
    enddo

    ! Bi-linear 2nd-order interstate coupling terms (mu)
    write(unit,'(/,a)') '# Bi-linear 2nd-order interstate &
         coupling constants (mu)'
    do s1=1,nsta-1
       write(as1,'(i2)') s1
       do s2=s1+1,nsta
          write(as2,'(i2)') s2
          do m1=1,nmodes-1
             write(am1,'(i2)') m1
             do m2=m1+1,nmodes
                write(am2,'(i2)') m2
                if (abs(mu(m1,m2,s1,s2)).lt.thrsh) cycle
                write(unit,'(a)') 'mu'//trim(adjustl(as1))&
                     //'_'//trim(adjustl(as2))&
                     //'_'//trim(adjustl(am1))//'_'//adjustl(am2)&
                     //'  |'//adjustl(am1)//'  q'&
                     //'  |'//adjustl(am2)//'  q'&
                     //'  |'//adjustl(afel)&
                     //'  S'//trim(adjustl(as1))//'&'//trim(adjustl(as2))
             enddo
          enddo
       enddo
    enddo

    write(unit,'(/,a)') 'end-hamiltonian-section'
    
!----------------------------------------------------------------------
! End line
!----------------------------------------------------------------------
    write(unit,'(/,a)') 'end-operator'

!----------------------------------------------------------------------
! Close the MCTDH operator file
!----------------------------------------------------------------------
    close(unit)
    return
    
  end subroutine wroper
    
!######################################################################
  
end program nadvibs2mctdh
