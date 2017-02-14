module ioqc

  implicit none

contains

!######################################################################

  subroutine getdim

    use constants
    use channels
    use iomod
    use global

    implicit none

    integer            :: unit,n
    character(len=150) :: filename,string

!----------------------------------------------------------------------
! Determine the no. states from ciudgsm.drt1.sp
!----------------------------------------------------------------------
    ! Open file
    call freeunit(unit)
    filename=trim(coldir)//'/LISTINGS/ciudgsm.drt1.sp'
    open(unit,file=filename,form='formatted',status='old')

    ! Read to the MRCI energy section
    n=0
5   read(unit,'(a)',end=111) string
    if (index(string,'final mr-sdci').ne.0) then
       n=n+1
       if (n.ne.2) goto 5
    else
       goto 5
    endif

    ! Determine the no. states
    nsta=0
10  read(unit,'(a)') string
    if (index(string,'mr-sdci').ne.0) then
       nsta=nsta+1
       goto 10
    endif

    ! Close file
    close(unit)

!----------------------------------------------------------------------
! Determine the no. atoms from the geom file
!----------------------------------------------------------------------
    ! Open file
    call freeunit(unit)
    filename=trim(coldir)//'/geom'
    open(unit,file=filename,form='formatted',status='old')

    ! Determine the no. atoms
    natm=0
15  read(unit,'(a)',end=20) string
    if (index(string,'X').eq.0) natm=natm+1
    goto 15
20  continue
    
    ! No. Cartesian coordinates
    ncoo=3*natm

    ! No. normal modes
    ! N.B. for now we assume that the molecule is not linear...
    nmodes=ncoo-6

    ! Close file
    close(unit)

!----------------------------------------------------------------------
! Output the system dimensions to the log file
!----------------------------------------------------------------------
    write(ilog,'(/,2x,a,2x,i2)') "No. states:",nsta
    write(ilog,'(/,2x,a,2x,i2)') "No. atoms:",natm

    return

111 continue
    errmsg='It looks like the MRCI calculation did not complete...'
    call error_control

  end subroutine getdim

!######################################################################

  subroutine rdgrad

    use constants
    use global
    use iomod

    implicit none

    integer            :: i,j,k,unit
    character(len=150) :: filename
    character(len=2)   :: asta

!----------------------------------------------------------------------
! Read the energy gradients in terms of the Cartesian coordinates
!----------------------------------------------------------------------
    call freeunit(unit)

    ! Loop over states
    do i=1,nsta

       ! Open the gradient file
       write(asta,'(i2)') i
       filename=trim(coldir)//'/GRADIENTS/cartgrd.drt1.state'&
            //trim(adjustl(asta))//'.sp'
       open(unit,file=filename,form='formatted',status='old')

       ! Read the gradient file
       do j=1,natm
          read(unit,*) (grad(k,i), k=j*3-2,j*3)
       enddo

       ! Close the gradient file
       close(unit)

    enddo

    return

  end subroutine rdgrad

!######################################################################

  subroutine rdnact

    use constants
    use global
    use iomod

    implicit none

    integer            :: s1,s2,i,j,unit
    character(len=150) :: filename
    character(len=2)   :: asta1,asta2

!----------------------------------------------------------------------
! Read the NACTS in terms of the Cartesian coordinates
!----------------------------------------------------------------------
    call freeunit(unit)

    ! Loop unique pairs of states
    do s1=1,nsta-1
       do s2=s1+1,nsta

          ! Open the NACT file
          write(asta1,'(i2)') s1
          write(asta2,'(i2)') s2
          filename=trim(coldir)//'/GRADIENTS/cartgrd.nad.drt1.state'&
            //trim(adjustl(asta1))//'.drt1.state'&
            //trim(adjustl(asta2))//'.sp'
          open(unit,file=filename,form='formatted',status='old')

          ! Read the NACT file
          do i=1,natm
             read(unit,*) (nact(j,s1,s2), j=i*3-2,i*3)
          enddo

          ! Fill in the lower triangle
          nact(:,s2,s1)=-nact(:,s1,s2)

          ! Close the gradient file
          close(unit)

       enddo
    enddo

    return

  end subroutine rdnact

!######################################################################

  subroutine rdfreqfile

    use constants
    use global
    use iomod

    implicit none

!----------------------------------------------------------------------
! Determine the quantum chemistry program used for the frequency
! calculation
!----------------------------------------------------------------------
    call freqtype

!----------------------------------------------------------------------
! Read the reference Cartesian coordinates
!----------------------------------------------------------------------
    call getxcoo

!----------------------------------------------------------------------
! Read the normal modes, frequencies and symmetry labels
!----------------------------------------------------------------------
    call getmodes

    return
    
  end subroutine rdfreqfile

!######################################################################
! freqtyp: determines the quantum chemistry program used for the
!          frequency calculation, and sets ityp accordingly:
!          
!          ityp = 1 <-> G98
!
!######################################################################
  subroutine freqtype

    use constants
    use global
    use iomod

    implicit none
    
    integer            :: unit
    character(len=120) :: string

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ityp=-1

!----------------------------------------------------------------------
! Open the frequency file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=freqfile,form='formatted',status='old')

!----------------------------------------------------------------------
! Determine the quantum chemistry program used for the frequency
! calculation
!----------------------------------------------------------------------    
    read(unit,'(a)') string

    ! G98
    if (index(string,'Entering Gaussian System').ne.0) then
       ityp=1
    endif

!----------------------------------------------------------------------
! Check that the quantum chemistry program used is supported
!----------------------------------------------------------------------
    if (ityp.eq.-1) then
       errmsg='The quantum chemistry program used for the &
            frequency calculation is not supported.'
       call error_control
    endif

!----------------------------------------------------------------------
! Close the frequency file
!----------------------------------------------------------------------
    close(unit)

    return

  end subroutine freqtype

!######################################################################

  subroutine getxcoo

    use constants
    use global
    use iomod

    implicit none

    if (ityp.eq.1) then
       ! G98
       call getxcoo_g98
    endif

    return

  end subroutine getxcoo

!######################################################################
  
  subroutine getxcoo_g98

    use constants
    use global
    use iomod

    implicit none
  
    integer            :: unit,i,j
    character(len=120) :: string

!----------------------------------------------------------------------
! Open the frequency file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=freqfile,form='formatted',status='old')

!----------------------------------------------------------------------
! Read the Cartesian coordinates
!----------------------------------------------------------------------
5   read(unit,'(a)',end=999) string
    if (index(string,'Z-Matrix orientation').eq.0) goto 5

    do i=1,4
       read(unit,*)
    enddo

    do i=1,natm
       read(unit,'(14x,i2,18x,3F12.6)') atnum(i),(xcoo0(j), j=i*3-2,i*3)
       atlbl(i)=num2lbl(atnum(i))
       mass(i*3-2:i*3)=num2mass(atnum(i))       
    enddo

    ! Convert to Bohr
    xcoo0=xcoo0*ang2bohr

!----------------------------------------------------------------------
! Close the frequency file
!----------------------------------------------------------------------
    close(unit)

    return

999 continue
    errmsg='The Cartesian coordinates could not be found in: '&
         //trim(freqfile)
    call error_control

  end subroutine getxcoo_g98

!######################################################################

  function num2lbl(num) result(lbl)

    use constants
    use iomod

    implicit none

    integer          :: num
    character(len=2) :: lbl

    if (num.eq.1) then
       lbl='H'
    else if (num.eq.6) then
       lbl='C'
    else
       write(errmsg,'(a,x,i2,x,a)') 'The atomic number',num,&
            'is not supported. See function num2lbl.'
       call error_control
    endif

    return

  end function num2lbl

!######################################################################

  function num2mass(num) result(mass)

    use constants
    use iomod

    implicit none

    integer :: num
    real(d) :: mass

    if (num.eq.1) then
       mass=1.00794d0
    else if (num.eq.6) then
       mass=12.0107d0
    else
       write(errmsg,'(a,x,i2,x,a)') 'The atomic number',num,&
            'is not supported. See function num2mass.'
       call error_control
    endif

    return

  end function num2mass

!######################################################################

  subroutine getmodes

    use constants
    use global
    use iomod

    implicit none

!----------------------------------------------------------------------
! Read the normal modes from file
!----------------------------------------------------------------------
    if (ityp.eq.1) then
       ! G98
       call getmodes_g98
    endif

!----------------------------------------------------------------------
! Orthogonalise the normal mode vectors
!----------------------------------------------------------------------
    call nmortho

    return

  end subroutine getmodes

!######################################################################

  subroutine getmodes_g98

    use constants
    use global
    use iomod

    implicit none

    integer                          :: unit
    integer                          :: i,j,n,nm,nrd,ncurr
    real(d), dimension(ncoo,ncoo)    :: matrix
    real(d), dimension(10)           :: xdum
    real(d)                          :: s2
    logical                          :: found
    character(len=120)               :: string
    character(len=16), dimension(10) :: buffer

!----------------------------------------------------------------------
! Open the frequency file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=freqfile,form='formatted',status='old')

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    matrix=0.0d0

!----------------------------------------------------------------------
! Get the frequencies
!----------------------------------------------------------------------
100 continue

    read(unit,'(a)',end=999) string

    if (string(2:30) .ne. 'Harmonic frequencies (cm**-1)') go to 100

    nrd = ((nmodes-1)/5) + 1

110 continue
    read(unit,'(a)',end=900) string
    if (string(8:22) .ne. 'Frequencies ---') go to 110

    found=.true.
    backspace(unit)
    backspace(unit)
    backspace(unit)

    nm=0
    do n=1,nrd
       ncurr=min(5,nmodes-5*(n-1))
       read(unit,*)
       read(unit,'(26x,5a10)') (buffer(j),j=1,ncurr)
       read(unit,'(23x,5f10.4)',err=900) (xdum(j),j=1,ncurr)       
       do j=1,ncurr
          nm=nm+1
          nmlab(nm)=buffer(j)
          freq(nm)=xdum(j)
       enddo
120    continue
       read(unit,'(a)',end=999) string
       if (string(2:20) .ne. 'Coord Atom Element:') go to 120
       do j=1,ncoo
          read(unit,'(a)',err=900) string
          read(string,20) (xdum(i),i=1,ncurr)
          nm=nm-ncurr
          do i=1,ncurr
             nm=nm+1
             matrix(j,nm)=xdum(i)
          enddo
       enddo
    enddo
20  format(23x,5f10.5)

!----------------------------------------------------------------------
! Close the frequency file
!----------------------------------------------------------------------
    close(unit)

!-----------------------------------------------------------------------
! Scale to mass-weighted x to obtain true normal mode vectors
!-----------------------------------------------------------------------
    do nm=1,nmodes
       do j=1,ncoo
          matrix(j,nm)=matrix(j,nm)*sqrt(mass(j))
       enddo
    enddo

!-----------------------------------------------------------------------
! Normalise the normal mode vectors
!-----------------------------------------------------------------------
    do i=1,nmodes
       matrix(:,i)=matrix(:,i) &
            /sqrt(dot_product(matrix(:,i),matrix(:,i)))
    enddo

!-----------------------------------------------------------------------
! Convert frequencies to ev
!-----------------------------------------------------------------------
    do i=1,nmodes
       freq(i)=freq(i)*invcm2ev
    enddo

!-----------------------------------------------------------------------
! Save the normal modes in the nmcoo array
!-----------------------------------------------------------------------
    nmcoo=matrix(1:ncoo,1:nmodes)

    return

900 continue
    errmsg='Incorrect normal mode format. Use IOP(7/8=11)'
    call error_control

999 continue
    errmsg='Problem reading the normal modes in: '//trim(freqfile)
    call error_control

  end subroutine getmodes_g98

!######################################################################

  subroutine nmortho

    use constants
    use global
    use channels

    implicit none

    integer :: i,i1,j
    real(d) :: x
      
    do i=1,nmodes
       
       do i1=1,i-1
          x=0.0d0
          do j=1,ncoo
             x=x+nmcoo(j,i)*nmcoo(j,i1)
          enddo
          do j=1,ncoo
             nmcoo(j,i)=nmcoo(j,i)-x*nmcoo(j,i1)
          enddo
       enddo
    
       x=0.0d0
       do j=1,ncoo
          x=x+nmcoo(j,i)*nmcoo(j,i)
       enddo
       
       x=1.0d0/sqrt(x)
       do j=1,ncoo
          nmcoo(j,i)=x*nmcoo(j,i)
       enddo
    
    enddo

    return

  end subroutine nmortho

!#######################################################################
! nmcoo transforms from nmodes to coo  x = nmcoo*Q
! coonm transforms from coo to nmodes  Q = coonm*x
!
! freq in ev, mass in amu
!#######################################################################

  subroutine nm2xmat

    use constants
    use global

    implicit none

    integer :: i,j
      
!-----------------------------------------------------------------------
! Transposition of nmcoo
!-----------------------------------------------------------------------
    do j=1,nmodes         
       coonm(j,:)=nmcoo(:,j)
    enddo
      
!-----------------------------------------------------------------------
! Multiply matrix by mass and frequency factors
! "frequencies" in eV and masses in amu
!-----------------------------------------------------------------------
    do j=1,ncoo
       do i=1,nmodes
          nmcoo(j,i)=nmcoo(j,i)/(15.4644d0*sqrt(freq(i)*mass(j)))
          coonm(i,j)=coonm(i,j)*(15.4644d0*sqrt(freq(i)*mass(j)))
       enddo
    enddo
    
    return
  
  end subroutine nm2xmat

!######################################################################

end module ioqc
