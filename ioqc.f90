module ioqc

  implicit none

  save

  integer :: idum
  logical :: ldum

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

  subroutine rdener

    use constants
    use iomod
    use global

    implicit none

    integer            :: unit,n
    character(len=150) :: filename,string

!----------------------------------------------------------------------
! Read the Davidson corrected MRCI energies from ciudgsm.drt1.sp
!----------------------------------------------------------------------
    ! Open file
    call freeunit(unit)
    filename=trim(coldir)//'/LISTINGS/ciudgsm.drt1.sp'
    open(unit,file=filename,form='formatted',status='old')

    ! Read the Davdson-corrected MRCI energies
    n=0
5   read(unit,'(a)',end=999) string
    if (index(string,'eci+dv3').ne.0) then
       n=n+1
       read(string,'(12x,F20.12)') ener(n)
       if (n.lt.nsta) goto 5
    else
       goto 5
    endif

    ! Close file
    close(unit)

    ! Excitation energies in eV
    do n=2,nsta
       ener(n)=(ener(n)-ener(1))*eh2ev
    enddo
    ener(1)=0.0d0

    return

999 continue
    errmsg='Not all MRCI energies could be found... Quitting!'
    call error_control

  end subroutine rdener
    
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
! Read the NACTs (multiplied by the energy difference!) in terms of 
! the Cartesian coordinates
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
          nact(:,s2,s1)=nact(:,s1,s2)

          ! Close the NACT file
          close(unit)

       enddo
    enddo

    return

  end subroutine rdnact

!######################################################################

  subroutine rddip

    use constants
    use global
    use iomod

    implicit none

    integer            :: s,s1,s2,i,j,unit
    character(len=150) :: filename,string
    character(len=2)   :: asta1,asta2

!----------------------------------------------------------------------
! Read the on-diagonal elements of the dipole matrix
!----------------------------------------------------------------------
    call freeunit(unit)

    ! Loop over states
    do s=1,nsta
       
       ! Open the propls file
       write(asta1,'(i2)') s
       filename=trim(coldir)//'/LISTINGS/propls.ci.drt1.state'&
            //trim(adjustl(asta1))//'.sp'
       open(unit,file=filename,form='formatted',status='old')
       
       ! Read the dipole moment
5      read(unit,'(a)',end=999) string
       if (index(string,'Dipole moments:').eq.0) goto 5
       do i=1,3
          read(unit,*)
       enddo
       read(unit,'(11x,3(2x,F14.8))') (dipole(i,s,s),i=1,3)

       ! Close the propls file
       close(unit)

    enddo

!----------------------------------------------------------------------
! Read the off-diagonal elements of the dipole matrix
!----------------------------------------------------------------------
    ! Loop unique pairs of states
    do s1=1,nsta-1
       do s2=s1+1,nsta

          ! Open the trncils file
          write(asta1,'(i2)') s1
          write(asta2,'(i2)') s2
          filename=trim(coldir)//'/LISTINGS/trncils.drt1.state'&
            //trim(adjustl(asta1))//'.drt1.state'&
            //trim(adjustl(asta2))
          open(unit,file=filename,form='formatted',status='old')

          ! Read the transition diole moment
10        read(unit,'(a)',end=999) string
          if (index(string,'Transition moment components:').eq.0) goto 10
          do i=1,3
             read(unit,*)
          enddo
          read(unit,'(13x,3(F13.6))') (dipole(i,s1,s2),i=1,3)
          dipole(:,s2,s1)=dipole(:,s1,s2)

          ! Close the trncils file
          close(unit)

       enddo
    enddo

    return

999 continue
    errmsg='No dipole moment matrix element could be found in: '&
         //trim(filename)
    call error_control

  end subroutine rddip

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
!                 2 <-> CFOUR
!######################################################################
  subroutine freqtype

    use constants
    use global
    use iomod

    implicit none
    
    integer            :: unit
    character(len=120) :: string
    logical            :: dir

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ityp=-1

!----------------------------------------------------------------------
! Determine the quantum chemistry program used for the frequency
! calculation
!----------------------------------------------------------------------
    if (isg98(freqfile)) then
       ! G98
       ityp=1
    else if (iscfour(freqfile)) then
       ! CFOUR
       ityp=2
    endif

!----------------------------------------------------------------------
! Check that the quantum chemistry program used is supported
!----------------------------------------------------------------------
    if (ityp.eq.-1) then
       errmsg='The quantum chemistry program used for the &
            frequency calculation is not supported.'
       call error_control
    endif

    return

  end subroutine freqtype

!######################################################################

  function isg98(filename) result(found)

    use constants
    use iomod

    implicit none

    integer            :: unit
    character(len=*)   :: filename
    character(len=120) :: string
    logical            :: found,dir

!----------------------------------------------------------------------
! First determine whether freqfile is actually a file.
! This is necessary as for certain programs, the name of a directory
! will actually be passed instead
!----------------------------------------------------------------------
    inquire(file=trim(filename)//'/.',exist=dir)

    if (dir) then
       found=.false.
       return
    endif

!----------------------------------------------------------------------
! Open the frequency file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='formatted',status='old')

!----------------------------------------------------------------------
! Check whether the calculation was performed using G98
!----------------------------------------------------------------------
    read(unit,'(a)') string
    if (index(string,'Entering Gaussian System').ne.0) then
       found=.true.
    else
       found=.false.
    endif

!----------------------------------------------------------------------
! Close the frequency file
!----------------------------------------------------------------------
    close(unit)

    return

  end function isg98

!######################################################################

  function iscfour(filename) result(found)

    use constants
    use iomod

    implicit none

    integer            :: unit
    character(len=*)   :: filename
    character(len=120) :: string
    logical            :: found,dir

!----------------------------------------------------------------------
! First determine whether freqfile is actually a file.
! This is necessary as for certain programs, the name of a directory
! will actually be passed instead
!----------------------------------------------------------------------
    inquire(file=trim(filename)//'/.',exist=dir)

    if (.not.dir) then
       found=.false.
       return
    endif

!----------------------------------------------------------------------
! Check to see if the cfour FCM file exists
!----------------------------------------------------------------------
    inquire(file=trim(filename)//'/FCM',exist=found)

    return

  end function iscfour

!######################################################################

  subroutine getxcoo

    use constants
    use global
    use iomod

    implicit none

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ldum=.false.
    idum=0

!----------------------------------------------------------------------
! Read in the reference Cartesian coordinates
!----------------------------------------------------------------------
    if (ityp.eq.1) then
       ! G98
       call getxcoo_g98
    else if (ityp.eq.2) then
       ! CFOUR
       call getxcoo_cfour
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

  subroutine getxcoo_cfour

    use constants
    use iomod
    use global

    implicit none

    integer            :: unit,i,n
    logical            :: found
    character(len=120) :: string
    character(len=2)   :: atmp

!----------------------------------------------------------------------
! First check to make sure that the cfour.log file exists
!----------------------------------------------------------------------
    inquire(file=trim(freqfile)//'/cfour.log',exist=found)

    if (.not.found) then
       errmsg='The CFOUR log file is assumed to be named cfour.log. &
            This file could not be found in: '//trim(freqfile)
       call error_control
    endif

!----------------------------------------------------------------------
! Open the CFOUR log file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=trim(freqfile)//'/cfour.log',form='formatted',&
         status='old')

!----------------------------------------------------------------------
! Read the reference Cartesian coordinates from the CFOUR log file
!----------------------------------------------------------------------
5   read(unit,'(a)',end=999) string
    if (index(string,'Coordinates (in bohr)').eq.0) goto 5

    do i=1,2
       read(unit,*)
    enddo

    n=0
10  read(unit,'(a)') string
    if (index(string,'---').eq.0) then
       if (string(6:6).ne.'X') then
          n=n+1
          read(string,'(5x,a2,8x,i2,3x,3(F15.8))') &
               atlbl(n),atnum(n),(xcoo0(i), i=n*3-2,n*3)
          mass(n*3-2:n*3)=num2mass(atnum(n))
       else
          ldum=.true.
          idum=n+1
       endif
       goto 10
    endif

!----------------------------------------------------------------------
! Close the CFOUR log file
!----------------------------------------------------------------------
    close(unit)

    return

999 continue
    errmsg='The Cartesian coordinates could not be found in: '&
         //trim(freqfile)//'/cfour.log'
    call error_control

  end subroutine getxcoo_cfour

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
    else if (num.eq.7) then
       mass=14.0067d0
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
    else if (ityp.eq.2) then
       ! CFOUR
       call getmodes_cfour
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
! Convert frequencies to eV
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

  subroutine getmodes_cfour

    use constants
    use iomod
    use global

    implicit none

    integer                           :: unit,i,j,k,lim1,lim2
    real(d), dimension(ncoo,ncoo)     :: hess
    real(d), dimension(ncoo+3,ncoo+3) :: dumhess
    real(d), dimension(nmodes)        :: cffreq
    character(len=120)                :: string

!**********************************************************************
! Note that cfour outputs the normal mode vectors to only a few
! decimal places. Therefore, we construct the normal modes ourself
! from the Hessian, which is written to a decent level of precision
! in the FCMFINAL file.
!**********************************************************************


!----------------------------------------------------------------------
! Read the Hessian from the FCMFINAL file
!----------------------------------------------------------------------
    ! Open the FCMFINAL file
    call freeunit(unit)
    open(unit,file=trim(freqfile)//'/FCMFINAL',form='formatted',&
         status='old')

    ! Read the FCMFINAL file
    read(unit,*)

    if (ldum) then
       ! Dummy atom present
       ! (1) Read the Hessian including the dummy atom
       do i=1,ncoo+3
          do j=1,natm+1
             read(unit,'(3F20.10)') (dumhess(i,k), k=j*3-2,j*3)
          enddo
       enddo
       ! (2) Remove the dummy atom contributions
       lim1=idum*3-2
       lim2=idum*3
       do i=1,ncoo+3
          if (i.ge.lim1.and.i.le.lim2) cycle
          do j=1,ncoo+3
             if (j.ge.lim1.and.j.le.lim2) cycle
             if (i.lt.lim1.and.j.lt.lim1) then
                hess(i,j)=dumhess(i,j)
             else if (i.lt.lim1.and.j.ge.lim1) then
                hess(i,j-3)=dumhess(i,j)
             else if (i.ge.lim1.and.j.lt.lim1) then
                hess(i-3,j)=dumhess(i,j)
             else if (i.ge.lim1.and.j.ge.lim1) then
                hess(i-3,j-3)=dumhess(i,j)
             endif
          enddo
       enddo
    else
       ! No dummy atom
       do i=1,ncoo
          do j=1,natm
             read(unit,'(3F20.10)') (hess(i,k), k=j*3-2,j*3)
          enddo
       enddo
    endif

    ! Close the FCMFINAL file
    close(unit)

!----------------------------------------------------------------------
! Calculate the normal mode vectors from the Hessian
!----------------------------------------------------------------------
    call hess2nm(hess)

!----------------------------------------------------------------------
! Consistency check: read the CFOUR log file frequencies and make sure
! that our frequencies match
!----------------------------------------------------------------------    
    call freeunit(unit)
    open(unit,file=trim(freqfile)//'/cfour.log')

5   read(unit,'(a)',end=888) string
    if (index(string,'rotational projection').eq.0) goto 5

    do i=1,7
       read(unit,*)
    enddo
       
    do i=1,nmodes
       read(unit,'(24x,F10.4)') cffreq(i)
       if (abs(cffreq(i)-freq(i)/invcm2ev).gt.1.0d0) then
          errmsg='Mismatch between calculated frequencies:'
          write(errmsg(len_trim(errmsg):),'(2(x,F10.4))') &
               cffreq(i),freq(i)/invcm2ev
          call error_control
       endif
    enddo

    close(unit)

!----------------------------------------------------------------------
! Read the symmetry labels from the CFOUR log file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=trim(freqfile)//'/cfour.log')

10  read(unit,'(a)',end=999) string
    if (index(string,'Normal Coordinate Analysis').eq.0) goto 10

    do i=1,12
       read(unit,*)
    enddo

    do i=1,nmodes
       read(unit,'(8x,a3)') nmlab(i)
    enddo

    close(unit)

    return

888 continue
    errmsg='The frequencies could not be found in: '&
         //trim(freqfile)//'/cfour.log'
    call error_control

999 continue
    errmsg='The Normal Coordinate Analysis section couldn''t &
         be found in: '//trim(freqfile)//'/cfour.log'
    call error_control

  end subroutine getmodes_cfour

!######################################################################

  subroutine hess2nm(hess)

    use constants
    use iomod
    use global
    use utils

    implicit none

    integer                       :: i,j,k,l,e2,error,unit
    integer, dimension(ncoo)      :: indx
    real(d), dimension(ncoo,ncoo) :: hess,proj,phess,eigvec
    real(d), dimension(ncoo)      :: eigval
    real(d), dimension(3*ncoo)    :: work

!----------------------------------------------------------------------
! Mass-weight the Hessian
!----------------------------------------------------------------------
    do i=1,ncoo
       do j=1,ncoo
          hess(i,j)=hess(i,j)/sqrt(mass(i)*mass(j))
       enddo
    enddo

!----------------------------------------------------------------------
! Project out the rotational and translational DOFs
!----------------------------------------------------------------------
    ! Construct the projector onto the complement of the
    ! translation-rotation subspace
    call rtproj(proj)

    ! Project out the rotational and translational DOFs
    phess=matmul(proj,matmul(hess,proj))

!----------------------------------------------------------------------
! Diagonalise the projected mass-weighted Hessian
!----------------------------------------------------------------------
    eigvec=phess
    e2=3*ncoo

    call dsyev('V','U',ncoo,eigvec,ncoo,eigval,work,e2,error)

    if (error.ne.0) then
       errmsg='Diagonalisation of the projected Hessian failed in &
            subroutine hess2nm'
       call error_control
    endif

!----------------------------------------------------------------------
! Sort the normal modes by absolute eigenvalues
!----------------------------------------------------------------------
    eigval=abs(eigval)
    call dsortindxa1('A',ncoo,eigval,indx)

!----------------------------------------------------------------------
! Frequencies in ev
!----------------------------------------------------------------------
    do i=7,ncoo
       freq(i-6)=sqrt(abs(eigval(indx(i))))*5140.8096d0*invcm2ev
    enddo

!-----------------------------------------------------------------------
! Normal mode vectors
! (Scale to mass-weighted x to obtain true normal mode vectors)
!-----------------------------------------------------------------------
    do i=7,ncoo
       nmcoo(:,i-6)=eigvec(:,indx(i))       
    enddo

    do i=1,ncoo
       nmcoo(i,:)=nmcoo(i,:)*sqrt(mass(j))
    enddo

    return

  end subroutine hess2nm

!######################################################################

  subroutine rtproj(proj)

    use constants
    use iomod
    use global

    implicit none

    integer                       :: i,j,k,l
    integer, dimension(6)         :: ipiv
    integer                       :: info
    real(d), dimension(ncoo,ncoo) :: proj,rmat
    real(d), dimension(6,ncoo)    :: rtvec
    real(d), dimension(6,6)       :: smat,invsmat
    real(d), dimension(ncoo,6)    :: bmat
    real(d), dimension(6)         :: work
    logical(kind=4)               :: lcheck

!------------------------------------------------------------------
! Initialise arrays
!------------------------------------------------------------------
    rtvec=0.0d0

!------------------------------------------------------------------
! Vectors 1-3: translation along the three Cartesian axes
!------------------------------------------------------------------
    ! Loop over the translational DOFs
    do i=1,3
       ! Construct the vector for the current DOF
       do j=1,natm
          k=j*3-3+i
          rtvec(i,k)=sqrt(mass(j*3))
       enddo
    enddo

!------------------------------------------------------------------
! Vectors 4-6: infinitesimal displacements corresponding to
!              rotation about the three Cartesian axes
!------------------------------------------------------------------
    ! Rotation about the x-axis
    do i=1,natm
       j=i*3-1
       k=i*3
       rtvec(4,j)=sqrt(mass(i*3))*xcoo0(k)
       rtvec(4,k)=-sqrt(mass(i*3))*xcoo0(j)
    enddo

    ! Rotation about the y-axis
    do i=1,natm
       j=i*3-2
       k=i*3
       rtvec(5,j)=-sqrt(mass(i*3))*xcoo0(k)
       rtvec(5,k)=sqrt(mass(i*3))*xcoo0(j)
    enddo

    ! Rotation about the z-axis
    do i=1,natm
       j=i*3-2
       k=i*3-1
       rtvec(6,j)=sqrt(mass(i*3))*xcoo0(k)
       rtvec(6,k)=-sqrt(mass(i*3))*xcoo0(j)
    enddo

!------------------------------------------------------------------
! Calculate the projector R onto the translational and rotational
! DOFs using R=b*S^-1*b^T, where S=vec^T*vec.
!
! Here, R <-> rmat, S <-> smat, b <-> bmat (matrix of vectors)
!------------------------------------------------------------------
    ! Construct the b-matrix
    bmat=transpose(rtvec)

    ! Calculate the S-matrix
    smat=matmul(transpose(bmat),bmat)

    ! Invert the S-matrix
    invsmat=smat
    call dgetrf(6,6,invsmat,6,ipiv,info)
    if (info.ne.0) then
       errmsg='LU factorisation of the S-matrix failed'
       call error_control
    endif
    call dgetri(6,invsmat,6,ipiv,work,6,info)
    if (info.ne.0) then
       errmsg='Diagonalisation of the S-matrix failed'
       call error_control       
    endif
    
    ! Calculate the projection matrix R <-> rmat
    rmat=matmul(bmat,matmul(invsmat,transpose(bmat)))

!------------------------------------------------------------------
! Construct the projector onto the complement of the
! translation-rotation subspace: P=1-R
!------------------------------------------------------------------
    ! 1
    proj=0.0d0
    do i=1,ncoo
       proj(i,i)=1.0d0
    enddo

    ! 1-R
    proj=proj-rmat

    return

  end subroutine rtproj

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
! freq in ev, mass in amu, x in Angstrom
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
