!######################################################################
! transw0: a program to transform the zeroth-order potential from one
!          set of normal modes to another. This simply involves a
!          shift of the point of expansion followed by a rotation of
!          the normal modes.
!######################################################################

program transw0
  
  implicit none

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
  call initialise

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
  call rdinput

!----------------------------------------------------------------------
! Determine the type or types of quantum chemistry calculations used
! for the two frequency calculations
!----------------------------------------------------------------------
  call get_freqtypes

!----------------------------------------------------------------------
! Determine the number of atoms and set the number of normal modes
! Note that we assume a non-linear molecule here.
!----------------------------------------------------------------------
  call set_dimensions
  
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  call alloc

!----------------------------------------------------------------------
! Read the frequency files to get the normal modes and reference
! geometries
!----------------------------------------------------------------------
  call get_coords

!----------------------------------------------------------------------
! Create the transformation matrices
!----------------------------------------------------------------------
  call get_transmat

!----------------------------------------------------------------------
! Put everything into the Eckart frame
!----------------------------------------------------------------------
  call eckart

!----------------------------------------------------------------------
! Calculate the shift DeltaQ in terms of normal modes A between
! reference point A and reference point B 
!----------------------------------------------------------------------
  call get_DeltaQ
  
!----------------------------------------------------------------------
! Calculate the normal mode A to normal mode B transformation matrix
!----------------------------------------------------------------------
  call get_overlap

!----------------------------------------------------------------------
! Calculate the transformed zeroth-order potential
!----------------------------------------------------------------------
  call get_transw0

!----------------------------------------------------------------------
! For testing purposes, output the minimum geometry in terms of the
! normal modes B
!----------------------------------------------------------------------
  call get_min_QB
  
!----------------------------------------------------------------------
! Output the transformed zeroth-order potential
!----------------------------------------------------------------------
  call wrtransw0
  
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
! Open the log file
!----------------------------------------------------------------------
    call freeunit(ilog)
    open(ilog,file='transw0.log',form='formatted',status='unknown')
    
    return
    
  end subroutine initialise

!######################################################################

  subroutine rdinput

    use constants
    use global
    use transmod
    use iomod
    use channels
    
    implicit none

    integer            :: n
    character(len=120) :: string1,string2

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ! Name of the frequency files
    freqfileA=''
    freqfileB=''

!----------------------------------------------------------------------
! Read the command line arguments
!----------------------------------------------------------------------
    n=0

5   continue
    n=n+1
    call getarg(n,string1)

    if (string1.eq.'-fa'.or.string1.eq.'-fA') then
       n=n+1
       call getarg(n,freqfileA)
    else if (string1.eq.'-fb'.or.string1.eq.'-fB') then
       n=n+1
       call getarg(n,freqfileB)
    else
       errmsg='Unknown keyword: '//trim(string1)
       call error_control
    endif

    if (n.lt.iargc()) goto 5

!----------------------------------------------------------------------
! Check that all the required information has been given
!----------------------------------------------------------------------
    ! Frequency file A
    if (freqfileA.eq.'') then
       errmsg='Frequency file/directory A has not been given'
       call error_control
    endif

    ! Frequency file B
    if (freqfileB.eq.'') then
       errmsg='Frequency file/directory B has not been given'
       call error_control
    endif
    
    return
    
  end subroutine rdinput
    
!######################################################################

  subroutine get_freqtypes

    use constants
    use global
    use transmod
    use ioqc
    
    implicit none

!----------------------------------------------------------------------
! Frequency calculation A
!----------------------------------------------------------------------
    freqfile=freqfileA

    call freqtype

    freqtypA=freqtyp

!----------------------------------------------------------------------
! Frequency calculation B
!----------------------------------------------------------------------
    freqfile=freqfileB

    call freqtype

    freqtypB=freqtyp

    return
    
  end subroutine get_freqtypes

!######################################################################

  subroutine set_dimensions

    use constants
    use global
    use transmod
    use ioqc
    
    implicit none

!----------------------------------------------------------------------
! Determine the number of atoms from one of the frequency files
!----------------------------------------------------------------------
    freqtyp=freqtypA
    freqfile=freqfileA
    call getnatm_freqfile

!----------------------------------------------------------------------
! Set the number of coordinates
! N.B. We are assuming a non-linear molecule here...
!----------------------------------------------------------------------
    ncoo=3*natm
    nmodes=ncoo-6
    
    return
    
  end subroutine set_dimensions
    
!######################################################################

  subroutine alloc

    use transmod
    use global
    
    implicit none

    ! Atomic masses
    allocate(mass(ncoo))
    mass=0.0d0

    ! Reference coordinates
    allocate(xcoo0(ncoo))
    allocate(xcoo0A(ncoo))
    allocate(xcoo0B(ncoo))
    xcoo0=0.0d0
    xcoo0A=0.0d0
    xcoo0B=0.0d0
    
    ! Atomic numbers
    allocate(atnum(natm))
    atnum=0

    ! Atom labels
    allocate(atlbl(natm))
    atlbl=''

    ! Normal mode labels
    allocate(nmlab(nmodes))
    allocate(nmlabA(nmodes))
    allocate(nmlabB(nmodes))
    nmlab=''
    nmlabA=''
    nmlabB=''
    
    ! Normal mode frequencies
    allocate(freq(nmodes))
    allocate(freqA(nmodes))
    allocate(freqB(nmodes))
    freq=0.0d0
    freqA=0.0d0
    freqB=0.0d0
    
    ! Transformation matrices
    allocate(nmcoo(ncoo,nmodes))
    allocate(coonm(nmodes,ncoo))
    allocate(nmcooA(ncoo,nmodes))
    allocate(coonmA(nmodes,ncoo))
    allocate(nmcooB(ncoo,nmodes))
    allocate(coonmB(nmodes,ncoo))
    nmcoo=0.0d0
    coonm=0.0d0
    nmcooA=0.0d0
    coonmA=0.0d0
    nmcooB=0.0d0
    coonmB=0.0d0

    ! Normal mode shift
    allocate(DeltaQ(nmodes))
    DeltaQ=0.0d0

    ! Rotation matrix
    allocate(Smatrix(nmodes,nmodes))
    Smatrix=0.0d0

    ! Transformed zeroth-order potential
    allocate(kappaB(nmodes))
    allocate(gammaB(nmodes,nmodes))
    kappaB=0.0d0
    gammaB=0.0d0
    
    return
    
  end subroutine alloc
  
!######################################################################

  subroutine finalise

    use channels
    use transmod
    use global

    implicit none

!----------------------------------------------------------------------
! Close the logfile
!----------------------------------------------------------------------
    close(ilog)

!----------------------------------------------------------------------
! Deallocation
!----------------------------------------------------------------------
    deallocate(mass)
    deallocate(xcoo0)
    deallocate(xcoo0A)
    deallocate(xcoo0B)
    deallocate(atnum)
    deallocate(atlbl)
    deallocate(nmlab)
    deallocate(nmlabA)
    deallocate(nmlabB)
    deallocate(freq)
    deallocate(freqA)
    deallocate(freqB)
    deallocate(nmcoo)
    deallocate(coonm)
    deallocate(nmcooA)
    deallocate(coonmA)
    deallocate(nmcooB)
    deallocate(coonmB)
    deallocate(DeltaQ)
    deallocate(Smatrix)
    deallocate(kappaB)
    deallocate(gammaB)
    
    return
    
  end subroutine finalise
    
!######################################################################

  subroutine get_coords

    use constants
    use global
    use transmod
    use ioqc
    
    implicit none

!----------------------------------------------------------------------
! Normal modes and reference geometry for point A
!----------------------------------------------------------------------
    ! Read the frequency file
    freqfile=freqfileA
    freqtyp=freqtypA
    call rdfreqfile

    ! Copy arrays
    xcoo0A=xcoo0
    nmlabA=nmlab
    freqA=freq
    nmcooA=nmcoo
    coonmA=coonm

    ! Rename the modes.xyz file
    call system('mv modes.xyz modesA.xyz')
    
!----------------------------------------------------------------------
! Normal modes and reference geometry for point B
!----------------------------------------------------------------------
    ! Read the frequency file
    freqfile=freqfileB
    freqtyp=freqtypB
    call rdfreqfile

    ! Copy arrays
    xcoo0B=xcoo0
    nmlabB=nmlab
    freqB=freq
    nmcooB=nmcoo
    coonmB=coonm

    ! Rename the modes.xyz file
    call system('mv modes.xyz modesB.xyz')
    
    return
    
  end subroutine get_coords

!######################################################################

  subroutine get_transmat

    use constants
    use global
    use transmod
    use ioqc

    implicit none

!----------------------------------------------------------------------
! Transformation matrices for geometry A
!----------------------------------------------------------------------
    nmcoo=nmcooA
    freq=freqA
    call nm2xmat

    nmcooA=nmcoo
    coonmA=coonm

!----------------------------------------------------------------------
! Transformation matrices for geometry B
!----------------------------------------------------------------------
    nmcoo=nmcooB
    freq=freqB
    call nm2xmat

    nmcooB=nmcoo
    coonmB=coonm
    
    return
    
  end subroutine get_transmat

!######################################################################

  subroutine eckart

    use constants
    use iomod
    use transmod
    use global
    
    implicit none

    integer                 :: i,j
    real(d), dimension(3,3) :: iteigA,iteigB,P
    
!----------------------------------------------------------------------
! Eckart frame transformation for Cartesian coordinates A
!----------------------------------------------------------------------
    call eckart_transmat(xcoo0A,iteigA)

!----------------------------------------------------------------------
! Eckart frame transformation for Cartesian coordinates B
!----------------------------------------------------------------------
    call eckart_transmat(xcoo0B,iteigB)

    ! TEMPORARY: only translate
    !return
    ! TEMPORARY: only translate
    
!----------------------------------------------------------------------
! Take the +/- combinations of the x-, y- and z-directions that puts
! the rotated xcoo0A into maximum coincidence with xcoo0B
!----------------------------------------------------------------------
    ! Determine the transformation that puts xcoo0A and xcoo0B into
    ! maximum coincidence
    call minrmsdAB(xcoo0A,xcoo0B,P)

    ! Transform xcoo0A
    do i=1,natm
       xcoo0A(i*3-2:i*3)=matmul(P,xcoo0A(i*3-2:i*3))
    enddo

    ! Update iteigA
    iteigA=matmul(P,iteigA)
    
!----------------------------------------------------------------------
! Transformation of the normal modes A
!----------------------------------------------------------------------
    do j=1,nmodes
       do i=1,natm
          nmcooA(i*3-2:i*3,j)=matmul(iteigA,nmcooA(i*3-2:i*3,j))
          coonmA(j,i*3-2:i*3)=matmul(iteigA,coonmA(j,i*3-2:i*3))
       enddo
    enddo
    
!----------------------------------------------------------------------
! Transformation of the normal modes B
!----------------------------------------------------------------------
    do j=1,nmodes
       do i=1,natm
          nmcooB(i*3-2:i*3,j)=matmul(iteigB,nmcooB(i*3-2:i*3,j))
          coonmB(j,i*3-2:i*3)=matmul(iteigB,coonmB(j,i*3-2:i*3))
       enddo
    enddo
    
    return
    
  end subroutine eckart

!######################################################################

  subroutine eckart_transmat(xcoo,iteig)

    use constants
    use iomod
    use transmod
    use global
    
    implicit none

    integer                  :: i,j,k,error
    real(d), dimension(ncoo) :: xcoo
    real(d), dimension(3,3)  :: iteig
    real(d), dimension(3)    :: com,eigval
    real(d), dimension(9)    :: work
    
!----------------------------------------------------------------------
! Calculate the centre of mass
!----------------------------------------------------------------------
    ! Centre of mass
    com=0.0d0
    do j=1,3
       do i=1,natm
          com(j)=com(j)+(mass(i*3-3+j)*xcoo(i*3-3+j))
       enddo
    enddo
    com=com/(sum(mass)/3)

!-----------------------------------------------------------------------
! Translate the coordinates to the centre of mass
!-----------------------------------------------------------------------
    do i=1,natm
       do j=1,3
          xcoo(i*3-3+j)=xcoo(i*3-3+j)-com(j)
       enddo
    enddo


    ! TEMPORARY: only translate
    !return
    ! TEMPORARY: only translate
    
!-----------------------------------------------------------------------
! Calculate the moment of inertia tensor
!-----------------------------------------------------------------------
    iteig=0.0d0
    do i=1,natm
       iteig(1,1)=iteig(1,1)+mass(i*3)*(xcoo(i*3-1)**2+xcoo(i*3)**2)
       iteig(2,2)=iteig(2,2)+mass(i*3)*(xcoo(i*3-2)**2+xcoo(i*3)**2)
       iteig(3,3)=iteig(3,3)+mass(i*3)*(xcoo(i*3-2)**2+xcoo(i*3-1)**2)
       iteig(1,2)=iteig(1,2)-mass(i*3)*(xcoo(i*3-2)*xcoo(i*3-1))
       iteig(1,3)=iteig(1,3)-mass(i*3)*(xcoo(i*3-2)*xcoo(i*3))
       iteig(2,3)=iteig(2,3)-mass(i*3)*(xcoo(i*3-1)*xcoo(i*3))
    enddo
    iteig(2,1)=iteig(1,2)
    iteig(3,1)=iteig(1,3)
    iteig(3,2)=iteig(2,3)

!-----------------------------------------------------------------------
! Diagonalise the moment of inertia tensor
!-----------------------------------------------------------------------
    call dsyev('V','U',3,iteig,3,eigval,work,9,error)
    if (error.ne.0) then
       errmsg='Diagonalisation of the moment of inertia tensor failed'
       call error_control
    endif

!-----------------------------------------------------------------------
! Rotate the coordinates
!-----------------------------------------------------------------------
    do i=1,natm
       xcoo(i*3-2:i*3)=matmul(iteig,xcoo(i*3-2:i*3))
    enddo
    
    return
    
  end subroutine eckart_transmat

!######################################################################

  subroutine minrmsdAB(xcooA,xcooB,P)

    use constants
    use iomod
    use transmod
    use global
    
    implicit none

    integer                  :: ix,iy,iz,k,l
    real(d), dimension(ncoo) :: xcooA,xcooB,tmp
    real(d), dimension(3,3)  :: P,P1
    real(d)                  :: rmsd,low

    low=1e+6_d

    do ix=1,2
       do iy=1,2
          do iz=1,2

             P1=0.0d0
             P1(1,1)=(-1)**ix
             P1(2,2)=(-1)**iy
             P1(3,3)=(-1)**iz
             do k=1,natm
                tmp(k*3-2:k*3)=matmul(P1,xcooA(k*3-2:k*3))
             enddo

             rmsd=sqrt(dot_product(xcooB-tmp,xcooB-tmp))

             if (rmsd.lt.low) then
                low=rmsd
                P=P1
             endif
             
          enddo
       enddo
    enddo

    return
    
  end subroutine minrmsdAB

!######################################################################

  subroutine get_DeltaQ

    use constants
    use iomod
    use transmod
    use global
    
    implicit none

    integer                  :: i,j
    real(d), dimension(ncoo) :: DeltaX,x

!----------------------------------------------------------------------
! Calculate the DeltaQ in terms of normal modes A between reference
! point A and reference point B 
!----------------------------------------------------------------------
    DeltaX=(xcoo0B-xcoo0A)/ang2bohr
    DeltaQ=matmul(coonmA,DeltaX)
        
    return
    
  end subroutine get_DeltaQ

!######################################################################

  subroutine get_overlap

    use constants
    use iomod
    use transmod
    use global

    implicit none

    Smatrix=matmul(coonmA,nmcooB)

    return
    
  end subroutine get_overlap

!######################################################################

  subroutine get_transw0

    use constants
    use iomod
    use transmod
    use global

    implicit none

    integer :: i,j,k
    
!----------------------------------------------------------------------
! First-order terms
!----------------------------------------------------------------------
    kappaB=0.0d0
    do i=1,nmodes
       do j=1,nmodes
          kappaB(i)=kappaB(i)+Smatrix(j,i)*freqA(j)*DeltaQ(j)
       enddo
    enddo

!----------------------------------------------------------------------
! Second-order terms
!----------------------------------------------------------------------
    gammaB=0.0d0
    do i=1,nmodes
       do j=1,nmodes
          do k=1,nmodes
             gammaB(i,j)=gammaB(i,j)+Smatrix(k,i)*Smatrix(k,j)*freqA(k)
          enddo
       enddo
    enddo
    
    return
    
  end subroutine get_transw0

!######################################################################

  subroutine get_min_QB

    use constants
    use iomod
    use channels
    use transmod
    use global
    
    implicit none

    integer                           :: i,j,info
    integer, dimension(nmodes)        :: ipiv
    real(d), dimension(nmodes)        :: Q0B
    real(d), dimension(ncoo)          :: X0
    real(d), dimension(nmodes,nmodes) :: tmp
    real(d), dimension(nmodes)        :: gradient
    
!----------------------------------------------------------------------
! Minumum point on the zeroth-order potential in terms of the normal
! modes B
!----------------------------------------------------------------------
    ! LU factorisation of gammaB
    tmp=gammaB
    call dgetrf(nmodes,nmodes,tmp,nmodes,ipiv,info)
    if (info.ne.0) then
       errmsg='LU decomposition of gammaB failed'
       call error_control
    endif

    ! Potential mimumum in terms of the normal modes B
    Q0B=-kappaB
    call dgetrs('N',nmodes,1,tmp,nmodes,ipiv,Q0B,nmodes,info)
    if (info.ne.0) then
       errmsg='Failed call to zgetrs'
       call error_control
    endif

!----------------------------------------------------------------------
! Conversion to Cartesian coordinates
!----------------------------------------------------------------------
    X0=(xcoo0B/ang2bohr)+matmul(nmcooB,Q0B)

!----------------------------------------------------------------------
! Output some information to the log file
!----------------------------------------------------------------------
    write(ilog,'(/,50a)') ('*',i=1,50)
    write(ilog,'(3x,a)') 'Minimum geometry in terms of normal modes B'
    write(ilog,'(50a)') ('*',i=1,50)
    do i=1,nmodes
       write(ilog,'(i2,2x,F10.7)') i,Q0B(i)
    enddo
    
    !write(6,'(i1,/)') natm
    !do i=1,natm
    !   write(6,'(a1,3(2x,F12.7))') atlbl(i),(X0(j),j=i*3-2,i*3)
    !enddo
    !
    !write(6,'(i1,/)') natm
    !do i=1,natm
    !   write(6,'(a1,3(2x,F12.7))') atlbl(i),(xcoo0A(j)/ang2bohr,j=i*3-2,i*3)
    !enddo
    !
    !stop
    
    return
    
  end subroutine get_min_QB
    
!######################################################################
  
  subroutine wrtransw0

    use constants
    use iomod
    use channels
    use transmod
    use global

    implicit none

    integer            :: i,m1,m2
    real(d), parameter :: thrsh=1e-5
    real(d)            :: par
    character(len=2)   :: am1,am2

!----------------------------------------------------------------------
! Section header
!----------------------------------------------------------------------
    write(ilog,'(/,50a)') ('*',i=1,50)
    write(ilog,'(8x,a)') 'Transformed zeroth-order potential'
    write(ilog,'(50a)') ('*',i=1,50)
    
!----------------------------------------------------------------------
! First-order terms
!----------------------------------------------------------------------
    write(ilog,'(/,a)') '# 1st-order transformed zeroth-order &
         potential terms'
    do m1=1,nmodes
       if (abs(kappaB(m1)).lt.thrsh) cycle
       write(am1,'(i2)') m1
       write(ilog,'(a,F9.6,a)') 'kappaB'//'_'//trim(adjustl(am1))// &
            ' = ',kappaB(m1),' , ev'
    enddo

!----------------------------------------------------------------------
! Second-order terms
!----------------------------------------------------------------------
    write(ilog,'(/,a)') '# 2nd-order transformed zeroth-order &
         potential terms'

    ! Loop over the unique pairs of modes
    do m1=1,nmodes
       do m2=m1,nmodes

          ! Take account of the fact that we will only write the
          ! unique pairs to the operator file
          if (m1.eq.m2) then
             par=gammaB(m1,m2)
          else
             par=2.0d0*gammaB(m1,m2)
          endif
          
          if (abs(par).lt.thrsh) cycle

          write(am1,'(i2)') m1
          write(am2,'(i2)') m2
          write(ilog,'(a,F9.6,a)') 'gammaB'//'_'&
               //trim(adjustl(am1))//'_'//trim(adjustl(am2))&
               //' = ',par,' , ev'
          
       enddo
    enddo
          
    return
    
  end subroutine wrtransw0
  
!######################################################################
  
end program transw0
