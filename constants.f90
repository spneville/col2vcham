module constants
  
  save

  integer, parameter    :: d=selected_real_kind(8)
  integer, parameter    :: lng=selected_int_kind(16)
  
  real(d), parameter    :: rzero=0._d
  real(d), parameter    :: rone=1._d
  real(d), parameter    :: pi=3.14159265358979_d
  real(d), parameter    :: ang2bohr=1.88972612d0
  real(d), parameter    :: invcm2ev=1.23985e-4_d
  complex(d), parameter :: ci=(0._d,1._d)
  complex(d), parameter :: czero=(0._d,0._d)
  complex(d), parameter :: cone=(1._d,0._d)
  
end module constants
