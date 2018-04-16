
module runge
implicit none

integer, parameter :: dp=selected_real_kind(15,300)
real(kind=dp), parameter, public :: pi=3.141592653589793238462643383279502884197_dp
real(kind=dp) :: k_1, k_11, k_2, k_22 

contains

!Runge-kutta Sub
  subroutine r2k(h, t, y, z,f1,f2,e)
    real (kind=dp) :: t, y, z, h, e
    real (kind=dp) :: f1, f2
    intent (in) :: h, e
    k_11 = (h)*f1(t,y,z)
    k_22 = (h)*f2(t,y,z, e)
    
    K_1 = (h)*f1(t+h, y+k_11, z+k_22)  
    k_2 = (h)*f2(t+h, y+k_11, z+k_22, e)
    
    y = y + 0.5*(k_11 + k_1)
    z = z + 0.5*(k_22 + k_2)
    t=t+h
  end subroutine r2k

end module runge


program schro
  use runge
  implicit none
  
  real(kind=dp) :: e, h, x, z, y,l,emin,emax, v_0 
  real(kind=dp) :: a,b,c, dz, zleft, zright, yleft, yright, v
  integer :: n, i, count, final
  
  !e=2.5_dp

 emin=5.0_dp
!Varaible to change, e is always > 0
  emax=6.0_dp
  e=0.5_dp*(emin+emax)

  v=10000.0_dp
  l=4.0_dp
  z = Sqrt(2.0_dp/l)

  open(unit=15, file='data.dat', status='replace')
  open(unit=20, file='data1.dat', status='replace')
  !when final =true stop loop. Solved
  final = 0
  n=1000

main:  do count=1, 100

   h=0.01_dp
   x=0.0_dp
   y=0.0_dp
   z=z
   a=2.0_dp
   v_0=10.0_dp
   call r2k(h,x,y,z,f1,f2,e)
   if (final ==1) then
      write (20,*) x,y
   end if
   
   if (final==1) then
      exit main
   end if
   
   else
      dz = y
      print*, 'diff =', dz
      if (dz > 0) then
         emin = e
         e = 0.5*(emin + emax)
      else if (dz < 0) then 
         emax=e
         e = 0.5*(emax + emin)
      end if
      print*, 'emax =', emax, 'emin =', emin
      if ((abs(dz) < 0.01)) then
         !exit main
         final = 1
      end if
      print*, 'e =', e
   end if
   
end do main

close(unit=15)
close(unit=20)

contains
  !z=d(psi)/dx
  function f1(t,y,z) 
    real (kind=dp) :: f1
    real (kind=dp), intent(in) :: t,y,z
    f1=z
  end function f1

  !d2(psi)/dx2
  function f2(t,y,z,e)
    real (kind=dp) :: f2
    real (kind=dp), intent(in) :: t,y,z,e,x,a
    f2=(v_0*(x/a)-e)*y
  end function f2
  
end program schro
