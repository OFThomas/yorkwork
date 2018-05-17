module runge
implicit none

integer, parameter :: dp=selected_real_kind(15,300)
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
  
  real(kind=dp) :: e, h, x1, z1, y1, x2, y2, z2, l,emin,emax,z 
  real(kind=dp) :: dz, zleft, zright, yleft, yright, v,dzr
  integer :: n, i, count, final, number, status
  
  !Varaible to change, e is always > 0
  emin=0.0_dp
  emax=1.0_dp
  !odd or even n state 
  number= 1

  e=0.5_dp*(emin+emax)

  v=100000000.0_dp
  l=4.0_dp
  z = Sqrt(2.0_dp/l)
  
  open(unit=15, file='data.dat', status='replace',iostat=status)
  if (status/=0) stop 'Error in oppening left file'
  open(unit=20, file='data1.dat', status='replace', iostat=status)
  if (status/=0) stop'Error in opening right file'  
  open(unit=25, file='dlsleft.dat', status='replace', iostat=status)
  if (status/=0) stop 'Error opening dls file'  
  open(unit=30, file='dlsright.dat', status='replace', iostat=status)
  if (status/=0) stop 'Error opening dls file'

  final =0
  n=100000

main:  do count=1, 1000

!right loop
   h=-0.001_dp
   x1=l
   y1=0.0_dp
 !changes gradient for odd solutions
   if (mod(number,2) /= 0) then 
      z1=-z
   else
      z1=z
   end if
   !outside potential is zero
      right: do i=1,n
        if (x1 >= 2.0_dp) then
           y1=0.0_dp
        else
           v=0
        end if
        !updating x, y, z
        if (x1 <= 0) then
           yright=y1
           zright=z1
           exit right
        end if
        !update position
        call r2k(h,x1,y1,z1,f1,f2,e)
        !check if e is within errror
        !if true print final trajectory
        if (final==1) then
           write (15,*) x1,y1
           write (25,*) x1, (y1/z1)
        end if
     end do right

!left loop
     h=0.001_dp
     x2=-l
     y2=0.0_dp
     z2=z
     left: do i=1,n
        !potential check again
        if (abs(x2) >= 2.0_dp) then
           y2=0.0_dp
        else
           v=0
        end if
        !updating x,y,z
        if(x2 >= 0) then
           yleft=y2
           zleft=z2
           exit left
        end if
        call r2k(h,x2,y2,z2,f1,f2,e)
        !check if e is within errror
        !if true print final trajectory
        if (final ==1) then
           write (20,*) x2,y2
           write(30,*) x2, (z2/y2)
        end if
     end do left
       !Exit loop with e in range
     if (final==1) then
        exit main
     end if

     !check if e is within range
        if (final==0 ) then
           !Difference for odd n
           if (mod(number,2) /= 0) then
              dz = (zleft/yleft) -  (zright/yright)
              if (dz > 0) then
                 emin = e
                 e = 0.5_dp*(emin + emax)
              else if (dz < 0) then 
                 emax=e
                 e = 0.5_dp*(emax + emin)
              end if
              if ((abs(dz) < 0.00000001_dp)) then
                 final = 1
              end if
              print*, 'e =', e
         
              !Dz for even n, reciprical
           else
              dzr = (yleft/zleft) - ( (yright/zright))
              if (dzr < 0) then
                 emin = e
                 e=0.5_dp*(emin+ emax)
              else if (dzr > 0) then
                 emax = e
                 e=0.5_dp*(emin + emax)
              end if
              if ((abs(dzr) < 0.00000001_dp)) then 
                 final =1
              end if
              print*, 'e=', e
           end if
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
    real (kind=dp), intent(in) :: t,y,z,e
    f2=(v-e)*y
  end function f2

end program schro
