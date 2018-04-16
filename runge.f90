
module runge
implicit none

integer, parameter :: dp=selected_real_kind(15,300)
real(kind=dp) :: k_1, k_11, k_2, k_22 

contains

!Runge-kutta Sub
  subroutine r2k(h, t, y, z,f1,f2)
    real (kind=dp) :: t, y, z, h
    real (kind=dp) :: f1, f2
    intent (in) :: h
    k_11 =h*f1(t,y,z)
    k_22 = h*f2(t,y,z)
    
    K_1 = h*f1(t+h, y+k_11, z+k_22)  
    k_2 = h*f2(t+h, y+k_11, z+k_22)
    
    y = y + 0.5*(k_11 + k_1)
    z = z + 0.5*(k_22 + k_2)
    t=t+h
  end subroutine r2k

end module runge

