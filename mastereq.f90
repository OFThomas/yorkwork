!Oliver Thomas 2016

program mastereq

use matrixfns
implicit none

real(kind=dp1) :: timestep, total_time

!number of states bosonic field
n_b=3

!number of states atom
n_a=2

total_time=1.0_dp1
timestep=1*1e-2_dp1
timesteps=nint(total_time/timestep)
print*, 'Time simulated 		', 'Timestep 		', 'Total steps '
print*, total_time, timestep, timesteps

!Matrix operators

! creation, annihilation, nummatrix
! sigmaz, sigmaplus, sigmaminus

! aident, bident

!-Make operator matrices
 call makeoperators

!initialising density matrix into groundstate-vacuum photon
rho = 0
rho(2,2,1)=1

!generate identities in respective spaces
 aident=identity(n_a)
 bident=identity(n_b)

!make operators only act on their own system, use tensor product of identity in other space
sigmaz=tproduct(bident,sigmaz)
sigmaplus=tproduct(bident,sigmaplus)
sigmaminus=tproduct(bident,sigmaminus)

nummatrix=tproduct(nummatrix,aident)
creation=tproduct(creation, aident)
annihilation=tproduct(annihilation,aident)

!- Print operators to terminal (testing only)

open(unit=11, file='rho.txt', status='replace', iostat=status)
  if (status/=0) stop 'Error in opening rho output file'

t=0
!Increment rho using runge-kutta, save results in rho array
do counter=1,timesteps-1
rho(:,:,counter+1)=rk4(timestep,t,rho(:,:,counter))
end do

write(11,*) rho(:,:,1)
write(11,*) rho(:,:,timesteps)
close(11)

! --------------- End of main program --------------------------------!
end program mastereq
