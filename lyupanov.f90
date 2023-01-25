program lyapunov 
 use precision
 use functions
 use solvers
 use gram_schmidt
 implicit none

 real(dp) :: h, t, tfin, abserr, error, s
 real(dp), dimension(6) :: u, unext, lambda, cum
 real(dp), dimension(6,6) :: VO, VA
 integer :: i, k, Nstep, Niters 
 character(10) :: arg, fname
 real(dp), parameter :: pi = acos(-1.0_dp)

  alfa1 = (8.0*sqrt(2.0_dp)/pi)*(1/3)*((b*b)/(b*b + 1))
  alfa2 = (8.0*sqrt(2.0_dp)/pi)*(4/15)*((b*b + 3)/(b*b + 4))

  delta1 = (64.0*sqrt(2.0_dp)/(15*pi))*((b*b)/(b*b + 1))
  delta2 = (64.0*sqrt(2.0_dp)/(15*pi))*((b*b - 3)/(b*b + 4))

  epsilon = 16*sqrt(2.0_dp)/(5.0*pi)

  beta1 = (1.25*b*b)/(b*b + 1)
  beta2 = (1.25*b*b)/(b*b + 4)

  gamma1s = 4/3 * (sqrt(2.0)*b*1)/(pi)
  gamma2s = 8/15 * (sqrt(2.0)*b*1)/(pi)

  gamma1 = 4/3 * (sqrt(2.0)*b*1)/(pi*(b*b + 1))
  gamma2 = 32/15 * (sqrt(2.0)*b*1)/(pi*(b*b + 4))
 
