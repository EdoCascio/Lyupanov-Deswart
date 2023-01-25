module functions
  use precision
  implicit none
  private

  
  public :: deswart
  public :: deswart_linear
  

  real(dp), public :: k1 
  real(dp), public :: k2
  real(dp), public :: Q
  real(dp), public :: A
  real(dp), public :: w
  real(dp), public :: rr
  real(dp), public :: ss
  real(dp), public :: bb
  real(dp), public :: u1,u2,u3,u4,u5,u6
  real(dp), public :: gamma1, gamma2,gamma1s,gamma2s,alfa1, alfa2, beta1,beta2,delta1,delta2, epsilon, us1, b
  real(dp), parameter :: C = 0.1_dp


 
  

contains

  

  function deswart(t, u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    allocate(up(size(u)))

    up(1) = gamma1s*u(3) - (0.1_dp)*(u(1)-us1)
    
    up(2) = -(alfa1*u(1)-beta1)*u(3) - (0.1_dp)*u(2) - delta1*u(4)*u(6)
    
    up(3) = (alfa1*u(1) - beta1)*u(2) - gamma1*u(1) - (0.1_dp)*u(3) + delta1*u(4)*u(5)
    
    up(4) = gamma2s*u(6) - (0.1_dp)*u(4) + epsilon*(u(2)*u(6) - u(3)*u(5))
    
    up(5) = -(alfa2*u(1) - beta2)*u(6) - (0.1_dp)*u(5) - delta2*u(3)*u(4)
    
    up(6) = (alfa2*u(1)-beta2)*u(5) - gamma2*u(4) - (0.1_dp)*u(6) + delta2*u(2)*u(4)


  end function deswart

  function deswart_linear(t, u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)

    allocate(up(size(u)))

    up(1) = -C*u(1) + gamma1s*u(3)
    
    up(2) = -alfa1*u3*u(1) - (0.1_dp)*u(2) - (alfa1*u1 - beta1)*u(3) - delta1*u6*u(4) - delta1*u4*u(6)
    
    up(3) = (alfa1*u2 - gamma1)*u(1) + (alfa1*u1 - beta1)*u(2) - (0.1_dp)*u(3) + delta1*u5*u(4) + delta1*u4*u(5)
    
    up(4) = epsilon*u6*u(2) - epsilon*u5*u(3) - (0.1_dp)*u(4) + epsilon*u3*u(5) + (gamma2s + epsilon*u2)*u(6)
    
    up(5) = -alfa2*u6*u(1)  - delta2*u4*u(3) - delta2*u3*u(4) - (0.1_dp)*u(5) - (alfa2*u1 - beta2)*u(6)
    
    up(6) = alfa2*u5*u(1) + delta2*u4*u(2) + (delta2*u2 - gamma2)*u(4) + (alfa2*u1 - beta2)*u(5) - (0.1_dp)*u(6)


  end function deswart_linear


end module functions


