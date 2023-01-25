module solvers
  use precision
  implicit none
  private

  public :: rk2 
  public :: rk4
  public :: dopri54

  interface
    function func(t,u) result(up)   
      use precision    
      real(dp), intent(in) :: t    
      real(dp), intent(in) :: u(:)    
      real(dp), allocatable :: up(:)    
    end function
  end interface

  contains

  subroutine rk2(f, t, dt, u0, u)
    procedure(func) :: f    
    real(dp), intent(in) :: t
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: u0(:)
    real(dp), intent(inout) :: u(:)

    real(dp), allocatable :: k1(:), k2(:)
 
    k1 = f(t,u0)
    k2 = f(t+dt, u0 + dt*k1)
    u = u0 + (k1+k2)*dt*0.5_dp  

  end subroutine rk2

  subroutine rk4(f, t, dt, u0, u)
    procedure(func) :: f    
    real(dp), intent(in) :: t
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: u0(:)
    real(dp), intent(inout) :: u(:)

    real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:)
 
    k1 = f(t          , u0             )
    k2 = f(t+dt*0.5_dp, u0+dt*k1*0.5_dp)
    k3 = f(t+dt*0.5_dp, u0+dt*k2*0.5_dp)
    k4 = f(t+dt       , u0+dt*k3       )
    u = u0 + (0.16666666666666666666_dp*k1+0.33333333333333333333_dp*k2+ &
          & 0.33333333333333333333_dp*k3+0.16666666666666666666_dp*k4)*dt

  end subroutine rk4
 
  ! --------------------------------------------------
  subroutine dopri54(f, t, h, u0, u, err)
    procedure(func) :: f    
    real(dp), intent(in) :: t, h
    real(dp), dimension(:), intent(in) :: u0 
    real(dp), dimension(:), intent(inout) :: u 
    real(dp), intent(inout) :: err
 
    real(dp), dimension(:,:), allocatable :: k 
 
    allocate(k(size(u0) ,7))
    
    k(:,1) = f(t, u0)*h
    k(:,2) = f(t+h/5.0d0, u0+k(:,1)/5.0d0)*h
    k(:,3) = f(t+3.d0*h/10.0d0, u0+3.d0/40.0d0*k(:,1)+9.d0/40.0d0*k(:,2))*h
    k(:,4) = f(t+h*4.d0/5.d0, u0+44.d0/45.d0*k(:,1)+160.d0/45.d0*k(:,3)-168.d0/45.d0*k(:,2))*h
    k(:,5) = f(t+h*8.d0/9.d0, u0+38744.d0/13122.d0*k(:,1)+128896.d0/13122.d0*k(:,3)-152160.d0/13122.d0*k(:,2) &
            -3816.d0/13122.d0*k(:,4))*h
    k(:,6) = f(t+h, u0+9017.d0/3168.d0*k(:,1)+46732.d0/5247.d0*k(:,3)+5194.d0/18656.d0*k(:,4) &
            -34080.d0/3168.d0*k(:,2)-5103.d0/18656.d0*k(:,5))*h

    u = u0+(35.d0/384.d0*k(:,1)+500.d0/1113.d0*k(:,3)+125.d0/192.d0*k(:,4)+&
              11.d0/84.d0*k(:,6)-2187.d0/6784.d0*k(:,5))

    k(:,7) = f(t+h, u)*h

    ! just estimate the error as  ( z - y ):
    err = maxval(abs( 71.d0/57600.d0*k(:,1)+22.d0/525.d0*k(:,6) &
        & +71.d0/1920.d0*k(:,4) &
        & -71.d0/16695.d0*k(:,3)-17253.d0/339200.d0*k(:,5)-1.d0/40.d0*k(:,7) ))
 
    deallocate(k)
 
  end subroutine dopri54
  ! --------------------------------------------------


  
end module solvers
