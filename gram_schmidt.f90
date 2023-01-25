module gram_schmidt
  use precision
  implicit none
  private

  public :: gs
    ! questa interface serve per chiamare con lo stesso nome due subroutine
    ! diverse a seconda del tipo di argomenti che si mettono nella chiamata
  interface gs
    module procedure gs_1
    module procedure gs_2
  end interface gs

  contains

  subroutine gs_1(U,y,n)
    real(dp), dimension(:,:), intent(in) :: U
    real(dp), dimension(:), intent(inout) :: y
    integer, intent(in) :: n

    integer :: i
    real(dp) :: norm
    
    do i = 1, n

       y(:) = y(:) - dot_product(U(:,i),y(:))*U(:,i)
       
    end do
    
    norm = sqrt(dot_product(y,y))
    
    if (norm > 1e-10) then
       y(:) = y(:) / norm 
    else
      write(*,*) "ERROR"
      stop
    endif  

  end subroutine gs_1


  subroutine gs_2(V,U,norm)
    real(dp), dimension(:,:), intent(in) :: V
    real(dp), dimension(:,:), intent(out) :: U
    real(dp), dimension(:), intent(out) :: norm
  
    integer :: k, j
    real(dp), dimension(:), allocatable :: x
    real(dp), dimension(:,:), allocatable :: R
  
    allocate(x(size(V,1)))
    allocate(R(size(V,2),size(V,2)))
    R = 0.0_dp
 
    ! G.S: 
    ! u_1 = v_1
    ! u_2 = v_2 - u_1 u_1*v_2
    ! u_3 = v_3 - u_1 u_1*v_2 - u_2 u_2*v_3

    !do k= 1, size(V,2)
     ! x(:) = V(:,k) 
      !do j= 1, k-1
       ! R(j,k) = dot_product(U(:,j),x) 
        !x(:) = x(:) - R(j,k) * U(:,j)
      !end do
      norm(k) = sqrt(dot_product(x,x))
     ! U(:,k) = x(:)/norm(k)
     ! R(k,k) = norm(k)
   ! end do  

    do k= size(V,2), 1, -1 ! size(V,2) è il numero di vettori
      x(:) = V(:,k) ! x è il vettore k-esimo in V
      do j= k+1, size(U,2)
        R(j,k) = dot_product(U(:,j),x) ! prod. tra j-esimo vett. di U e k-esimo di V
        x(:) = x(:) - R(j,k) * U(:,j) ! comb. lin. dei vettori
      end do
      norm(k) = sqrt(dot_product(x,x))
      U(:,k) = x(:)/norm(k)
      R(k,k) = norm(k)
    end do  
  
  
    deallocate(x)
    deallocate(R)
  
  end subroutine gs_2
  

end module gram_schmidt


