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
  
   if(iargc()<5) then
  print*,'lyapunov h Nstep Niters b us1'  
  stop  
 endif

 call getarg(1,arg)
 read(arg,*) h
 
 call getarg(2,arg)
 read(arg,*) Nstep
 
 call getarg(3,arg)
 read(arg,*) Niters 

 call getarg(4,arg)
 read(arg,*) b
 
 call getarg(5,arg)
 read(arg,*) us1
   
 !call getarg(6,arg)
 !read(arg,*) bb
 
 ! ------------------------------------------------------------------
 open(194, file='values.dat')

 t=0.0_dp
 u(1) = 1.0_dp
 u(2) = 0.0_dp
 u(3) = 0.0_dp
 u(4) = 0.0_dp
 u(5) = 0.0_dp
 u(6) = 0.0_dp
 s = 1.0_dp

 tfin = 0.0_dp
 VO = 0.0_dp
 VO(1,1) = 1.0_dp 
 VO(2,2) = 1.0_dp 
 VO(3,3) = 1.0_dp
 VO(4,4) = 1.0_dp
 VO(5,5) = 1.0_dp
 VO(6,6) = 1.0_dp
 cum = 0.0_dp
 
 ! remove transient
 do k = 1, Niters*Nstep 
   call dopri54(deswart, t, h, u, unext, abserr)
   t = t + h
   u = unext
 end do

 do k = 1, Niters ! numero di iterazioni tra una call di gs e l'altra 
    !print*,"iter:",k
    VA=VO
    ! Propago per un po' di step
    do i = 1, Nstep
       call dopri54(deswart, t, h, u, unext, abserr)
       u = unext
       !print*,t,y,abserr
       u1=u(1)
       u2=u(2)
       u3=u(3)
       u4=u(4)
       u5=u(5)
       u6=u(6)
       call dopri54(deswart_linear, t, h, VA(:,1), VA(:,1), abserr)
       call dopri54(deswart_linear, t, h, VA(:,2), VA(:,2), abserr)
       call dopri54(deswart_linear, t, h, VA(:,3), VA(:,3), abserr)
       call dopri54(deswart_linear, t, h, VA(:,4), VA(:,4), abserr)
       call dopri54(deswart_linear, t, h, VA(:,5), VA(:,5), abserr)
       call dopri54(deswart_linear, t, h, VA(:,6), VA(:,6), abserr)
       
       t = t + h
    end do
    
    call gs(VA,VO,lambda)
    
    tfin = tfin + Nstep*h
    if (any(lambda<0)) then
        stop 'lambda<0: reduce Nstep'
    elseif (any(isnan(lambda))) then
        stop 'lambda=NaN: reduce Nstep'
    else
        cum = cum + log(lambda)
    end if
    write(194,*) tfin, cum/tfin

 end do

 close(194)

  !Faccio la media sugli ultimo 1/3 del campione (convergenza)
 open(194,file='values.dat')
 do i = 1, int(2*Niters/3)
   read(194,*) tfin, cum
end do  
 lambda=0.0_dp
 do i = int(2*Niters/3)+1, Niters
   read(194,*) tfin, cum
   lambda=lambda+cum
 end do
 lambda(1) = lambda(1)/(real(Niters,dp)/3.0_dp)
 lambda(2) = lambda(2)/(real(Niters,dp)/3.0_dp)
 lambda(3) = lambda(3)/(real(Niters,dp)/3.0_dp)
 lambda(4) = lambda(4)/(real(Niters,dp)/3.0_dp)
 lambda(5) = lambda(5)/(real(Niters,dp)/3.0_dp)
 lambda(6) = lambda(6)/(real(Niters,dp)/3.0_dp)
 print*,'lyapunov exponents: '
 print*, lambda(1)
 print*, lambda(2)
 print*, lambda(3)
 print*, lambda(4)
 print*, lambda(5)
 print*, lambda(6)
 print*, sum(lambda), -(0.1_dp)*(6.0_dp)
 !print*, "D_kap_yor=", 2.d0+(lambda(3)+lambda(2))/abs(lambda(1))
  
 close(194)

end program lyapunov 

 
