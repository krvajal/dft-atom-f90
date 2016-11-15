module dft
    use double
    use ode, only: numerov
    use quad, only: trapz
    implicit none
        
    interface
         function callback(x, params)
            
            use double

            real(kind=dp), intent(in) :: x, params(:)
            real(kind=dp) :: callback
        end function callback
    end interface 
    
    private

    public solve_radial, compute_Vhartree, compute_radial_func


contains


    ! solve the radial equation for a give radial function and a given energy
    ! the radial equation is -2*(Vext + E - Veff)

    subroutine compute_radial_func(r,E,V,radial,Z)
        real(dp),intent(in) :: r(:),E,V(:)
        integer :: Z
        real(dp),intent(out) :: radial(:)

        radial(2:) = -2*(E + Z/r(2:)  - V(2:))


    end subroutine compute_radial_func


    subroutine  solve_radial(radial_func, u, r,params)
        implicit none
    
        real(dp), intent(out)    :: u(:)
        real(dp), intent(in)       :: r(:), params(:)
        real(dp),dimension(size(r))::w , radial_func
        integer  :: n , i
        real(dp) :: h
        real(dp) :: rmax 
        real(dp) :: norm


        


        ! print *,'solving for E = ' , params(1)

        n= size(u)
        h = r(2) - r(1)
        rmax = r(n)
        ! solve using numerov stepper
        u(n-1:n) =  r(n-1:n) * exp(-r(n-1:n))
        
        w(n) = (1- h*h/12_dp * radial_func(n))* u(n)
        w(n-1) = (1- h*h/12_dp * radial_func(n-1))* u(n-1)
        do i = n-2,2,-1
              w(i) = 2*w(i+1) - w(i+2) + h*h*radial_func(i + 1)*u(i+1)
              u(i) = w(i)/(1 - h*h/12_dp * radial_func(i))
        enddo
        
        
        u(1) = 2 * u(2) - u(3) +  h*h*radial_func(2)*u(2)
        norm = trapz(u**2,h)*4* pi
        ! if (norm  < 1e-15) stop "norm equal to zero for wavefunction for energy "//params(1)
        u = u / sqrt(norm)
    end subroutine solve_radial

    ! This routine solve the Poison equation for the give electronic density rho 
    ! and returns the value of the Hartee potentia in Vhartree
    ! h is the step in the discretization of the space
    subroutine compute_Vhartree ( rho, Vhartree, h)
        implicit none
        real(dp), intent(in)  :: rho(:)
        real(dp), intent(out) :: Vhartree(:)
        real(dp), intent(in)  :: h
        integer               :: nsteps, i
        real(dp),dimension(size(rho)) :: U,r
        real(dp) :: qmax, alpha
        ! solve here the problem for U using the density
        nsteps = size(Vhartree)
        r= (/(i*h,i=0,nsteps-1)/)

        ! solve the auxiliar problem for U(r) = r Vh(r)
        U(1) = 0.0_dp
        U(2) = h

        do i = 3, nsteps
            U(i)  =  2  * U(i-1) - U(i-2) + h*h*(-4_dp * pi*r(i)*rho(i)) 
        end do

        ! check this later on can be trouble
        qmax = trapz(r(2:nsteps)**2*rho(2:nsteps),h) * 4 * pi;

        alpha = (qmax - U(nsteps))/r(nsteps)
        U = U + r*alpha;
        Vhartree = U/r
        Vhartree(1) = 0.0_dp

    end subroutine compute_Vhartree

end module dft
