module solvers
! This module contains the subroutines for solving
! the Schrödinger radial equation and 
! the Poisson equation using a linear mesh
use types
use constants
use quad
implicit none



contains

    ! =============================================
    ! This routine solve the Poisson equation for the
    ! give electronic density rho 
    ! and returns the value of the Coulomb potential 
    ! h is the step in the discretization of the space
    ! The verlet algorith is used with a linear mesh
    ! =============================================

    subroutine SolvePoisson (rpoints, rho, Potential, h)
        implicit none
        real(dp), intent(in)  :: rho(:),rpoints(:)
        real(dp), intent(out) :: Potential(:)
        real(dp), intent(in)  :: h
        integer               :: nsteps, i
        real(dp),dimension(size(rho)) :: U,r
        real(dp) :: qmax, alpha
        ! solve here the problem for U using the density
        nsteps = size(rho)


        r= (/(i*h,i=0,nsteps-1)/)
        

        ! solve the auxiliar problem for U(r) = r V(r)
        U(1) = 0.0_dp
        U(2) = h
        
        do i = 3, nsteps
            U(i)  =  2  * U(i-1) - U(i-2) + h*h*(-4_dp * pi*r(i-1)*rho(i-1)) 
        end do

        qmax = trapz(rho * rpoints**2,h) * 4 * pi;
        print *, "qmax", qmax   
        print *, "U(rmax)",U(nsteps)    
        alpha = (qmax - U(nsteps))/r(nsteps)
        print *, "alpha",alpha
        U = U + r*alpha;

        Potential(2:) = U(2:)/r(2:)
        Potential(1) = 0.0_dp

    end subroutine SolvePoisson

    ! =============================================
    ! Solve the Schodinger radial equation for a given
    ! function f  = -2(E  - V) in a lienear mesh
    ! Numerov method is used
    ! =============================================
    subroutine  SolveSchodingerRadial(RadialFunc, u, r)
        implicit none
    
        real(dp), intent(out)    :: u(:)
        real(dp), intent(in)       :: r(:)
        real(dp),dimension(size(r))::w , RadialFunc
        integer  :: n , i
        real(dp) :: h
        real(dp) :: rmax 
        real(dp) :: norm


        n= size(u)
        h = r(2) - r(1)
        rmax = r(n)
        ! solve using numerov stepper
        u(n-1:n) =  r(n-1:n) * exp(-r(n-1:n))
        
        w(n) = (1- h*h/12_dp * RadialFunc(n))* u(n)
        w(n-1) = (1- h*h/12_dp * RadialFunc(n-1))* u(n-1)
        
        do i = n-2,2,-1
              w(i) = 2*w(i+1) - w(i+2) + h*h*RadialFunc(i + 1)*u(i+1)
              ! print *, RadialFunc(i)
              u(i) = w(i)/(1 - h*h/12_dp * RadialFunc(i))
        enddo
        
        
        
        u(1) = 2 * u(2) - u(3) +  h*h*RadialFunc(2)*u(2)
        norm = trapz(u**2,h)*4* pi
        if (abs(norm )  < tiny(1.0_dp)) stop "norm equal to zero for wavefunction for energy "

        u = u / sqrt(norm)

    end subroutine SolveSchodingerRadial


    ! =============================================
    ! Returns the value of the radial 
    ! part of the Schrodinger equation for a give 
    ! energy
    ! ============================================
    real(dp) function UStartFunctor(E,params)

        real(dp),intent(in) :: E,params(:)
        integer :: n
        real(dp),dimension(size(params)/2) :: radial, r, V,u
        ! print *, "called with E = ", E
        n = size(params)/2
    
        r = params(1:n)
        V= params(n+1:)
        radial = -2.0_dp *(E - V )
        call SolveSchodingerRadial(radial,u,r)
        UStartFunctor = u(1)
        ! print *, "u0= ", u0f
    end function UStartFunctor

    ! =============================================
    ! Construct the radial fuction 
    ! Parameters:
    !       E: the energy
    !       V: the effective potential
    ! Output 
    !       radial:  -2*(E-Veff)
    ! =============================================
    subroutine ComputeRadialFunct(E,Veff,radial)
        real(dp),intent(in) :: E,Veff(:)
        real(dp),intent(out) :: radial(:)
        radial  = -2*(E - Veff)
    end subroutine ComputeRadialFunct

end module solvers
