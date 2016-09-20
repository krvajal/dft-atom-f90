
program dft
! --------------------------------------------------------------
! Program to solve the radial part of Schroedinger equation 
! for the hidrogen atom
! --------------------------------------------------------------
use quad, only: trapz
use double
use ode, only: verlet
use dft
use io_utils
use roots

implicit none


    
    real(dp)                :: r_step = 1e-2,params(2), currentr
    real(dp)                    :: aux
    real(dp),parameter      :: rmax = 30
    integer,parameter               :: steps=1000
    integer :: i
    real(dp),allocatable    :: r(:), u(:)
    real(dp)                :: stiffness =  1.0_dp
    real(dp)                :: x0, x1, x2, t
    integer                 :: n
    real(dp)                :: E(steps),u0(steps)   
    real(dp),parameter                :: E_start= -0.6, E_end=-0.01, E_step = (E_end - E_start)/(steps-1)
    real(dp) :: E0 
    real(dp),parameter,dimension(1) :: EMPTY_PARAMS = (/0.0_dp/)
    integer :: intervals(3),interv_found



    n = rmax/r_step + 1
    allocate(r(n),u(n))
    r = (/(i*r_step,i=0,n-1)/)



    E = (/(E_start + i*E_step,i=0,steps-1)/)
    call energies(E,u0,E_step)

    call bracket_zeros(u0,3, intervals,interv_found)
    do i = 1, interv_found
        E0 = bisec(u0f,EMPTY_PARAMS, E(intervals(i)),E(intervals(i)+1),0.001_dp)    
        print *, E0
    end do
    
    
    ! params = -0.125
    ! call solve_radial(radial,u,r,params)
!   
    ! Find the numerical root in the first interval
    !call energies(-0.6_dp, -0.1_dp, 0.001_dp)
!   E =  bisec(solve_radial, params, -0.2d+0, -0.1d+0, 1d-6)

!   !get the values of the points in that interval
!   currentr = rmax
!   do i=1,steps
!       u(i) = solve_radial(E, currentr)
!       r(i) = currentr
!       currentr = currentr - h
        
! !         if ( currentr < h) exit 
!   end do
!   aux = trapz(u**2, h)
!   u = u/( sqrt(aux * 4 * 3.1415))
!   do i =1,steps
!       print *,r(i),u(i)
!   end do
        
contains 

    real(dp) function osc_ham(t, x, params)
        real(dp), intent(in) :: x,t, params(:)
        ! here the t is  not used but it can be used
        ! for other second order functions
        osc_ham  = - params(1)* x1

    end function osc_ham



    real(dp) function u2_r(t, x, params)
        real(dp), intent(in) :: x,t, params(:)
        ! here the t is  not used but it can be used
        ! for other second order functions
        u2_r  = - 4.0_dp * t * exp( - 2 * t)

    end function u2_r


    ! solves the radial Schroedinger equation 
    ! using numerov method and returns the value
    ! of u(0,E) for the given enery E



    real(dp) function radial(r, params)
        implicit none
        real(dp), intent(in) :: r, params(:)
        real(dp)                :: E

        E =  params(1)

        radial  = -2 * (E + 1/r)
    end function radial

    real(dp) function u0f(E,params)
        real(dp),intent(in) :: E,params(:)
        real(dp) :: p(1)
        ! print *, "called with E = ", E
        p = E    
        call solve_radial(radial,u,r,p)
        u0f = u(1)
        ! print *, "u0= ", u0f
    end function u0f

    real(dp) function bisec(func, params, x1, x2, xacc)
        implicit none
        real(dp), intent(in)  :: x1, x2, xacc, params(:)
        integer, parameter   :: jmax = 40
        integer              :: j
        real(dp)              :: dx, f, fmid, xmid
        interface   
        real(dp) function ff(x, params)
                use double
                real(dp), intent(in) :: x, params(:)
            end function ff
        end interface
        procedure(ff) :: func


        fmid = func(x2,params)
        f = func(x1,params)
        if ((fmid * f) > 0)  stop 'root must be bracketed in bisec'
        if (f < 0) then
            bisec = x1
            dx = x2 - x1
        else 
            bisec = x2
            dx = x1 - x2
        end if
        do j = 1, jmax
            dx = dx * 0.5
            xmid  = bisec + dx
            fmid = func(xmid, params)
            if (fmid < 0) bisec = xmid
            if (abs(dx) < xacc .or. fmid == 0) return

        end do

    end function bisec

    ! compute the value of u(0,E) in the range 
    ! from E1 to E2 in steps of h
    subroutine energies(E,u0,h)
        implicit none
        real(dp), intent(in) :: E(:),h
        real(dp), intent(out) :: u0(:)
        integer :: N = 1001 
        integer      :: steps, i

        steps = size(E)
        
        do i = 1, steps
            params = E(i) ! the value of E is pased as a param to radial
        call solve_radial(radial,u,r,params)
        u0(i) = u(1)
        end do
    end subroutine energies




end program  dft


