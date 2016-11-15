
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
    real(dp) :: param(1)
    character(1024) :: filename
    integer :: intervals(3),interv_found



    call prob4()
    stop


    ! n = rmax/r_step + 1
    ! allocate(r(n),u(n))
    ! r = (/(i*r_step,i=0,n-1)/)



    ! E = (/(E_start + i*E_step,i=0,steps-1)/)
    ! call energies(E,u0,E_step)
    ! ! do i =1,steps
    ! !     print *,E(i), u0(i)
    ! ! enddo
    
    ! call bracket_zeros(u0,3, intervals,interv_found)
    ! do i = 1, interv_found
    !     E0 = bisec(u0f,EMPTY_PARAMS, E(intervals(i)),E(intervals(i)+1),0.001_dp) 
    !     param = E0
    !     call solve_radial(radial, u,r,param)
    !     print *, "E0 =", E0
    !     write(filename,"(A1,I1,A4)") "U", i,".dat"
    !     print *,trim(filename)
    !     call write_to_file(filename,r,u)
    ! end do
    
    
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

    subroutine prob2()
        integer,parameter :: n = 5001 
        real(dp) :: u(n),r(n), rho(n)
        real(dp) :: rmax = 30
        real(dp) :: Vh(n)
        integer :: i
        real(dp) :: h
        character(len=20) :: filename = "Vh_u1.dat"

        h =  rmax/(n-1)
        !solve the poisson equation for a given electron density
        r = (/(i*h,i=0,n-1)/)
        u = r * exp(-r)/sqrt(pi)

        rho = (u/r )**2
        call compute_Vhartree(rho,Vh,h)
        call write_to_file(filename, r, Vh)


    end subroutine prob2

    subroutine prob2b()

        integer,parameter :: n = 5001 
        real(dp) :: u(n),r(n), rho(n)
        real(dp) :: rmax = 30
        real(dp) :: Vh(n), radial(n)
        integer :: i,Z = 2
        real(dp) :: h, energy0_old,energy0, ehartree, etotal
        integer,parameter :: steps = 1000
        real(dp) :: E(steps), u0(steps),param(1)
        integer :: interv_found
        integer :: intervals(2), iter = 1
        real(dp) :: tol = 1e-6

        real(dp),parameter :: E_start= -5, E_end=-0.01,&
                              E_step = (E_end - E_start)/(steps-1)
        Vh = 0
        h =  rmax/(n-1)
        r = (/(i*h,i=0,n-1)/)

        E = (/ (E_start + i * E_step, i=0,steps-1)/)
        do 
            print *, "iter =",iter
            call energies(E,u0,E_step,Vh,r)

            ! search for the lowest energy
            call bracket_zeros(u0,1, intervals,interv_found)
        
            if (interv_found == 0) stop 'no roots founds'
            energy0_old = energy0
            energy0 = bisec(u0f,(/r,Vh /), E(intervals(1)),E(intervals(1)+1),0.001_dp) 
            
            param = energy0
            print *, "e0 =", energy0
            
            call compute_radial_func(r,energy0,Vh,radial,Z)
            call solve_radial(radial, u,r,param)
            rho(2:) = (u(2:)/r(2:) )**2
            !solve the poisson equation for a given electron density
            call compute_Vhartree(rho,Vh,h)
            ehartree = trapz(Vh*u**2,h) * 4 * pi
            
            print *,"ehartree =", ehartree

            etotal = 2* energy0 - ehartree
            print *,"etotal =", etotal
            iter = iter +  1
            if( abs(energy0 - energy0_old) < tol) exit

        enddo
        call write_to_file("prob2b_u_r.dat",r,u)



    end subroutine prob2b


    subroutine prob3()


        integer,parameter :: n = 10001
        real(dp) :: u(n),r(n), rho(n)
        real(dp) :: rmax = 10
        real(dp) ::  radial(n), Vx (n),V(n),Vh(n)
        integer :: i,Z = 2
        real(dp) :: h, energy0_old,energy0, ehartree,eexchange, etotal, etotal_old
        integer,parameter :: steps = 1000
        real(dp) :: E(steps), u0(steps),param(1)
        integer :: interv_found
        integer :: intervals(2), iter = 1
        real(dp) :: tol = 1e-8

        real(dp),parameter :: E_start= -5, E_end=-0.01,&
                              E_step = (E_end - E_start)/(steps-1)
        Vh = 0;
        Vx= 0
        h =  rmax/(n-1)
        print *,h
        
        r = (/(i*h,i=0,n-1)/)

        E = (/ (E_start + i * E_step, i=0,steps-1)/)
        do 
            print *, "iter =",iter
            V = Vx + Vh
            call energies(E,u0,E_step,V,r)

            ! search for the lowest energy
            call bracket_zeros(u0,1, intervals,interv_found)
            call write_to_file("tmp.dat", E,u0)

            if (interv_found == 0) stop 'no roots founds'
            energy0_old = energy0
            energy0 = bisec(u0f,(/r, V /), E(intervals(1)),E(intervals(1)+1),0.001_dp) 
            
            param = energy0
            print *, "e0 =", energy0
            
            call compute_radial_func(r,energy0,V,radial,Z)
            call solve_radial(radial, u,r,param)
            rho(2:) = (u(2:)/r(2:) )**2
            rho(1) = 0
            !solve the poisson equation for a given electron density
            call compute_Vhartree(2*rho,Vh,h)
            call compute_Vx(r,2*rho,Vx) !the total density is twice the density 
            ehartree = trapz(Vh*u**2,h) * 4 * pi
            eexchange = trapz(Vx*u**2,h) *4 * pi
            print *,"ehartree =", ehartree
            print *,"eexchange =", eexchange
            

            etotal_old = etotal
            etotal = 2* energy0 - ehartree  - 0.5*eexchange
            print *,"etotal =", etotal
            iter = iter +  1
            if( abs(energy0 - energy0_old) < tol .and. abs(etotal - etotal_old) < tol) exit

        enddo
        call write_to_file("prob3_u_r.dat",r,u)


    end subroutine prob3

    subroutine prob4()
         integer,parameter :: n = 2001 
        real(dp) :: u(n),r(n), rho(n)
        real(dp) :: rmax = 10
        real(dp) ::  radial(n), Vx (n),V(n),Vh(n),Vc(n),ec(n)
        integer :: i,Z = 2
        real(dp) :: h, energy0_old,energy0, ehartree,eexchange, etotal, etotal_old, ecorrelacion
        integer,parameter :: steps = 1000
        real(dp) :: E(steps), u0(steps),param(1)
        integer :: interv_found
        integer :: intervals(2), iter = 1
        real(dp) :: tol = 1e-5
        integer :: PotFileHandler



        real(dp),parameter :: E_start= -5, E_end=-0.01,&
                              E_step = (E_end - E_start)/(steps-1)


        open(newunit=PotFileHandler, file="potentials.dat",status="replace")

        Vh = 0; Vx= 0; Vc = 0
        ec = 0
        h =  rmax/(n-1)
        r = (/(i*h,i=0,n-1)/)

        E = (/ (E_start + i * E_step, i=0,steps-1)/)
        do 
            print *, "iter =",iter
            V = Vx + Vh + Vc
            call energies(E,u0,E_step,V,r)

            ! search for the lowest energy
            call bracket_zeros(u0,1, intervals,interv_found)
            call write_to_file("tmp.dat", E,u0)

            if (interv_found == 0) stop 'no roots founds'
            energy0_old = energy0
            energy0 = bisec(u0f,(/r, V /), E(intervals(1)),E(intervals(1)+1),0.000001_dp) 
            
            param = energy0
            print *, "e0 =", energy0
            
            call compute_radial_func(r,energy0,V,radial,Z)
            call solve_radial(radial, u,r,param)
            rho(2:) = (u(2:)/r(2:) )**2
            rho(1) = 0
            !solve the poisson equation for a given electron density
            call compute_Vhartree(2*rho,Vh,h)
            call compute_Vx(r,2*rho,Vx)
            call compute_Vc_and_ec(r(2:), 2*rho(2:),Vc(2:),ec(2:))



            ehartree = trapz(Vh* u**2,h) * 4 * pi
            eexchange = trapz(Vx*u**2,h) *4 * pi
            ecorrelacion =  (trapz(ec*u**2,h)- trapz(Vc *u **2 ,h))*2*4 *pi
            print *,"ehartree =", ehartree
            print *,"eexchange =", eexchange
            print *, "ecorrelacion =", ecorrelacion

            etotal_old = etotal
            etotal = 2* energy0 - ehartree -  0.5 * eexchange - ecorrelacion
            print *,"etotal =", etotal
            iter = iter +  1
            if( abs(energy0 - energy0_old) < tol .and. abs(etotal - etotal_old) < tol) exit

        enddo
        print *, "Convergence achieved"
        call write_to_file("prob4_u_r.dat",r,u)
        do i = 1,n
            write(PotFileHandler,*) r(i), Vh(i), Vx(i), Vc(i)
        enddo
        close(PotFileHandler)



    end subroutine prob4


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
        integer :: i
        real(dp),dimension(size(params)-1) :: Vh
        real(dp) :: h = 0.01
        E =  params(1)
        Vh = params(2:size(params))
        i = r/h + 1
        radial  = -2 * (E + 1/r -Vh(i))
    end function radial

    real(dp) function u0f(E,params)

        real(dp),intent(in) :: E,params(:)
        real(dp) :: p(1)
        integer :: n
        real(dp),dimension(size(params)/2) :: radial, r, V,u
        ! print *, "called with E = ", E
        n = size(params)/2
        r = params(1:n)
        V= params(n+1:)
        call compute_radial_func(r,E,V,radial,2)    
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

    subroutine compute_Vx(r,rho,Vx)
        real(dp),intent(in) :: r(:), rho(:)
        real(dp),intent(out):: Vx(:)
        
        
       
        Vx =  - (3.0_dp/pi*rho)**(1.0/3)

    end subroutine compute_Vx


    ! compute the value of u(0,E) in the range 
    ! from E1 to E2 in steps of h
    subroutine energies(E,u0,h,V,r)
        implicit none
        real(dp), intent(in) :: E(:),h,V(:),r(:)
        real(dp), intent(out) :: u0(:)
        real(dp),dimension(size(r)) :: u
        integer      :: steps, i
        real(dp),dimension(size(r)) :: radial

        steps = size(E)
        
        do i = 1, steps
            call compute_radial_func(r,E(i),V,radial,2)
            call solve_radial(radial,u,r,params)
            u0(i) = u(1)
        end do
    end subroutine energies

    subroutine compute_Vc_and_ec(r,rho,Vc,ec)
        implicit none
        real(dp),intent(in) ::r(:), rho(:)
        real(dp), intent(out) :: Vc(:),ec(:)
        real(dp),dimension(size(r)) :: rs
        real(dp),parameter :: gamma = -0.1423, beta1 = 1.0529,&
                              beta2 =  0.3334, A = 0.0311, &
                              B     =  -0.048,C = 0.002,&
                              D = -0.0116 
        integer :: loc,n


        rs=  (3./(4*pi*rho))**(1./3)
        loc = minloc(abs(rs - 1.0_dp),1)
        if (rs(loc) > 1) loc = loc -1
        ! print *, rs(loc  -1)
        ! print *, rs(loc)
        ! print *, rs(loc + 1)
   
        n = size(r)
        ! rs >= 1
        Vc(loc+1:) = gamma/(1 + beta1*sqrt(rs(loc+1:)) + beta2*rs(loc+1:))
        Vc(loc+1:) = Vc(loc+1:)* (1 + 7./6 *beta1 * sqrt(rs(loc+1:))+ 4./3 * beta2 * rs(loc+1:))
        Vc(loc+1:) = Vc(loc+1:) / (1 + beta1*sqrt(rs(loc+1:)) + beta2*rs(loc+1:))
        !rs < 1 
        Vc(:loc) = A *log(rs(:loc)) + B - A/3 + 2./3*C*rs(:loc) + (2*D-C)*rs(:loc)/3.0_dp 

        !rs >=1
        ec(loc+1:) = gamma/(1 + beta1*sqrt(rs(loc+1:)) + beta2*rs(loc+1:))
        !rs < 1
        ec(:loc) = A *log(rs(:loc)) + B + C*rs(:loc) * log(rs(:loc)) + D* rs(:loc)


    end subroutine compute_Vc_and_ec




end program  dft


