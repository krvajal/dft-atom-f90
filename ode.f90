module ode

use double
implicit none

private
public verlet, numerov

interface
	!define the callback function f
	real(dp) function  func2(t,x, params)
		use double, only: dp
		implicit none
		real(dp), intent(in) :: t,x
		real(dp), intent(in) :: params(:)
	end function func2

end interface 



contains

	real(dp) function  verlet(f, x0, x1, t, h , params) result(x2) 
		real(dp), intent(in)  :: x0,x1,t,h
		real(dp), intent(in)  :: params(:)
		procedure(func2) f
			x2 = 2*x1 -  x0 + h * h * f(t, x1, params)
	end function verlet


	! -----------------------------------
	! advance one step the value of x 
	! using the numerov method 
	! the ecuation is x''(t) = f(t)x(t)
	! t is the actual time
	! x is the value at time t 
	! xout will containt the value at x(t+h)
    subroutine numerov( x, t, h, f, params, xout, w)
        use double
        real(dp), intent(in)  :: x, t, h, params(:)
        real(8), intent(inout) :: w(2)
        real(8), intent(out) :: xout

        interface 
            function f(x, params)
            use double
             real(dp),intent(in) :: x, params(:)
             real(dp) :: f
			end function f
        end interface


        real(dp)              :: w1, w2
        real(dp)              :: tmp

        tmp = 2 * w(2) - w(1) + h * h * f(t,params) * x
        w(1) = w(2)
        w(2) = tmp
        xout = tmp/ (1.0 - h * h * f(t + h, params)/12.0 )

	end subroutine numerov


end module ode