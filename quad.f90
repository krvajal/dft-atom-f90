module quad

use double, only: dp

implicit none
private
public simpson, trapz

interface
    !define the callback function f
    real(dp) function callback(x, data)
        use double, only: dp
        implicit none
        real(dp), intent(in) :: x
        real(dp), intent(in) :: data(:)
    end function callback

end interface 

contains
    real(dp) function simpson(f,a,b , data) result(s)

        real(dp), intent(in) :: a,b, data(:)
        procedure(callback)  :: f
        ! function simpson body begins
        s = (b-a)/ 6.* (f(a,data) + 4 * f(a + b,data)/2 + f(b,data))
    end function simpson


    real(dp) function trapz(f,h) result(s)
        real(dp), intent(in) :: f(:)
        real(dp),optional :: h
        
        integer                ::   n
        real(dp)               :: accum = 0
        !default h to one
        if (.not. present(h) ) h = 1.0_dp
        n = size(f)
        accum = sum(f(2:n-1))

        accum = accum + (f(1) + f(n)) *0.5_dp
        accum = accum * h
        s = accum
    end function trapz

end module quad