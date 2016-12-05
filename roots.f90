module roots
    use double
implicit none
private
 public bracket_zeros, bisec
contains

  subroutine bracket_zeros(v, num,intervals,found)
    integer,optional :: num
    real(dp) :: v (:)
    integer :: intervals(:)
    real(dp) :: ff 
    integer :: i,n
    integer,intent(out) :: found 
    found  = 0  
    n = size(v)
    ff = v(1)
    i  = 2
    if (.not. present(num)) num = 1
    
    do while(found < num .and. i < n)
         if((v(i)*ff) < 0 ) then
            intervals(found + 1) = i-1 ! start of the interval where a zero was found
            found = found + 1
            ff = v(i)
         endif  
         i = i + 1
    enddo

  end subroutine bracket_zeros
  
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

end module roots