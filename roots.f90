module roots
    use double
implicit none
private
 public bracket_zeros
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

end module roots