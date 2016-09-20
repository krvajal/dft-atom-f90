module io_utils
    use double
implicit none
private 
public print_vector
contains
    ! print a vector in column form 
    ! to the console
    subroutine print_vector(v)
          implicit none
          real(dp) :: v(:)
          integer :: i, n
          n = size(v)
          do i=1,n
            print *,v(i)
          enddo
            
      end subroutine print_vector  

end module io_utils