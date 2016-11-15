module io_utils
    use double
implicit none
private 
public print_vector, write_to_file
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

      subroutine write_to_file(filename,x,y)
        implicit none
        character(len=*) ::filename
        real(dp),intent(in) ::  x(:),y(:)
        integer :: u, i

        open(newunit = u, file=filename,status="replace")
        
        do i=1,size(x)
          write (u,*) x(i),y(i)
        enddo
        close(u)
      end subroutine
      integer function newunit(unit) result(n)
        ! returns lowest i/o unit number not in use
        integer, intent(out), optional :: unit
        logical inuse
        integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
        integer, parameter :: nmax=999  ! may be system-dependent
        do n = nmin, nmax
            inquire(unit=n, opened=inuse)
            if (.not. inuse) then
                if (present(unit)) unit=n
                return
            end if
          end do
          stop "newunit ERROR: available unit not found."
        end function
end module io_utils