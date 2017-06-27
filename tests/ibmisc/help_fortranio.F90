module help_fortranio

! Functions used by test_fortranio.cpp
! This uses real Fortran to write out a file, so it can be read in C++

use iso_c_binding

implicit none

contains

subroutine write_sample() bind(C)
    character(len=80) :: str0, str1
    integer :: vals(17)
    integer :: i

    do i=1,17
        vals(i) = i
    end do
    str0 = 'Hello World'
    str1 = 'Goodbye World'

    open (1, File='sample_fortranio', Form='Unformatted', Status='REPLACE')
    write (1) str0, vals
    write (1) vals, str0

    do i=1,17
        vals(i) = 17-i
    end do
    write (1) str0,vals,str1

    close (1)

end subroutine write_sample

end module help_fortranio
