module help_fortranio

! Functions used by test_fortranio.cpp
! This uses real Fortran to write out a file, so it can be read in C++

use iso_c_binding

implicit none

contains

subroutine write_single_sample(fname,convert)
    character(*), intent(in) :: fname, convert
    ! -------- Local Vars
    character(len=80) :: str0, str1
    integer :: vals(17)
    integer :: i,j
    integer :: data(4,3)


    do i=1,17
        vals(i) = i
    end do
    str0 = 'Hello World'
    str1 = 'Goodbye World'

    open (1, File=fname, Form='Unformatted', Status='REPLACE', CONVERT=convert)
    write (1) str0, vals
    write (1) vals, str0

    do i=1,17
        vals(i) = 17-i
    end do
    write (1) str0,vals,str1

    ! [4] Write a 2D array
    do j=1,3
    do i=1,4
        data(i,j) = (i-1) * j;
    end do
    end do
    write (1) str0, data, str1

    close (1)

end subroutine write_single_sample

subroutine write_sample() bind(C)
    call write_single_sample('sample_fortranio_be', 'BIG_ENDIAN')
    call write_single_sample('sample_fortranio_le', 'LITTLE_ENDIAN')
end subroutine write_sample


end module help_fortranio
