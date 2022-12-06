program maths
    implicit none
    real :: r, area
    real :: r_rd, area_rd
    real, parameter :: pi = 4.*atan(1.)
    character(len=20) :: filename
    r = 2.0
    area = pi*r**2
    filename='results.out'


    ! write to ascii file
    open(unit=1,file=filename,status='replace',form='formatted') !rewrite the file if exists
    write(1,11) r,area
    11 format(F3.1,F5.2)
    close(1)
    
    ! read from ascii file
    open(unit=2,file=filename,status='old') !read from the existing file
    read(2,11) r_rd, area_rd
    close(2)
    
    write(*,12) r_rd, area_rd
    12 format(F3.1,2X,F5.2)
end program maths