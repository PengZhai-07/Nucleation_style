program geometric
    implicit none     
    integer :: iter
    real :: a0, r, res
    write(*, *) 'What is the vlalue of a0?'
    read(*, *) a0
    write(*, *) 'What is the vlalue of r?'
    read(*, *) r

    do iter=1, 10
    write(*, *)  iter, a0
    a0 = a0*r
    enddo
        
end program geometric
