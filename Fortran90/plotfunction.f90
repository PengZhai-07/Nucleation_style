program plotfunction
    implicit none
    integer :: i
    real::x,f
    real,parameter :: xmin=0., xmax=10.,a=-2.     

    open(10, file='myplot.dat')
    do i = 1,100
        x = xmin + xmax*(i-1.0)/(100.0-1.0)
        write(10,*) x, f(x,a)
    enddo
    close(10)
    
end program plotfunction


