program geometry
    implicit none
    type point
        real :: x,y
    END type point

    type line
        real :: a,b,c
    end type line

    type(point)::p1,p2
    type(line)::p1_to_p2

    write(*,*) "Please type co-ordinates of first point"
    read(*,*) p1
    write(*,*) "Please type co-ordinates of second point"
    read(*,*) p2

    p1_to_p2%a = p2%y - p1%y
    p1_to_p2%b = p1%x - p2%x
    p1_to_p2%c = p1%y*p2%x - p2%y*p1%x

    write(*,*) "where a = ", p1_to_p2%a
    write(*,*) "        b=", p1_to_p2%b
    write(*,*) "         c", p1_to_p2%c

end program geometry
