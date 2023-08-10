program solbybisection
    implicit none
    real, parameter :: xmin=0.0, xmax=2.0
    integer, parameter :: iter = 10
    real :: sol, err, fcosx
    external fcosx

    call bisect(xmin, xmax, fcosx, sol, iter, err)

    write(*,*) sol, err

end program solbybisection