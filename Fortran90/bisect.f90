subroutine bisect(xmin,xmax,func,sol,iter,err)
    implicit none

    real :: xmin,xmax,sol,func
    integer :: iter
    real:: x1,x2,res,err
    integer :: i
    external func

    x1 = xmin
    x2 = xmax
    res = func(x1)*func(x2)
    if(res.gt.0.0) then
    write(*,*)    'Try again with a different interval.'
    endif

    do i=1,iter
    sol = (x1+x2)/2
    res = func(sol)*func(x2)
    if(res.gt.0.0) then
    x2 = sol
    else
    x1 = sol
    endif
    enddo

    sol = (x1+x2)/2
    err = (x2-x1)/2

end subroutine bisect

