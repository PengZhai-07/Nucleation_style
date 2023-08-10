program addmats 
    implicit none

    integer, parameter :: dimmat = 3
    real, dimension(dimmat,dimmat) :: a,b,c
    integer :: i,j

    ! This creates the matrices.
    a(1,2) = 2.0
    do i =2, dimmat-1
    a(i,i+1) = 2.0
    b(i,i-1) = 1.0
    enddo
    b(dimmat, dimmat-1) = 1.0

    ! This adds the matrices a and b
    ! do i = 1, dimmat
    ! do j = 1, dimmat
    ! c(i,j) = a(i,j) + b(i,j)
    ! enddo
    ! enddo
    c = a+b
    c(1,1)=0
    c(1,3)=0
    c(2,2)=0
    c(3,1)=0
    c(3,3)=0
    ! This prints c
    do i = 1, dimmat
        write(*,*) ( c(i,j), j=1, dimmat )
    enddo

end program addmats