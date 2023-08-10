program myprogram
implicit none

integer :: i, j, k
real :: x, y, z

x = 3.61
y = cos(x)
z = x + y

i = 3
j = i**2
k = i - j

write(*,*) i,j,k,x,y,z

open(10, file='mydata.dat')
write(10,*) 'The value of x is ',x,' and the value of y is', y
close(10) 

! write(*,*) 'What is the value of x?'
! read(*,*) x
! write(*,*) 'x is equal to  ', x

open(11, file='input.dat')
read(11,*) x
read(11,*) i
write(*,*) 'x is equal to  ', x
write(*,*) 'i is equal to  ', i

close(11)

end program myprogram