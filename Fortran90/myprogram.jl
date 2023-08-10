using FortranFiles
using OffsetArrays
using Parameters
using Printf

function myprogram()
#implicit none

#integer :: i, j, k
#REAL :: x, y, z

x = 3.61
y = cos(x)
z = x + y

i = 3
j = i^2
k = i - j

println(stdout, i, j, k, x, y, z)

f = open("mydata.dat", "r")
println(stdout, "The value of x is ", x, " and the value of y is", i)
close(f)

# write(*,*) 'What is the value of x?'
# read(*,*) x
# write(*,*) 'x is equal to  ', x

OPEN(11, file = "input.dat")
READ(11, x)
READ(11, i)
println(stdout, "x is equal to  ", x)
println(stdout, "i is equal to  ", i)

CLOSE(11)

end