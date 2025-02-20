! Exersice 2
! want to compute the sum

program variables
    implicit none

    integer*2 :: int1, int2
    integer*4 :: int3, int4
    real*4 :: pi1 , real1 , real2
    real*8 :: pi2 , real3, real4

    int1 = 2000000
    int2 = 1
    print *, "With INTEGER*2:", int1+int2
    ! we want to sum in two different ways
    int3 = 2000000
    int4 = 1
    print *, "With INTEGER*4:", int3+int4

    pi1 = acos(-1.)
    real1 = pi1*10e32
    real2 = sqrt(2.)*10e21
    print *, "With single precision:", real1+real2

    pi2 = 4.D0*datan(1.D0)
    real3 = pi2*10e32
    real4 = sqrt(2.)*10e21
    print *, "With double precision:", real3+real4

    stop

end program variables
