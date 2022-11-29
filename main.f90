program main
    implicit none

    integer :: ni
    integer :: i
    double precision :: x
    double precision :: phimin, phimax
    double precision :: dx
    double precision :: xl
    double precision :: t, dt
    double precision :: u
    double precision, dimension(:), allocatable :: phi
    integer :: maxstep, step
    character(32) fname

    maxstep = 10
    ni = 128
    dx = 1.0d0
    phimin = 0.265
    phimax = 0.405
    xl = dx*dble(ni)
    dt = 2.5e-2
    u = 0.5d0
    step = 0
    include'allocate.h'

    call init(ni, phi, phimin, phimax, 32)

    include'mkphi.h'

    do step = 1, maxstep
        call bundset(ni, phi)

    end do

end program main
