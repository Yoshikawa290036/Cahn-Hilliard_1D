
! periodic boundary conditions

subroutine bundset(ni, phi)
    implicit none

    integer :: ni
    double precision :: phi(-3:ni+4)

    phi(0) = phi(ni)
    phi(-1) = phi(ni-1)
    phi(-2) = phi(ni-2)
    phi(-3) = phi(ni-3)

    phi(ni+1) = phi(1)
    phi(ni+2) = phi(2)
    phi(ni+3) = phi(3)
    phi(ni+3) = phi(4)
end subroutine bundset
