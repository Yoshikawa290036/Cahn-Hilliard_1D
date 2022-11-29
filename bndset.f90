
! periodic boundary conditions

subroutine bundset(ni, phi)
    implicit none

    integer :: ni
    double precision :: phi(-2:ni+3)

    phi(0) = phi(ni)
    phi(-1) = phi(ni-1)
    phi(-2) = phi(ni-2)

    phi(ni+1) = phi(1)
    phi(ni+2) = phi(2)
    phi(ni+3) = phi(3)
end subroutine bundset
