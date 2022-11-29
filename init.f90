subroutine init(ni, phi, phimin, phimax, width)
    implicit none

    integer :: ni, width
    double precision :: phimin, phimax
    double precision :: phi(-2:ni+3)

    integer :: i

    do i = -2, ni+3
        if (i < ni/2-width/2 .or. i > ni/2+width/2) then
            phi(i) = phimax
        else
            phi(i) = phimin
        end if
    end do
end subroutine init
