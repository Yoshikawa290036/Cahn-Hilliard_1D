! cal phi by Cahn-Hilliard equation

subroutine calphi(ni, u, dxinv, phi, dt, a, b, temperature, kappa)
    implicit none
    integer :: ni
    double precision :: u
    double precision :: phi(-3:ni+4), nphi(-3:ni+4)
    double precision :: dt
    double precision :: a, b, temperature, kappa
    double precision :: dxinv
    integer :: i
    double precision :: adv, dif, J(-1:ni+2), zeta

    do i = -1, ni+2
        zeta = temperature/((1.0d0-b*phi(i))**2)-2.0d0*a*phi(i)
        J(i) = -zeta*(0.5d0*(-phi(i-1)+phi(i+1))*dxinv) &
          & +kappa*phi(i)*(0.5d0*dxinv**3*(-phi(i-2)+2.0d0*phi(i-1)-2.0d0*phi(i+1)+phi(i+2)))
    end do

    do i = 0, ni+1
        adv = -0.5d0*dxinv*(-u*phi(i-1)+u*phi(i+1))
        dif = -0.5d0*dxinv*(-12.0d0*phi(i-1)*J(i-1)+12.0d0*phi(i+1)*J(i+1))
        nphi(i) = phi(i)+dt*(adv+dif)
    end do

    do i = 0, ni+1
        phi(i) = nphi(i)
    end do

end subroutine calphi
