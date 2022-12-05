! cal phi by Cahn-Hilliard equation

subroutine calphi(ni, u, dxinv, phi, dt, a, b, temperature, kappa)
    implicit none
    integer :: ni
    double precision :: u
    double precision :: phi(-3:ni+4), nphi1(-3:ni+4), nphi2(-3:ni+4)
    double precision :: dt
    double precision :: a, b, temperature, kappa
    double precision :: dxinv
    integer :: i
    double precision :: adv1(-1:ni+1), dif1(-1:ni+1), adv2(-1:ni+1), dif2(-1:ni+1)
    double precision :: J(-1:ni+2), zeta
    double precision :: Gamma

    Gamma = 64.0d0

    do i = -1, ni+2
        zeta = (temperature/((1.0d0-b*phi(i))**2))-2.0d0*a*phi(i)
        J(i) = -zeta*(0.5d0*(-phi(i-1)+phi(i+1))*dxinv) &
             & +kappa*phi(i)*(0.5d0*dxinv**3*(-phi(i-2)+2.0d0*phi(i-1)-2.0d0*phi(i+1)+phi(i+2)))
    end do

    do i = 0, ni+1
        adv1(i) = 0.5d0*dxinv*(-u*phi(i-1)+u*phi(i+1))
        dif1(i) = 0.5d0*dxinv*(-Gamma*phi(i-1)*J(i-1)+Gamma*phi(i+1)*J(i+1))
        nphi1(i) = phi(i)-dt*(adv1(i)+dif1(i))
    end do

    do i = -1, ni+2
        zeta = (temperature/((1.0d0-b*nphi1(i))**2))-2.0d0*a*nphi1(i)
        J(i) = -zeta*(0.5d0*(-nphi1(i-1)+nphi1(i+1))*dxinv) &
             & +kappa*nphi1(i)*(0.5d0*dxinv**3*(-nphi1(i-2)+2.0d0*nphi1(i-1)-2.0d0*nphi1(i+1)+nphi1(i+2)))
    end do

    do i = 0, ni+1
        adv2(i) = 0.5d0*dxinv*(-u*nphi1(i-1)+u*nphi1(i+1))
        dif2(i) = 0.5d0*dxinv*(-Gamma*nphi1(i-1)*J(i-1)+Gamma*nphi1(i+1)*J(i+1))
        nphi2(i) = phi(i)-dt*0.5d0*(adv1(i)+dif1(i)+adv2(i)+dif2(i))
    end do

    do i = 0, ni+1
        phi(i) = nphi2(i)
    end do

end subroutine calphi
