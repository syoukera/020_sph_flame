subroutine enth_trans
!
    use main_variables
    real*8 a_i(nmax), b_i(nmax), c_i(nmax), d_i(nmax)
    real*8 phi(nmax), S_i(nmax), G_i(nmax)
    integer n_up
    real*8, parameter :: x_L = 0.0257      !!temporary heat trans coef (N2@300K[W/mK]) 
    real*8, parameter :: C_p = 1006.0      !!temporary specific heat   (N2@300K[J/kgK]) 
!   --------- n_up: upwind switch u_up=1 ----------------
    n_up = 2
!       
    do n=1, nmax
        phi(n) = temp(n)
        S_i(n) = 0.0d0
        G_i(n) = x_L/C_p
    end do
!
    call calc_scl_coef(phi, G_i, S_i, a_i, b_i, c_i, d_i,n_up)
!
    call calc_tdma(phi, a_i, b_i, c_i, d_i)
!
    do n=1, nmax
        temp(n) = phi(n)
    end do
!
end subroutine