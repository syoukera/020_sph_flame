subroutine grid_set
!
    use main_variables
!
!   grid setting
    delta_x = 0.05e-3
    do n=1, nmax
        xscl(n) = float(n-1)*delta_x
        xvel(n) = (0.5d0 + float(n-1))*delta_x
    end do
end subroutine