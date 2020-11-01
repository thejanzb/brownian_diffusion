module parameters
    implicit none

    !! Precision Variables !!
    integer, parameter ::  sRp = selected_real_kind(  6, 37  ) ! Single    Real precision
    integer, parameter ::  dRp = selected_real_kind( 15, 307 ) ! Double    Real precision
    integer, parameter ::  mRp = selected_real_kind( 18, 4931) ! Medium    Real precision
    integer, parameter ::  qRp = selected_real_kind( 33, 4931) ! Quadruple Real precision
    
    ! Choosing double precision for real here !
    integer, parameter ::  uRp = dRp                           ! User selected Real precision

    !! Constants !!
    real(kind = uRp), parameter                      :: pi = 4.0_uRp*atan(1.0_uRp)
    real(kind = uRp), parameter                      :: kb = 1.38e-23_uRp
    integer(kind = selected_int_kind(23)), parameter :: av_num = 6023*10**20

    !! Particle properties !!
    real(kind = uRp), parameter :: rad_part = 0.5e-6_uRp
    real(kind = uRp), parameter :: vol_part = (4.0_uRp/3.0_uRp)*pi*rad_part**3
    real(kind = uRp), parameter :: rho_part = 2000
    real(kind = uRp), parameter :: mass_part = rho_part*vol_part

    !! Gas properties !!
    real(kind = uRp), parameter :: temp_gas = 300
    real(kind = uRp), parameter :: molar_m_gas = 40.0e-3_uRp ! Argon (Ar) molar mass 40 gram
    real(kind = uRp), parameter :: m_gas_cule = molar_m_gas/real(av_num, kind = uRp)
    real(kind = uRp), parameter :: mean_speed_gas = sqrt(8.0_uRp*kb*temp_gas/(pi*m_gas_cule)) ! <v> from M.B. distribution

end module parameters

program bro_force
    use parameters
    implicit none
    !! Simulation Variables !!
    real(kind = uRp), allocatable, dimension(:) :: pos, vel, acc
    !! Scenario Parameters !!
    ! integer, parameter                          :: num_part = 7*10**3
    integer, parameter                          :: num_part = 2*10**3
    !! Simulation Parameters !!
    ! integer, parameter                          :: num_steps = 10**6
    integer, parameter                          :: num_steps = 10**5
    real(kind = uRp), parameter                 :: del_time = 1.0e-8_uRp
    real(kind = uRp), parameter                 :: tot_time = num_steps*del_time
    
        ! real(kind = uRp), parameter                 :: del_time = 1.0e-8_uRp
        ! real(kind = uRp), parameter                 :: tot_time = 1.0e-3_uRp
        ! integer, parameter                          :: num_steps = ceiling(tot_time/del_time)

    logical, parameter                          :: an_int = .true. ! Analytical integration
    !! Browninan Parameters !!
    real(kind = uRp), parameter :: a = 0.50_uRp ! diffuse reflection
    real(kind = uRp), parameter :: den_gas = 0.42e-1_uRp
    real(kind = uRp), parameter :: fric_fac = (4.0_uRp/3.0_uRp)*pi*(rad_part**2)*den_gas*mean_speed_gas*(1.0_uRp + (pi*a/8.0_uRp))
    real(kind = uRp), parameter :: brown_int = fric_fac/mass_part
    real(kind = uRp), parameter :: s = 2.0_uRp*kb*temp_gas*brown_int/(pi*mass_part)
    real(kind = uRp), parameter :: s1 = sqrt(pi*s/del_time)

    real(kind = uRp), dimension(1:2) :: u, g
    real(kind = uRp) :: time
    !! Lin Reg !!
    real(kind = uRp), dimension(:), allocatable :: x_lin, y_lin
    !! Dummy Variables !!
    integer :: i, j
    integer, allocatable, dimension(:) :: seed

    allocate( pos(1:num_part), vel(1:num_part), acc(1:num_part) )
    pos = 0.0_uRp; vel = 0.0_uRp; acc = 0.0_uRp;
    
    !! Initialize the randomness - random_seed() etc
    call random_seed(size=i)
    allocate(seed(1:i))
    seed = (/ (j, j=1,i) /)
    call random_seed(put=seed)
    call random_seed(get=seed)
    write(*,*) seed
    write(*,*) i
    ! STOP 'After writing seed'

    time_loop:do i = 1,num_steps
        time = del_time*i
        part_loop:do j = 1, num_part
            !! Random for each particle !!
            call random_number(u)
            g = sqrt( -2.0_uRp*log(u(1)) ) * (/ cos( 2.0_uRp*pi*u(2) ), sin( 2.0_uRp*pi*u(2) ) /) 
            !! Time integration - Analytical or Numerical (Euler) !!
            if(an_int) then
                pos(j) = pos(j) + (mass_part**2)*g(1)*s1*(exp(-fric_fac*del_time/mass_part) - 1.0_uRp)/(fric_fac**2) &
                              & + mass_part*g(1)*s1*del_time/fric_fac &
                              & - mass_part*vel(j)*(exp(-fric_fac*del_time/mass_part) - 1.0_uRp)/fric_fac
                
                vel(j) = (vel(j) - mass_part*g(1)*s1/fric_fac)*exp(-fric_fac*del_time/mass_part) &
                               & + mass_part*g(1)*s1/fric_fac
            else
                acc(j) = g(1)*s1 - fric_fac*vel(j)/mass_part
                pos(j) = pos(j) + vel(j)*del_time + acc(j)*(del_time**2)
                vel(j) = vel(j) + acc(j)*del_time
            end if
        end do part_loop
        !! making time and mds arrays for linear regression y = m*x + c
        call append(time, x_lin)
        call append(sum(pos**2)/size(pos), y_lin)
        if( mod(i, num_steps/10) == 0) write(*,'(I0XA)',advance='no') int(i*100/num_steps), '% ' 
    end do time_loop
    write(*,*)
    write(*,'(AXES7.1)') 'Number of particles =', real(num_part)
    write(*,'(AXES8.2)') 'Time step =', del_time
    write(*,'(AXES7.1)') 'Number of time steps =', real(num_steps)
    write(*,'(AXL)') 'Analytical integration is', an_int

    u = lin_reg(x=x_lin, y=y_lin, report=.true.)
    write(*,'(AXES11.4XAXES11.4)') 'Linear fit with slope, m =', u(1), 'and y-intercept, c =', u(2)
    write(*,'(AXES11.4)') 'D = m/2 =', u(1)/2.0_uRp

    contains
    
    subroutine append(element,  array)
        implicit none
        real(kind = uRp), intent(in)  :: element
        real(kind = uRp), dimension(:), allocatable, intent(inout) :: array
        real(kind = uRp), dimension(:), allocatable                :: copy_array
        integer :: lb, ub
        if( allocated(array) ) then
            lb = lbound(array,1) ; ub = ubound(array,1) ! get array bounds lower and upper
            allocate(copy_array(lb:ub+1))
            copy_array(lb:ub) = array(lb:ub)
            copy_array(ub+1) = element
            
            deallocate(array)
            call move_alloc(copy_array, array)
        else
            allocate(array(1))
            array(1) = element
        end if 
        ! array = [array, element]
    end subroutine append
    !! Least Squares linear regression y = m*x + c. Takes x, y arrays gives out(m, c) !!
    function lin_reg(x, y, report) result(out)
        implicit none
        real(kind = uRp), dimension(:), allocatable, intent(in) :: x, y
        logical, optional , intent(in) :: report
        real(kind = uRp) :: out(1:2)
        real(kind = uRp) :: sumx, sumy, sumx2, sumy2, sumxy
        real(kind = uRp) :: n, m, c, r

        sumx  = sum(x); sumx2 = sum(x**2)
        sumy  = sum(y); sumy2 = sum(y**2)
        sumxy = sum(x*y)

        n = size(x,1)

        m = (n*sumxy - sumx*sumy) / (n*sumx2 - sumx**2)
        c = (sumy*sumx2 - sumx*sumxy) / (n*sumx2 - sumx**2)

        out = (/ m, c /)

        if(report) then
            r = (sumxy - sumx*sumy/n ) / sqrt((sumx2 - sumx**2/n)*(sumy2 - sumy**2/n))
            write(*,'(AXES9.3)') 'Correlation r value for linear fit is', r
        end if

        return
    end function lin_reg

end program bro_force