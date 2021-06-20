program blasius
    use parameters, only: wp => PARM_DP, intlength => PARM_SI, charlength => PARM_SCL
    implicit none

    ! Declaring the integer variables
    integer(intlength) :: eta_points, unidadeleitura, unidadeescrita, rc1, rc2, unidadeescrita2
    ! Declarando variaveis da linha de comando
    character(charlength) :: input_name, output_name
    ! Declaring real variables
    real(wp) :: eta, eta_initial, eta_final, eta_step, f, df, d2f, f_new, &
    df_new, d2f_new, tol_1, resultado(3)
    real(wp), allocatable :: blasius_parameters(:) 
    
    call get_command_argument(1,input_name)     ! Entrada
    call get_command_argument(2,output_name)    ! Saida

    open(newunit=unidadeleitura, & ! Abrir o arquivo csv pela linha de comando
    file=trim(input_name), &
    action='read', &
    iostat=rc1 &
    )

    open(newunit=unidadeescrita, & ! Escrever neste arquivo chamado na linha de comando
    status = 'replace', &
    file=trim(output_name), &
    action = 'write', &
    position = 'rewind', &
    iostat = rc2 &
    )

    ! Condicoes iniciais do sistema
    call read_properties
    if (rc1 .ne. 0) then     
        print *, "The file has not been read correctly."
        stop
    end if
    
    eta_points = 1 + int((eta_final-eta_initial)/eta_step)
    allocate(blasius_parameters(3))
    
    ! Using values of prior results to find real f''(0)
    eta = eta_initial
    blasius_parameters = newton_raphson_rk4(eta_initial,d2f,rk4(eta,(/f,df,d2f/)))

    deallocate(blasius_parameters)
    close(unidadeleitura)
    close(unidadeescrita)

    contains

    subroutine read_properties()
        read(unidadeleitura,*) ! Just the variable names
        read(unidadeleitura,*) eta_initial, eta_final, eta_step
        read(unidadeleitura,*) ! Just the variable names
        read(unidadeleitura,*) f, df, d2f, tol_1
    end subroutine read_properties
    
    function order_reduction_blasius_rk4(eta,bp) result(derivatives)
        real(wp), intent(in) :: eta, bp(3)
        real(wp) :: derivatives(3)

        derivatives(1) = bp(2)
        derivatives(2) = bp(3)
        derivatives(3) = -0.5d0*bp(1)*bp(3)

    end function order_reduction_blasius_rk4

    function rk4(eta,bp) result(solver)

        integer(intlength) :: iter=1
        real(wp) :: eta, eta_new, bp(3), solver(3)
        real(wp), dimension(3) :: c1, c2, c3, c4

        eta_new = eta
        write(unidadeescrita,'(a18,3x,a18,3x,a18,3x,a18)') "eta","f","f'","f''"
        do iter = 1, eta_points, 1
            write(unidadeescrita,'(f18.10,3x,f18.10,3x,f18.10,3x,f18.10)') eta_new, bp
            c1 = order_reduction_blasius_rk4(eta_new,bp)
            c2 = order_reduction_blasius_rk4(eta_new+0.5*eta_step,bp+0.5*eta_step*c1)
            c3 = order_reduction_blasius_rk4(eta_new+0.5*eta_step,bp+0.5*eta_step*c2)
            c4 = order_reduction_blasius_rk4(eta_new+1.0*eta_step,bp+1.0*eta_step*c3)
            bp = bp + eta_step*(c1+2.0d0*(c2+c3)+c4)/6.0d0
            eta_new = eta_new + eta_step
        end do

        solver = bp
    
    end function rk4

    function newton_raphson_rk4(eta, d2f, bp) result(bp_new)
        integer(intlength) :: iter=1, imax=100
        real(wp) :: d2f, eta, eta_new, f_1stderiv_eta, f_2ndderiv_eta, &
        f_2ndderiv_new, bp(3), bp_new(3)
        real(wp), dimension(3) :: c1, c2, c3, c4
        
        f_1stderiv_eta = bp(2)
        f_2ndderiv_eta = d2f

        ! 100 iterations to prevent infinite loop
        do iter = 1,eta_points,1
            if (abs(f_1stderiv_eta-1.0d0) < tol_1) stop
            ! Find the root
            f_2ndderiv_new = f_2ndderiv_eta - (f_1stderiv_eta-1.0d0)/finite_difference(0.0d0,bp)
            ! Using the root to calculate rk4 again with new value
            rewind(unidadeescrita)
            bp = rk4(0.0d0,(/0.0d0, 0.0d0, f_2ndderiv_new/))
            f_1stderiv_eta = bp(2)
            f_2ndderiv_eta = f_2ndderiv_new
        end do
        
        bp_new = bp

    end function newton_raphson_rk4

    function finite_difference(eta, bp) result(cent_diff_1st) 
        real(wp) :: eta, bp(3), cent_diff_1st, breg(3), bpro(3), &
        breg_new(3), bpro_new(3)
        
        breg = (/0.0d0, 0.0d0, bp(3) - tol_1/)
        bpro = (/0.0d0, 0.0d0, bp(3) + tol_1/)

        breg_new = rk4(eta,breg)
        bpro_new = rk4(eta,bpro)

        cent_diff_1st = (bpro_new(2)-breg_new(2))/(2.0d0*tol_1)
    end function finite_difference

end program blasius