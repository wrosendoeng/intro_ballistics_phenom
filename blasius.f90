program blasius
    use parameters, only: wp => PARM_DP, intlength => PARM_SI, charlength => PARM_SCL
    implicit none

    ! Declaring the integer variables
    integer(intlength) :: eta_points, unidadeleitura, unidadeescrita, rc1, rc2
    ! Declarando variaveis da linha de comando
    character(charlength) :: input_name, output_name
    ! Declaring real variables
    real(wp) :: eta_initial, eta_final, eta_step, f, df, d2f, f_new, df_new, d2f_new
    real(wp), allocatable :: blasius_parameters(:) 
    
    call get_command_argument(1,input_name)  ! Entrada
    call get_command_argument(2,output_name) ! Saida

    open(newunit=unidadeleitura, & ! Abrir o arquivo csv pela linha de comando
    file=trim(input_name), &
    action='read', &
    iostat=rc1 &
    )

    if (rc1 .ne. 0) then     
        print *, "The file has not been read correctly."
        stop
    end if

    open(newunit=unidadeescrita, & ! Escrever neste arquivo chamado na linha de comando
    status = 'replace', &
    file = trim(output_name), &
    action = 'write', &
    iostat = rc2 &
    )

    allocate(blasius_parameters(3))

    ! Condicoes iniciais do sistema
    call read_properties
    
    ! Preparando um vetor que armazene todos os dados para realizar o RK4
    blasius_parameters = (/f,df,d2f/)
    eta_points = 1 + int((eta_final-eta_initial)/eta_step)

    ! Solving considering ficticious f''(0)
    call rk4(eta_initial,blasius_parameters) 

    ! Using values of prior results to find real f''(0)
    d2f_new = newton_raphson(blasius_parameters)

    ! Recreating the blasius parameters
    f_new = 0.0e0; df_new = f_new
    blasius_parameters = (/f_new,df_new,d2f_new/)
    
    ! Solving considering with real f''(0)
    call rk4(eta_initial,blasius_parameters)

    allocate(blasius_parameters(3))
    close(unidadeleitura)
    close(unidadeescrita)

    contains

    subroutine read_properties()
        read(unidadeleitura,*) ! Just the variable names
        read(unidadeleitura,*) eta, eta_step, f, df, d2f, tol_1
    end subroutine read_properties
    
    function order_reduction_blasius_rk4(eta,blasius_parameters) result(derivatives)
        
        ! Input arguments
        real(wp), intent(in) :: eta, blasius_parameters(3) 
        ! f_rk4, df_rk4, d2f_rk4
        real(wp) :: derivatives(3), f_rk4, df_rk4, d2f_rk4
    
        ! Separando os termos para a posicao e velocidade
        f_rk4 = blasius_parameters(1)
        df_rk4 = blasius_parameters(2)
        d2f_rk4 = blasius_parameters(3)

        derivatives(1) = df_rk4
        derivatives(2) = d2f_rk4
        derivatives(3) = -0.5e0*f_rk4*d2f_rk4

    end function

    subroutine rk4(eta,blasius_parameters)
        
        integer(intlength) :: eta_points
        real(wp), intent(in) :: eta
        real(wp), intent(inout) :: blasius_parameters(3)
        real(wp), dimension(3) :: const1, const2, const3, const4

        do eta_points = 1,eta_points,1
            const1 = order_reduction_blasius_rk4(eta,blasius_parameters)
            const2 = order_reduction_blasius_rk4(eta+0.5*eta_step,blasius_parameters+0.5*eta_step*const1)
            const3 = order_reduction_blasius_rk4(eta+0.5*eta_step,blasius_parameters+0.5*eta_step*const2)
            const4 = order_reduction_blasius_rk4(eta+1.0*eta_step,blasius_parameters+1.0*eta_step*const3)

            blasius_parameters = blasius_parameters + eta_step/6.0e0*(const1+2.0e0*(const2+const3)+const4)
            !write(unidadeescrita,*) eta, blasius_parameters
            eta = eta + eta_step	
        end do
    
    end subroutine rk4

    function newton_raphson(blasius_parameters) result(x1)
        integer(intlength) :: i, imax=100
        real(wp) :: f_1stderiv_eta, f_2ndderiv_eta, f_2ndderiv_new, support(3)
        
        f_1stderiv_eta = blasius_parameters(2)
        f_2ndderiv_eta = blasius_parameters(3)

        ! 100 iterations to prevent infinite loop
        do while (i < imax .or. abs(f_1stderiv_eta-1.0) > tol_1) 
            f_2ndderiv_new = f_2ndderiv_eta - (f_1stderiv_eta-1.0e0)/cent_diff_1st(blasius_parameters) !COMO CALCULAR A DIFERENCA FINITA ??
            support = (/0.0e0,0.0e0,f_2ndderiv_new/)
            call rk4(eta_initial,support)
            f_1stderiv_eta = support(2)
            f_2ndderiv_eta = support(3)
            i = i+1
        end do

    end function newton_raphson

    function finite_difference(blasius_parameters) result(cent_diff_1st) 
        real(wp) :: cent_diff_1st, f2, f1, blasius_progressive, blasius_regressive
        
        !Create numerical method

        blasius_progressive(3) = blasius_parameters(3) + tol_1
        blasius_regressive(3) = blasius_parameters(3) - tol_1
        
        call rk4(eta_initial,blasius_progressive)
        call rk4(eta_initial,blasius_regressive)

        f1 = blasius_parameters1(1)
        f2 = blasius_parameters2(1)
        
        cent_diff_1st = (f2-f1)/(2.0e0*tol_1)
    end function finite_difference 

end program blasius