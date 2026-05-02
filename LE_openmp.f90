program EspectroLyapunov
    use, intrinsic :: iso_fortran_env, only: qp=>real128
    use omp_lib
    implicit none
    integer , parameter :: N_equ  = 3
    integer , parameter :: N_equ2 = 12
    !------------------------------------------------------------------
    ! Parametros del modelo (solo lectura tras parametros())
    !------------------------------------------------------------------
    real(qp) :: a, b, c, d, e, ff, gg, h
    real(qp) :: dt, dr
    real(qp) :: R_max, R_min
    real(qp) :: x00, y00, z00
    integer  :: Ntrans, Nr, ntau, IO, RRR

    !------------------------------------------------------------------
    ! Variables de control
    !------------------------------------------------------------------
    integer  :: g
    real(qp) :: RR_val   ! valor local del parametro para cada hilo

    !------------------------------------------------------------------
    ! Lectura de parametros y apertura de archivo de salida
    !------------------------------------------------------------------
    call parametros()

    RRR = int((R_max - R_min) / dr)
    open(3, file='Espectro1.dat')

    !------------------------------------------------------------------
    ! Loop paralelo sobre el parametro RR
    ! Cada iteracion es independiente
    !------------------------------------------------------------------
    !$omp parallel do                                                   &
    !$omp&   private(g, RR_val)                                         &
    !$omp&   schedule(dynamic)                                          &
    !$omp&   ordered
    do g = 0, RRR

        RR_val = R_min + dr * g
        call calcular_exponentes(RR_val, g)

    enddo
    !$omp end parallel do

    close(3)

!**********************************************************************
contains
!**********************************************************************

    subroutine calcular_exponentes(RR_loc, idx)
        real(qp), intent(in) :: RR_loc
        integer,  intent(in) :: idx

        !--------------------------------------------------------------
        ! Variables locales — cada hilo tiene su propia copia
        !--------------------------------------------------------------
        real(qp) :: yL(N_equ), y(N_equ2)
        real(qp) :: cum(N_equ), cum_traza, traza
        real(qp) :: v1(N_equ), v2(N_equ), v3(N_equ)
        real(qp) :: V11(N_equ), V22(N_equ), V33(N_equ)
        real(qp) :: NN(N_equ)
        real(qp) :: t
        integer  :: i, j, k, m

        !--------------------------------------------------------------
        ! Transitorio
        !--------------------------------------------------------------
        yL = [x00, y00, z00]
        t  = 0._qp

        do i = 1, Ntrans
            yL = yL + rk4_3(yL, t, dt, RR_loc)
            t  = t + dt
        enddo

        !--------------------------------------------------------------
        ! Condiciones iniciales para el calculo de exponentes
        !--------------------------------------------------------------
        y(1) = yL(1)
        y(2) = yL(2)
        y(3) = yL(3)

        do k = 4, N_equ2
            y(k) = 0._qp
        enddo
        y(4)  = 1._qp   ! M = Identidad
        y(8)  = 1._qp
        y(12) = 1._qp

        cum       = 0._qp
        cum_traza = 0._qp
        t         = 0._qp

        !--------------------------------------------------------------
        ! Loop principal de calculo de exponentes
        !--------------------------------------------------------------
        do m = 1, Nr

            ! Integrar ntau pasos
            do j = 1, ntau
                y = y + rk4_12(y, t, dt, RR_loc)
                t = t + dt
            enddo

            ! Acumular traza del Jacobiano para verificacion
            traza = (1._qp - 2._qp*y(1) - a*y(2) - b*y(3))    &
                  + (RR_loc*(1._qp - 2._qp*y(2)) - d*y(1))     &
                  + (e*y(1)/(y(1)+ff) - gg*y(1) - h)
            cum_traza = cum_traza + traza * dt * ntau

            ! Extraer vectores tangentes de la matriz M
            do k = 1, 3
                v1(k) = y(3*k+1)
                v2(k) = y(3*k+2)
                v3(k) = y(3*k+3)
            enddo

            ! Ortonormalizacion de Gram-Schmidt
            call OGS(v1, v2, v3, V11, V22, V33, NN)

            ! Acumular logaritmos de las normas
            cum = cum + log(NN)

            ! Imprimir progreso (solo hilo 0 para no mezclar salidas)
            if (omp_get_thread_num() == 0) then
                if (mod(m, IO) == 0 .or. m == Nr) then
                    print '(A,F6.3,A,F10.5,A,F10.5,A,F10.5)', &
                        'r2=', RR_loc,                          &
                        '  L1:', cum(1)/t,                      &
                        '  L2:', cum(2)/t,                      &
                        '  L3:', cum(3)/t
                endif
            endif

            ! Reiniciar M con vectores ortonormales
            do k = 1, 3
                y(3*k+1) = V11(k)
                y(3*k+2) = V22(k)
                y(3*k+3) = V33(k)
            enddo

        enddo

        !--------------------------------------------------------------
        ! Verificacion: suma lambdas vs media traza
        !--------------------------------------------------------------
        if (omp_get_thread_num() == 0) then
            print '(A,F6.3,A,F10.5,A,F10.5,A,ES10.3)',          &
                'r2=', RR_loc,                                     &
                '  Suma L=', sum(cum)/t,                           &
                '  Traza= ', cum_traza/t,                          &
                '  Dif=',    abs(sum(cum)/t - cum_traza/t)
        endif

        !--------------------------------------------------------------
        ! Escritura ordenada al archivo (evita mezcla entre hilos)
        !--------------------------------------------------------------
        !$omp ordered
        write(3,*) RR_loc, cum/t
        !$omp end ordered

    end subroutine calcular_exponentes

!**********************************************************************
!  Funcion del sistema: 3 ecuaciones
!**********************************************************************
    pure function f3(r, t, RR_loc) result(res)
        real(qp), intent(in) :: r(N_equ), t, RR_loc
        real(qp) :: res(N_equ)
        real(qp) :: u, v, w

        u = r(1); v = r(2); w = r(3)

        res(1) = u*(1._qp - u) - a*u*v - b*u*w
        res(2) = RR_loc*v*(1._qp - v) - d*u*v
        res(3) = (e*u*w)/(u + ff) - gg*u*w - h*w

    end function f3

!**********************************************************************
!  Funcion del sistema ampliado: 12 ecuaciones
!**********************************************************************
    pure function f12(r, t, RR_loc) result(res)
        real(qp), intent(in) :: r(N_equ2), t, RR_loc
        real(qp) :: res(N_equ2)
        real(qp) :: u, v, w
        real(qp) :: J(3,3), M(3,3), P(3,3)
        integer  :: aa, bb, cc

        u = r(1); v = r(2); w = r(3)

        ! Ecuaciones del modelo
        res(1) = u*(1._qp - u) - a*u*v - b*u*w
        res(2) = RR_loc*v*(1._qp - v) - d*u*v
        res(3) = (e*u*w)/(u + ff) - gg*u*w - h*w

        ! Jacobiano — RR_loc pasado explicitamente, sin variable global
        J(1,1) = 1._qp - 2._qp*u - a*v - b*w
        J(1,2) = -a*u
        J(1,3) = -b*u

        J(2,1) = -d*v
        J(2,2) = RR_loc*(1._qp - 2._qp*v) - d*u
        J(2,3) = 0._qp

        J(3,1) = (e*w/(u+ff))*(1._qp - (u/(u+ff))) - gg*w
        J(3,2) = 0._qp
        J(3,3) = (e*u/(u+ff)) - gg*u - h

        ! Matriz de variacion M
        M(1,1) = r(4);  M(1,2) = r(5);  M(1,3) = r(6)
        M(2,1) = r(7);  M(2,2) = r(8);  M(2,3) = r(9)
        M(3,1) = r(10); M(3,2) = r(11); M(3,3) = r(12)

        ! P = J * M
        do aa = 1, 3
            do bb = 1, 3
                P(aa,bb) = 0._qp
                do cc = 1, 3
                    P(aa,bb) = P(aa,bb) + J(aa,cc)*M(cc,bb)
                enddo
            enddo
        enddo

        res(4)  = P(1,1); res(5)  = P(1,2); res(6)  = P(1,3)
        res(7)  = P(2,1); res(8)  = P(2,2); res(9)  = P(2,3)
        res(10) = P(3,1); res(11) = P(3,2); res(12) = P(3,3)

    end function f12

!**********************************************************************
!  RK4 para 3 ecuaciones (transitorio)
!**********************************************************************
    pure function rk4_3(r, t, dt, RR_loc)
        real(qp), intent(in) :: r(N_equ), t, dt, RR_loc
        real(qp) :: rk4_3(N_equ)
        real(qp) :: k1(N_equ), k2(N_equ), k3(N_equ), k4(N_equ)

        k1 = dt * f3(r,               t,             RR_loc)
        k2 = dt * f3(r + 0.5_qp*k1,  t + 0.5_qp*dt, RR_loc)
        k3 = dt * f3(r + 0.5_qp*k2,  t + 0.5_qp*dt, RR_loc)
        k4 = dt * f3(r + k3,          t + dt,         RR_loc)

        rk4_3 = (k1 + 2._qp*k2 + 2._qp*k3 + k4) / 6._qp
    end function rk4_3

!**********************************************************************
!  RK4 para 12 ecuaciones (calculo de exponentes)
!**********************************************************************
    pure function rk4_12(r, t, dt, RR_loc)
        real(qp), intent(in) :: r(N_equ2), t, dt, RR_loc
        real(qp) :: rk4_12(N_equ2)
        real(qp) :: k1(N_equ2), k2(N_equ2), k3(N_equ2), k4(N_equ2)

        k1 = dt * f12(r,               t,             RR_loc)
        k2 = dt * f12(r + 0.5_qp*k1,  t + 0.5_qp*dt, RR_loc)
        k3 = dt * f12(r + 0.5_qp*k2,  t + 0.5_qp*dt, RR_loc)
        k4 = dt * f12(r + k3,          t + dt,         RR_loc)

        rk4_12 = (k1 + 2._qp*k2 + 2._qp*k3 + k4) / 6._qp
    end function rk4_12

!**********************************************************************
!  Ortonormalizacion de Gram-Schmidt
!**********************************************************************
    subroutine OGS(v1, v2, v3, V11, V22, V33, NN)
        real(qp), intent(in)  :: v1(N_equ), v2(N_equ), v3(N_equ)
        real(qp), intent(out) :: V11(N_equ), V22(N_equ), V33(N_equ)
        real(qp), intent(out) :: NN(N_equ)
        real(qp) :: Coef21, Coef31, Coef32

        NN(1) = norma(v1)
        V11   = v1 / NN(1)

        Coef21 = producto_punto(v2, V11)
        V22    = v2 - Coef21*V11
        NN(2)  = norma(V22)
        V22    = V22 / NN(2)

        Coef31 = producto_punto(v3, V11)
        Coef32 = producto_punto(v3, V22)
        V33    = v3 - Coef31*V11 - Coef32*V22
        NN(3)  = norma(V33)
        V33    = V33 / NN(3)

    end subroutine OGS

!**********************************************************************
!  Lectura de parametros desde archivo
!**********************************************************************
    subroutine parametros()
        integer :: unit_par
        unit_par = 10
        open(unit_par, file='LE_c1.par', status='old')
        read(unit_par,*) a, b, c, d, e, ff, gg, h
        read(unit_par,*) R_min, R_max
        read(unit_par,*) dt, dr
        read(unit_par,*) Nr, ntau
        read(unit_par,*) IO
        read(unit_par,*) x00, y00, z00
        read(unit_par,*) Ntrans
        close(unit_par)

        print '(A,F10.5,A,F10.5,A,F10.5)', 'a=',a,  ' b=',b,  ' c=',c
        print '(A,F10.5,A,F10.5,A,F10.5)', 'd=',d,  ' e=',e,  ' f=',ff
        print '(A,F10.5,A,F10.5)',          'g=',gg, ' h=',h
        print '(A,F10.5,A,F10.5)',          'R_min=',R_min, ' R_max=',R_max
        print '(A,F10.5,A,F10.5)',          'dt=',dt, ' dr=',dr
        print '(A,I10,A,I10)',              'Nr=',Nr, ' ntau=',ntau
        print '(A,I10)',                    'IO=',IO
        print '(A,F10.5,A,F10.5,A,F10.5)', 'x0=',x00,' y0=',y00,' z0=',z00
        print '(A,I10)',                    'Ntrans=',Ntrans
        print '(A,I4,A)',                   'Hilos disponibles: ', &
                                             omp_get_max_threads(), ' (OpenMP)'
    end subroutine parametros

!**********************************************************************
!  Norma euclidiana
!**********************************************************************
    pure function norma(vector)
        real(qp), intent(in) :: vector(N_equ)
        real(qp) :: norma
        integer  :: zz
        norma = 0._qp
        do zz = 1, N_equ
            norma = norma + vector(zz)**2
        enddo
        norma = sqrt(norma)
    end function norma

!**********************************************************************
!  Producto punto
!**********************************************************************
    pure function producto_punto(vector1, vector2)
        real(qp), intent(in) :: vector1(N_equ), vector2(N_equ)
        real(qp) :: producto_punto
        integer  :: zzz
        producto_punto = 0._qp
        do zzz = 1, N_equ
            producto_punto = producto_punto + vector1(zzz)*vector2(zzz)
        enddo
    end function producto_punto

!**********************************************************************
end program EspectroLyapunov
