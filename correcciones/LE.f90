program EspectroLyapunov
    use, intrinsic :: iso_fortran_env, only: qp=>real128
    implicit none
    integer , parameter :: N_equ = 3    ! Numero de ecuaciones
    integer , parameter :: N_equ2 = 12
    real(qp) :: a,b,c,d,e,ff,gg,h, dt, dr, t, x0, y0, z0, Ntot, RR
    real(qp) :: R_max,R_min
    real(qp) :: x00, y00, z00
    real(qp) :: cum(N_equ),yL(N_equ),y(N_equ2)
    real(qp) :: v1(N_equ),v2(N_equ),v3(N_equ),V11(N_equ),V22(N_equ)
    real(qp) :: V33(N_equ),N(N_equ)
    integer  :: i, j, k, Ntrans, m, Nr, ntau, IO, rrr, g
    
    real(qp) :: cum_traza, traza
!**********************************************************************
    call parametros()                           ! entrada de parametros
!**********************************************************************
    Ntot=Nr*ntau                                          ! transitorio
!    open(1,file='LE_trans1.dat') !-
    RRR=int((R_max-R_min)/dr)
    open(3,file='Espectro1.dat')
    do g = 0, RRR
    yL = [ x00 , y00 , z00 ]
    t=0._qp
    RR=R_min + dr*g
        do i = 1, Ntrans
            yL = yL + rk4( yL, t, dt, N_equ )
            t = t + dt
!            write(1,*) t, yL !-
        end do
!        close(1)           !-
    !    call system('gnuplot -c LE_trans.p')
        x0=yL(1)                                  ! Fin del transitorio
        y0=yL(2)
        z0=yL(3)       ! Guardando valores finales
        t=0._qp
!        print '(AX,F10.5x,F10.5,F10.5,F10.5)','iniciaremos el calculo en',x0,y0,z0,t
!**********************************************************************
        y(1)=x0                                 ! Calculo de Exponentes
        y(2)=y0  ! valores iniciales
        y(3)=z0
        cum=0._qp
        cum_traza=0._qp
        do k=4, N_equ2
            y(k)=0._qp
        enddo
        y(4)=1._qp; y(8)=1._qp; y(12)=1._qp ! Cambiando a sis de 12 dim

!        open(2,file='LyaExp1.dat')    !- Avanzando ntau pasos con rk4
        do m=1,Nr
            do j=1,ntau
                y = y + rk4( y, t, dt, N_equ2 )
                t = t + dt
            enddo
            
    ! Calcular y acumular tr(J) en el punto actual
    traza = (1._qp - 2._qp*y(1) - a*y(2) - b*y(3)) &
          + (RR*(1._qp - 2._qp*y(2)) - d*y(1))      &
          + (e*y(1)/(y(1)+ff) - gg*y(1) - h)
    cum_traza = cum_traza + traza * dt * ntau

            do k=1,3
                v1(k)=y(3*k+1)  ! v1(1)=y(4)  v1(2)=y(7)  v1(3)=y(10)
                v2(k)=y(3*k+2)  ! v2(1)=y(5)  v2(2)=y(8)  v1(3)=y(11)
                v3(k)=y(3*k+3)  ! v3(1)=y(6)  v3(2)=y(9)  v3(3)=y(12)
            enddo
            call OGS(v1,v2,v3,V11,V22,V33,N)  ! ortonormalizacion de GS
            !print'(F10.5X,F10.5X,F10.5)',N
            cum = cum + log(N)
            if(mod(m,IO) .eq. 0 .or. m .eq. Nr) then  
!                write(2,*) t-dt, cum/t            
              !  print   *, t-dt, cum/t            
            endif                                 

            do k=1,3
               y(3*k+1)=V11(k) ! y(4)=V11(1)  y(7)=V11(2)  y(10)=V11(3)
               y(3*k+2)=V22(k) ! y(5)=V22(1)  y(8)=V22(2)  y(11)=V22(3)
               y(3*k+3)=V33(k) ! y(6)=V33(1)  y(9)=V33(2)  y(12)=V33(3)
            enddo
        enddo
    print'(A,F10.5,F10.5,F10.5)','Exponentes',cum/t
    print'(A,F10.5)','Suma lambdas    = ', sum(cum)/t
    print'(A,F10.5)','Media de traza  = ', cum_traza/t
    print'(A,F10.5)','Diferencia      = ', abs(sum(cum)/t - cum_traza/t)
!    close(2)           !-
    write(3,*) rr, cum/t
    enddo
    close(3)
   ! call system('gnuplot -c LyaExp.p')
!    call system('gnuplot -c Espectro.p')
!**********************************************************************
contains
!**********************************************************************
    pure function f(r, t, Nequ)
        real(qp), intent(in) :: r(Nequ) ! Valores
        real(qp), intent(in) :: t    ! Paso
        integer , intent(in) :: Nequ
        real(qp) :: f(Nequ)
        real(qp) :: u, v, w
        real(qp) :: J(3,3),M(3,3),P(3,3)
        integer  :: aa, bb, cc

        u = r(1)
        v = r(2)
        w = r(3)

  f(1) = u*(1.0-u) - a*u*v - b*u*w
  f(2) = RR*v*(1.0-v) - d*u*v
  f(3) = (e*u*w)/(u+ff) - gg*u*w - h*w

        if(Nequ.gt.3) then

        J(1,1) = 1._qp - 2._qp*u - a*v - b*w;            J(1,2)= -a*u;                       J(1,3)= -b*u
        J(2,1) = -d*v;                                   J(2,2)= RR*(1._qp - 2._qp*v) - d*u; J(2,3)= 0.0_qp
        J(3,1) = (e*w/(u+ff))*(1._qp-(u/(u+ff))) - gg*w; J(3,2)= 0._qp;                      J(3,3)= (e*u/(u+ff)) - gg*u - h

        M(1,1) =  r(4); M(1,2) =  r(5); M(1,3) =  r(6)
        M(2,1) =  r(7); M(2,2) =  r(8); M(2,3) =  r(9)
        M(3,1) = r(10); M(3,2) = r(11); M(3,3) = r(12)

        do aa=1,3
            do bb=1,3
                P(aa,bb) = 0._qp
                do cc=1,3
                    P(aa,bb) = P(aa,bb) + J(aa,cc)*M(cc,bb)
                enddo
            enddo
        enddo

         f(4) = P(1,1);  f(5) = P(1,2);  f(6) = P(1,3)
         f(7) = P(2,1);  f(8) = P(2,2);  f(9) = P(2,3)
        f(10) = P(3,1); f(11) = P(3,2); f(12) = P(3,3)

        endif

    end function f
!**********************************************************************
    pure function rk4(r, t, dt, Nequ)                   ! Runge-Kutta 4
        real(qp), intent(in) :: r(Nequ) ! Valores
        real(qp), intent(in) :: t    ! Paso
        real(qp), intent(in) :: dt   ! Tamano de paso
        integer , intent(in) :: Nequ
        real(qp)  :: rk4(Nequ)
        real(qp)  :: k1(Nequ), k2(Nequ), k3(Nequ), k4(Nequ)

        k1 = dt * f( r              , t              , Nequ )
        k2 = dt * f( r + 0.5_qp * k1, t + 0.5_qp * dt, Nequ )
        k3 = dt * f( r + 0.5_qp * k2, t + 0.5_qp * dt, Nequ )
        k4 = dt * f( r + k3         , t + dt         , Nequ )

        rk4 = ( k1 + ( 2._qp * k2 ) + ( 2._qp * k3 ) + k4 ) / 6._qp
    end function rk4
!**********************************************************************
    subroutine OGS(v1,v2,v3,V11,V22,V33,N)
        real(qp), intent(in)    :: v1(N_equ),v2(N_equ),v3(N_equ)
        real(qp), intent(out)   :: V11(N_equ),V22(N_equ),V33(N_equ)
        real(qp), intent(inout) :: N(N_equ)
        real(qp) :: Coef21, Coef31, Coef32

        N(1) = norma(v1)
        V11 = v1/N(1)

        Coef21 = producto_punto(v2,V11)
        V22 = v2 - Coef21*V11
        N(2) = norma(V22)
        V22 = V22/N(2)

        Coef31 = producto_punto(v3,V11)
        Coef32 = producto_punto(v3,V22)
        V33 = v3 - Coef31*V11 - Coef32*V22
        N(3) = norma(V33)
        V33 = V33/N(3)
    end subroutine OGS
!**********************************************************************
    subroutine parametros()
    integer :: unit_par
    unit_par = 10
    open(unit_par, file='LE_c1.par', status='old')
    read(unit_par,*) a,b,c,d,e,ff,gg,h
    read(unit_par,*) R_min, R_max
    read(unit_par,*) dt, dr
    read(unit_par,*) Nr, ntau
    read(unit_par,*) IO
    read(unit_par,*) x00, y00, z00
    read(unit_par,*) Ntrans
    close(unit_par)

        print '(A,F10.5,A,F10.5,A,F10.5)','a',a,'b',b,'c',c
        print '(A,F10.5,A,F10.5,A,F10.5)','d',d,'e',e,'f',ff
        print '(A,F10.5,A,F10.5)','g',gg,'h',h
        print '(A,F10.5,A,F10.5)', 'R_min', R_min, 'R_max', R_max
        print '(A,F10.5,A,F10.5)', 'dt', dt, 'dr', dr
        print '(A,I20,A,I20)', 'Nr',Nr,'ntau', ntau
        print '(A,I20)','IO', IO
        print '(A,F10.5,A,F10.5,A,F10.5)','x0',x00,'y0',y00,'z0',z00
        print '(A,I20)','ntrans',Ntrans
    end subroutine parametros
!**********************************************************************
    pure function norma(vector)
        real(qp), intent(in) :: vector(N_equ)
        real(qp) :: norma
        integer :: zz
        norma=0
        do zz =1,N_equ
            norma = norma + vector(zz)**2
        enddo
        norma = sqrt(norma)
    end function norma
!**********************************************************************
    pure function producto_punto(vector1,vector2)
        real(qp), intent(in) :: vector1(N_equ),vector2(N_equ)
        real(qp) :: producto_punto
        integer  :: zzz
        producto_punto=0
        do zzz=1,N_equ
            producto_punto = producto_punto + vector1(zzz)*vector2(zzz)
        enddo
    end function producto_punto
!**********************************************************************
end program EspectroLyapunov
