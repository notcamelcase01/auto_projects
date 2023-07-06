! pure bending
    SUBROUTINE FUNC(NDIM,  U, ICP, PAR, IJAC, F, DFDU, DFDP)
    IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
        DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
        DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM, NDIM), DFDP(NDIM,*)
        DOUBLE PRECISION, DIMENSION(3) :: v, k, v1, k1, v2, k2, r, n, m, neff, meff, nbar, mbar
        DOUBLE PRECISION, DIMENSION(4) :: q
        DOUBLE PRECISION, DIMENSION(4, 3) :: A
        DOUBLE PRECISION, DIMENSION(3, 3) :: R0
        INTEGER i, j
        DOUBLE PRECISION C1, C2, C3, CC1, CC2, CC3, D1, D2, D3, DD1, DD2, DD3, l, E, nu, G, dia, MO, ar, Ii, D, Jj, kappa, pi
        pi = 22.0 / 7.0
		dia = 0.05
        l = 0.05
        E = 10.0 ** 6
        nu = 0.0
        Ar = 0.25 * pi * dia ** 2
        Ii = pi * dia ** 4 / 64
        Jj = 2 * Ii
  	    G = E / (2 * (1 + nu))
        C1 = G * Ar
        C2 = G * Ar
        C3 = E * Ar
        
        CC1 = l * l * G * Ar
        CC2 = l * l * G * Ar
        CC3 = l * l * E * Ar

        D1 = E * Ii + l * l * E * Ar
        D2 = E * Ii + l * l * E * Ar
        D3 = G * Jj + 2 * l * l * G * Ar

        DD1 = E * Ii * l * l
        DD2 = E * Ii * l * l
        DD3 = G * Jj * l * l

        DO i = 1, 3
            v1(i) = U(i)
            k1(i) = U(3 + i)
            r(i) = U(6 + i)
            q(i) = U(9 + i)
            v(i) = U(13 + i)
            k(i) = U(16 + i)
            v2(i) = U(19 + i)
            k2(i) = U(22 + i)
        END DO
        q(4) = U(13)

        R0(1, 1) = 2.0 * (q(1) * q(1) + q(2) * q(2) - 0.5)
        R0(2, 2) = 2.0 * (q(1) * q(1) + q(3) * q(3) - 0.5)
        R0(3, 3) = 2.0 * (q(1) * q(1) + q(4) * q(4) - 0.5)
        R0(1, 2) = 2.0 * (q(2) * q(3) - q(1) * q(4))
        R0(2, 1) = 2.0 * (q(2) * q(3) + q(1) * q(4))
        R0(1, 3) = 2.0 * (q(2) * q(4) + q(1) * q(3))
        R0(3, 1) = 2.0 * (q(2) * q(4) - q(1) * q(3))
        R0(2, 3) = 2.0 * (q(3) * q(4) - q(1) * q(2))
        R0(3, 2) = 2.0 * (q(3) * q(4) + q(1) * q(2))

        A(2, 1) =  0.5 * q(1)
        A(3, 2) =  0.5 * q(1)
        A(4, 3) =  0.5 * q(1)
        A(1, 1) = -0.5 * q(2)
        A(3, 3) = -0.5 * q(2)
        A(4, 2) =  0.5 * q(2)
        A(1, 2) = -0.5 * q(3)
        A(2, 3) =  0.5 * q(3)
        A(4, 1) = -0.5 * q(3)
        A(1, 3) = -0.5 * q(4)
        A(2, 2) = -0.5 * q(4)
        A(3, 1) =  0.5 * q(4)

        F(1) = v2(1)
        F(2) = v2(2)
        F(3) = v2(3)

        F(4) = k2(1)
        F(5) = k2(2)
        F(6) = k2(3)

        DO i = 1, 3
            F(6 + i) = 0.0
            F(9 + i) = PAR(2) * q(i)
            DO j = 1, 3
                F(6 + i) = F(6 + i) + R0(i, j) * v(j)
                F(9 + i) = F(9 + i) + A(i, j) * k(j)
            END DO
        END DO

        F(13) = PAR(2) * q(4)

        DO j = 1, 3
            F(13) = F(13) + A(4, j) * k(j)
        END DO

        F(14) = v1(1)
        F(15) = v1(2)
        F(16) = v1(3)

        F(17) = k1(1)
        F(18) = k1(2)
        F(19) = k1(3)

        n(1) = C1 * v(1)
        n(2) = C2 * v(2)
        n(3) = C3 * (v(3) - 1)

        m(1) = D1 * k(1)
        m(2) = D2 * k(2)
        m(3) = D3 * k(3)

        nbar(1) = CC1 * v2(1)
        nbar(2) = CC2 * v2(2)
        nbar(3) = CC3 * v2(3)

        mbar(1) = DD1 * k2(1)
        mbar(2) = DD2 * k2(2)
        mbar(3) = DD3 * k2(3)

        DO i = 1, 3
            neff(i) = n(i) - nbar(i)
            meff(i) = m(i) - mbar(i)
        END DO

        F(20) = (C1 * v1(1) - (neff(2) * k(3) - neff(3) * k(2))) / CC1
        F(21) = (C2 * v1(2) - (neff(3) * k(1) - neff(1) * k(3))) / CC2
        F(22) = (C3 * v1(3) - (neff(1) * k(2) - neff(2) * k(1))) / CC3

        F(23) = (D1 * k1(1) - (meff(2) * k(3) - meff(3) * k(2) + neff(2) * v(3) - neff(3) * v(2))) / DD1
        F(24) = (D2 * k1(2) - (meff(3) * k(1) - meff(1) * k(3) + neff(3) * v(1) - neff(1) * v(3))) / DD2
        F(25) = (D3 * k1(3) - (meff(1) * k(2) - meff(2) * k(1) + neff(1) * v(2) - neff(2) * v(1))) / DD3
    END SUBROUTINE FUNC

! ------------------------------------------------------------ !

    SUBROUTINE STPNT(NDIM, U, PAR, T)
    IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
        DOUBLE PRECISION, INTENT(IN) :: T
        INTEGER i
        PAR(1) = 0.0 ! load
        PAR(2) = 0.0 
        DO i = 1, 25
            U(i) = 0.0
        END DO
        U(10) = 1
        U(9) = 1
        U(16) = 1.0
    END SUBROUTINE STPNT

! ------------------------------------------------------------ !

    SUBROUTINE BCND(NDIM, PAR, ICP, NBC, U0, U1, FB, IJAC, DBC)
    IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
        DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
        DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
        DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
        DOUBLE PRECISION, DIMENSION(3, 3) :: R0
        DOUBLE PRECISION, DIMENSION(3) :: v, k, v1, k1, v2, k2, r,  n, m, neff, meff, nbar, mbar, C, CC, D, DD, rv, r_dash
        DOUBLE PRECISION, DIMENSION(4) :: q
        DOUBLE PRECISION, DIMENSION(4, 3) :: A

        INTEGER i, j
               DOUBLE PRECISION C1, C2, C3, CC1, CC2, CC3, D1, D2, D3, DD1, DD2, DD3, l, E, nu, G, dia, MO, Ar, Ii, Jj, kappa, pi
        pi = 22 / 7
		nu = 0.0
		dia = 0.05
        l = 0.05
        E = 10.0 ** 6
        Ar = 0.25 * pi * dia ** 2
        Ii = pi * dia ** 4 / 64
        Jj = 2 * Ii
  	    G = E / (2 * (1 + nu))
        C1 = G * Ar
        C2 = G * Ar
        C3 = E * Ar
        
        CC1 = l * l * G * Ar
        CC2 = l * l * G * Ar
        CC3 = l * l * E * Ar

        D1 = E * Ii + l * l * E * Ar
        D2 = E * Ii + l * l * E * Ar
        D3 = G * Jj + 2 * l * l * G * Ar

        DD1 = E * Ii * l * l
        DD2 = E * Ii * l * l
        DD3 = G * Jj * l * l
        C(1) = C1
        C(2) = C2
        C(3) = C3
        D(1) = D1
        D(2) = D2
        D(3) = D3
        CC(1) = CC1
        CC(2) = CC2
        CC(3) = CC3
        DD(1) = DD1
        DD(2) = DD2
        DD(3) = DD3
        
        DO i = 1, 3
            v1(i) = U0(i)
            k1(i) = U0(3 + i)
            r(i) = U0(6 + i)
            q(i) = U0(9 + i)
            v(i) = U0(13 + i)
            k(i) = U0(16 + i)
            v2(i) = U0(19 + i)
            k2(i) = U0(22 + i)
        END DO
        q(4) = U0(13)

            FB(1) = r(1)
            FB(2) = r(2)
            FB(3) = r(3) 
            FB(4) = q(1) - 1.0
            FB(5) = q(2)
            FB(6) = q(3)
            FB(7) = q(4)
            FB(8) = v(1)
            FB(9) = v(2)
            FB(10) = v(3) - 1.0
            FB(11) = k(1)
            FB(12) = k(2)
            FB(13) = k(3)

        DO i = 1, 3
            v1(i) = U1(i)
            k1(i) = U1(3 + i)
            r(i) = U1(6 + i)
            q(i) = U1(9 + i)
            v(i) = U1(13 + i)
            k(i) = U1(16 + i)
            v2(i) = U1(19 + i)
            k2(i) = U1(22 + i)
        END DO
        q(4) = U1(13)

		R0(1, 1) = 2.0 * (q(1) * q(1) + q(2) * q(2) - 0.5)
        R0(2, 2) = 2.0 * (q(1) * q(1) + q(3) * q(3) - 0.5)
        R0(3, 3) = 2.0 * (q(1) * q(1) + q(4) * q(4) - 0.5)
        R0(1, 2) = 2.0 * (q(2) * q(3) - q(1) * q(4))
        R0(2, 1) = 2.0 * (q(2) * q(3) + q(1) * q(4))
        R0(1, 3) = 2.0 * (q(2) * q(4) + q(1) * q(3))
        R0(3, 1) = 2.0 * (q(2) * q(4) - q(1) * q(3))
        R0(2, 3) = 2.0 * (q(3) * q(4) - q(1) * q(2))
        R0(3, 2) = 2.0 * (q(3) * q(4) + q(1) * q(2))



		

            DO i = 1, 3
				FB(13 + i) = 0.0

            	DO j = 1, 2
 					FB(13 + i) = FB(13 + i) + R0(i, j) * C(j) * v(j) - R0(i, j) * CC(j) * v2(j)
            	END DO
            	FB(13 + i) = FB(13 + i) + R0(i, 3) * C(3) * (v(3) - 1) - R0(i, 3) * CC(3) * v2(3)
            END DO
            
            DO i = 1, 3
            	FB(16 + i) = 0.0
            	r_dash(i) = 0.0
            	rv(i) = 0.0
            	DO j = 1, 3
            		rv(i) = rv(i) + R0(i, j) * v1(j) * DD(j)
            	    r_dash(i) = r_dash(i) + R0(i, j) * v(j)
            		FB(16 + i) = FB(16 + i) + R0(i, j) * D(j) * k(j) - R0(i, j) * DD(j) * k2(j)
            	END DO
            END DO
            
                        
            FB(17) = FB(17) - (rv(2) * r_dash(3) - rv(3) * r_dash(2))
            FB(18) = FB(18) - (rv(3) * r_dash(1) - rv(1) * r_dash(3)) - PAR(1)
            FB(19) = FB(19) - (rv(1) * r_dash(2) - rv(2) * r_dash(1))
            
            FB(20) = U1(10) * U1(10) + U1(11) * U1(11) + U1(12) * U1(12) + U1(13) * U1(13) - 1.0
		
		!FB(21) = v1(1) 
        !FB(22) = v1(2)
        !FB(23) = v1(3)
            
        !FB(24) = k1(1)
        !FB(25) = k1(2) 
        !FB(26) = k1(3)
		FB(20) = U1(10) * U1(10) + U1(11) * U1(11) + U1(12) * U1(12) + U1(13) * U1(13) - 1.0
        FB(21) = R0(1, 1) * CC1 * v1(1) + R0(1, 2) * CC2 * v1(2) + R0(1, 3) * CC3 * v1(3)
        FB(22) = R0(2, 1) * CC1 * v1(1) + R0(2, 2) * CC2 * v1(2) + R0(2, 3) * CC3 * v1(3)
        FB(23) = R0(3, 1) * CC1 * v1(1) + R0(3, 2) * CC2 * v1(2) + R0(3, 3) * CC3 * v1(3)
            
        FB(24) = R0(1, 1) * DD1 * k1(1) + R0(1, 2) * DD2 * k1(2) + R0(1, 3) * DD3 * k1(3)
        FB(25) = R0(2, 1) * DD1 * k1(1) + R0(2, 2) * DD2 * k1(2) + R0(2, 3) * DD3 * k1(3)
        FB(26) = R0(3, 1) * DD1 * k1(1) + R0(3, 2) * DD2 * k1(2) + R0(3, 3) * DD3 * k1(3)
        
 
    END SUBROUTINE BCND


! ------------------------------------------------------------ !

    SUBROUTINE ICND
    END SUBROUTINE ICND

! ------------------------------------------------------------ !

    SUBROUTINE FOPT
    END SUBROUTINE FOPT

! ------------------------------------------------------------ !

    !SUBROUTINE PVLS
    !END SUBROUTINE PVLS

! ------------------------------------------------------------ !

     SUBROUTINE PVLS(NDIM,U,PAR)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

      DOUBLE PRECISION, EXTERNAL :: GETP
      !out of plane displacement at s = 1
      PAR(3)=GETP('BV1',8,U)
      !out of place displacement MAX
      PAR(4)=GETP('MAX',8,U)

      END SUBROUTINE PVLS

