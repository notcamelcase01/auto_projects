! ruler with strain gradients
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
        DOUBLE PRECISION C1, C2, C3, CC1, CC2, CC3, D1, D2, D3, DD1, DD2, DD3, l

        l = 0.0005
        
        C1 = 10.0
        C2 = 10.0
        C3 = 20.0
        
        CC1 = l * l * C1
        CC2 = l * l * C2
        CC3 = l * l * C3

        D1 = 100000.0 + l * l * C1
        D2 = 1.0 + l * l * C2
        D3 = 1.0 + 2 * l * l * C3

        DD1 = 100000.0 * l * l
        DD2 = 1.0 * l * l
        DD3 = 1.0 * l * l

        DO i = 1, 3
            v1(i) = U(i)
            k1(i) = U(3 + i)
            r(i) = U(6 + i)
            q(i) = U(9 + i)
            nbar(i) = U(13 + i)
            mbar(i) = U(16 + i)
            n(i) = U(19 + i)
            m(i) = U(22 + i)
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
        
	    v(1) = n(1) / C1
        v(2) = n(2) / C2
        v(3) = n(3) / C3 + 1

        k(1) = m(1) / D1
        k(2) = m(2) / D2
        k(3) = m(3) / D3

        F(1) = nbar(1) / CC1
        F(2) = nbar(2) / CC2 
        F(3) = nbar(3) / CC3

        F(4) = mbar(1) / DD1
        F(5) = mbar(2) / DD2
        F(6) = mbar(3) / DD3

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


    
     

        DO i = 1, 3
            neff(i) = n(i) - nbar(i)
            meff(i) = m(i) - mbar(i)
        END DO
	
        F(14) = C1 * v1(1) - (neff(2) * k(3) - neff(3) * k(2)) 
        F(15) = C2 * v1(2) - (neff(3) * k(1) - neff(1) * k(3)) 
        F(16) = C3 * v1(3) - (neff(1) * k(2) - neff(2) * k(1)) 

        F(17) = D1 * k1(1) - (meff(2) * k(3) - meff(3) * k(2) + neff(2) * v(3) - neff(3) * v(2))
        F(18) = D2 * k1(2) - (meff(3) * k(1) - meff(1) * k(3) + neff(3) * v(1) - neff(1) * v(3))
        F(19) = D3 * k1(3) - (meff(1) * k(2) - meff(2) * k(1) + neff(1) * v(2) - neff(2) * v(1))
        
        F(20) = C1 * v1(1)
    	F(21) = C2 * v1(2)
        F(22) = C3 * v1(3)
        F(23) = D1 * k1(1)
    	F(24) = D1 * k1(2)
        F(25) = D1 * k1(3)
    END SUBROUTINE FUNC

! ------------------------------------------------------------ !

    SUBROUTINE STPNT(NDIM, U, PAR, T)
    IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM), PAR(*)
        DOUBLE PRECISION, INTENT(IN) :: T
        INTEGER i
        PAR(1) = 0.0 ! load
        PAR(2) = 0.0 !
        DO i = 1, 25
            U(i) = 0.0
        END DO
        U(10) = 1
        U(9) = 1

    END SUBROUTINE STPNT

! ------------------------------------------------------------ !

    SUBROUTINE BCND(NDIM, PAR, ICP, NBC, U0, U1, FB, IJAC, DBC)
    IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
        DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
        DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
        DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
        DOUBLE PRECISION, DIMENSION(3, 3) :: R0
        DOUBLE PRECISION, DIMENSION(3) :: v, k, v1, k1, v2, k2, r,  n, m, neff, meff, nbar, mbar, C, CC, D, DD, r_dash, rv
        DOUBLE PRECISION, DIMENSION(4) :: q
        DOUBLE PRECISION, DIMENSION(4, 3) :: A

        INTEGER i, j
        DOUBLE PRECISION C1, C2, C3, CC1, CC2, CC3, D1, D2, D3, DD1, DD2, DD3, l

        l = 0.0005
        C1 = 10.0
        C2 = 10.0
        C3 = 20.0
        
        CC1 = l * l * C1
        CC2 = l * l * C2
        CC3 = l * l * C3

        D1 = 100000.0 + l * l * C1
        D2 = 1.0 + l * l * C2
        D3 = 1.0 + 2 * l * l * C3

        DD1 = 100000.0 * l * l
        DD2 = 1.0 * l * l
        DD3 = 1.0 * l * l
	
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
            nbar(i) = U0(13 + i)
            mbar(i) = U0(16 + i)
            n(i) = U0(19 + i)
            m(i) = U0(22 + i)
        END DO
        q(4) = U0(13)

            FB(1) = r(1)
            FB(2) = r(2)
            FB(3) = r(3) 
            FB(4) = q(1) - 1.0
            FB(5) = q(2)
            FB(6) = q(3)
            FB(7) = q(4)
            FB(8) = v1(1)
            FB(9) = v1(2)
            FB(10) = v1(3)
            FB(11) = k1(1)
            FB(12) = k1(2)
            FB(13) = k1(3)

        DO i = 1, 3
            v1(i) = U1(i)
            k1(i) = U1(3 + i)
            r(i) = U1(6 + i)
            q(i) = U1(9 + i)
            nbar(i) = U1(13 + i)
            mbar(i) = U1(16 + i)
            n(i) = U1(19 + i)
            m(i) = U1(22 + i)
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

        FB(14) = n(1) - 2.0 * PAR(1) * (U1(11) * U1(12) + U1(10) * U1(13))
        FB(15) = n(2) - 2.0 * PAR(1) * (U1(10) * U1(10) + U1(12) * U1(12) - 0.5)
        FB(16) = n(3) - 2.0 * PAR(1) * (U1(12) * U1(13) - U1(10) * U1(11))
        
        FB(17) = nbar(1) 
        FB(18) = nbar(2) 
        FB(19) = nbar(3) 
        FB(20) = mbar(1) 
        FB(21) = mbar(2) 
        FB(22) = mbar(3) 
        FB(23) = m(1) 
        FB(24) = m(2) 
        FB(25) = m(3) 
 	    FB(26) = U1(10) * U1(10) + U1(11) * U1(11) + U1(12) * U1(12) + U1(13) * U1(13) - 1.0
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
      PAR(3)=GETP('BV1',7,U)
      !out of place displacement MAX
      PAR(4)=GETP('MAX',7,U)

      END SUBROUTINE PVLS

