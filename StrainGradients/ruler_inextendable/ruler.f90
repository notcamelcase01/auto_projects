! ruler without strain gradients
SUBROUTINE FUNC(NDIM,  U, ICP, PAR, IJAC, F, DFDU, DFDP)
    IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
        DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
        DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM, NDIM), DFDP(NDIM,*)
        DOUBLE PRECISION, DIMENSION(3) :: mu, kappa, n, m, r
        DOUBLE PRECISION, DIMENSION(4) :: q
        DOUBLE PRECISION, DIMENSION(4, 3) :: A
        DOUBLE PRECISION, DIMENSION(3, 3) :: R0
        INTEGER i, j
        DOUBLE PRECISION  D1, D2, D3. DD1, DD2, DD3, l
        
		l = 0.05
		
		D1 = 10.0
		D2 = 1.0
		D3 = 1.0
		
		DD1 = 10.0 * l * l
		DD2 = l * l
		DD3 = l * l
		
        DO i = 1, 3
            n(i) = U(i)
            kappa(i) = U(3 + i)
            r(i) = U(6 + i)
            q(i) = U(9 + i)
        END DO

        q(4) = U(13)
        
        
        n(1) = C1 * mu(1)
        n(2) = C2 * mu(2)
        n(3) = C3 * (mu(3) - 1.0)
        
        m(1) = D1 * kappa(1)
        m(2) = D2 * kappa(2)
        m(3) = D3 * kappa(3) 

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

        F(1) = (n(2) * kappa(3) - n(3) * kappa(2)) / C1
        F(2) = (n(3) * kappa(1) - n(1) * kappa(3)) / C2
        F(3) = (n(1) * kappa(2) - n(2) * kappa(1)) / C3
        F(4) = (m(2) * kappa(3) - kappa(2) * m(3) + mu(3) * n(2) - n(3) * mu(2)) / D1
        F(5) = (m(3) * kappa(1) - kappa(3) * m(1) + mu(1) * n(3) - n(1) * mu(3)) / D2
        F(6) = (m(1) * kappa(2) - kappa(1) * m(2) + mu(2) * n(1) - n(2) * mu(1)) / D3

        DO i = 1, 3
            F(6 + i) = 0.0
            F(9 + i) = PAR(2) * q(i)
            DO j = 1, 3
                F(6 + i) = F(6 + i) + R0(i, j) * mu(j)
                F(9 + i) = F(9 + i) + A(i, j) * kappa(j)
            END DO
        END DO

        F(13) = PAR(2) * q(4)

        DO j = 1, 3
            F(13) = F(13) + A(4, j) * kappa(j)
        END DO
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
        DO i = 1, 13
            U(i) = 0.0
        END DO
        U(10) = 1.0
        U(9) = 1.0
        U(3) = 1.0
    END SUBROUTINE STPNT

! ------------------------------------------------------------ !

    SUBROUTINE BCND(NDIM, PAR, ICP, NBC, U0, U1, FB, IJAC, DBC)
    IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
        DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
        DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
        DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)
        DOUBLE PRECISION C1, C2, C3, D1, D2, D3
		
		C1 = 10.0 ** 10
		C2 = 10.0 ** 10
		C3 = 10.0 ** 10
		
		D1 = 10.0
		D2 = 1.0
		D3 = 1.0
        
            FB(1) = U0(7)
            FB(2) = U0(8)
            FB(3) = U0(9)
            FB(4) = U0(10) - 1.0
            FB(5) = U0(11)
            FB(6) = U0(12)
            FB(7) = U0(13)
            
            FB(8) = C1 * U1(1) - 2.0 * PAR(1) * (U1(11) * U1(12) + U1(10) * U1(13))
            FB(9) = C2 * U1(2) - 2.0 * PAR(1) * (U1(10) * U1(10) + U1(12) * U1(12) - 0.5)
            FB(10) = C3 * (U1(3) - 1.0) - 2.0 * PAR(1) * (U1(12) * U1(13) - U1(10) * U1(11))

            
            FB(11) = U1(4)
            FB(12) = U1(5)
            FB(13) = U1(6)
            FB(14) = U1(10) * U1(10) + U1(11) * U1(11) + U1(12) * U1(12) + U1(13) * U1(13) - 1.0
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
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
