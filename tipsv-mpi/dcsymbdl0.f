*
      SUBROUTINE DCSYMBDL0(A,M,N,NN,EPS,Z,W,L,LI,LJ,IER)
***********************************************************************
*  GAUSS METHOD FOR A SYMMETRIC BAND MATRIX WITH A SHORT WIDTH.       *
*  THIS ROUTINE IS FOR A VECTOR COMPUTER.                             *
*                                                                     *
*  PARAMETERS:                                                        *
*   ON ENTRY:                                                         *
*     A      THE ARRAY WITH DIMENSION (M+1)*N WHICH CONTAINS          *
*            THE LOWER BAND MATRIX IN ROW-WISE.                       *
*     M      THE HALF BAND WIDTH OF THE MATRIX A, EXCLUDING DIADONALS.*
*     N      THE ORDER OF THE MATRIX A.                               *
*     NN     THE ORDER OF WORKING ARRAYS L, LI, AND LJ.               *
*     EPS    THE TOLERANCE FOR PIVOTAL ELEMENTS.                      *
*   ON RETURN:                                                        *
*     A      THE LOWER TRIANGULAR MATRIX L AFTER DECOMPOSITION.       *
*     IER    ERROR CODE. IF IER=0, NORMAL RETURN.                     *
*   OTHERS: WORKING PARAMETERS.                                       *
*                                                                     *
*  COPYRIGHT:       FUMIKO NAGAHORI    1 SEP. 1991      VER. 1        *
***********************************************************************
        INTEGER N,M,MM,NN
        INTEGER L(NN),LI(NN),LJ(NN),IER
        REAL*8 EPS
        COMPLEX*16 A((M+1)*N),Z(M+1),W(M+1)
        INTEGER I,J,K,IJ,KK,NK,NKK,NKI
        COMPLEX*16 PIV
C
        IER = 0
        IJ = 0
        DO 20 I=1,M
          IJ = IJ + 1
          LI(IJ) = I + 1
          LJ(IJ) = 2
          L(IJ) = M * (I-1) + 1
          DO 10 J=I+1,M
            IJ = IJ + 1
            L(IJ) = L(IJ-1) + M + 1
            LI(IJ) = LI(IJ-1) + 1
            LJ(IJ) = LJ(IJ-1) + 1
   10     CONTINUE
   20   CONTINUE
        MM = (M+1) * M / 2
        DO 100 K=1,N-M
          NK = (M+1) * (K-1) + 1
          IF (cdabs(A(NK+M)) .LT. EPS) THEN
            WRITE(*,*) '(SUBR. SYMBDL) SINGULAR AT STEP = ', K
            IER = 1
            RETURN
          ENDIF
          PIV = DCMPLX(1.0D0) / A(NK+M)
          DO 30 J=2,M+1
            Z(J) = - A(NK+M*J)
            W(J) = A(NK+M*J) * PIV
            A(NK+M*J) = W(J)
   30     CONTINUE
C
          KK = NK + M + M
*VORTION VEC
          DO 50 I=1,MM
            A(KK+L(I)) = A(KK+L(I)) + W(LJ(I)) * Z(LI(I))
   50     CONTINUE
C
  100   CONTINUE
C
        DO 200 K=N-M+1,N-1
          NK = (M+1) * (K-1) + 1
          NKK = (M+1) * K - 1
          IF (cdabs(A(NK+M)) .LT. EPS) THEN
            WRITE(*,*) '(SUBR. SYMBDL) SINGULAR AT STEP = ', K
            IER = 1
            RETURN
          ENDIF
          PIV = DCMPLX(1.0D0) / A(NK+M)
          DO 150 J=2,N-K+1
            Z(J) = - A(NK+M*J)
            W(J) = A(NK+M*J) * PIV
            A(NK+M*J) = W(J)
C
            NKI = NKK + M * (J-1)
            DO 130 I=2,J
              A(NKI+I) = A(NKI+I) + W(J) * Z(I)
  130       CONTINUE
  150     CONTINUE
  200   CONTINUE
C
        RETURN
C  END OF SYMBDL
      END
*
c      SUBROUTINE DCSBDLV0(A,B,M,N,NP,EPS,Z,IER)
      SUBROUTINE DCSBDLV0(A,B,M,N,EPS,Z,IER)
***********************************************************************
*  GAUSS METHOD FOR A SYMMETRIC BAND MATRIX WITH A SHORT WIDTH.       *
*  THIS ROUTINE IS FOR A VECTOR COMPUTER.                             *
*                                                                     *
*  PARAMETERS:                                                        *
*   ON ENTRY:                                                         *
*     A      THE ARRAY WITH DIMENSION (M+1),N WHICH CONTAINS          *
*            THE LEFT HAND SIDE LOWER BAND MATRIX.                    *
*     M      THE HALF BAND WIDTH OF THE MATRIX A, EXCLUDING DIADONALS.*
*     N      THE ORDER OF THE MATRIX A.                               *
*     EPS    THE TOLERANCE FOR PIVOTAL ELEMENTS.                      *
*     B      RIGHT HAND SIDE VECTOR                                   *
*   ON RETURN:                                                        *
*     B      SOLUTION.                                                *
*     IER    ERROR CODE. IF IER=0, NORMAL RETURN.                     *
*   OTHERS: WORKING PARAMETERS.                                       *
*                                                                     *
*  COPYRIGHT:       FUMIKO NAGAHORI    1 SEP. 1991      VER. 1        *
***********************************************************************
        IMPLICIT NONE
        INTEGER M,N,NP,IER
        REAL*8 EPS
        COMPLEX*16 A(M+1,N),B(N),Z(N)
        INTEGER MM,J,K,I1,J1
        COMPLEX*16 SUM
C
C  FORWARD SUBSTITUTION
        MM = M + 1
        IF (MM .LT. 3) THEN
          Z(1) = B(1)
          DO 40 J=2,N
            Z(J) = B(J) - A(1,J) * Z(J-1)
   40     CONTINUE
          DO 55 J=1,N
            Z(J) = Z(J) / A(M+1,J)
   55     CONTINUE
          B(N) = Z(N)
          DO 60 J=1,N-1
            B(N-J) = Z(N-J) - A(1,N-J+1) * B(N-J+1)
   60     CONTINUE
        ELSE
          Z(1) = B(1)
          Z(2) = B(2) - A(M-1,2) * Z(1)
          DO 80 J=3,N
            IF (J .GT. MM) THEN
              I1 = 1
            ELSE
              I1 = MM - J + 1
            ENDIF
            SUM = DCMPLX(0.0D0)
            DO 70 K=I1,MM-1
              SUM = SUM + A(K,J) * Z(J-MM+K)
   70       CONTINUE
            Z(J) = B(J) - SUM
   80     CONTINUE
          DO 90 J=1,N
            Z(J) = Z(J) / A(M+1,J)
   90     CONTINUE
C
          B(N) = Z(N)
          B(N-1) = Z(N-1) - A(MM-1,N) * Z(N)
          DO 110 J=3,N
            J1 = N - J + 1
            I1 = 1
            IF (J .LT. MM) I1 = MM - J + 1
            SUM = DCMPLX(0.0D0)
            DO 105 K=I1,MM-1
              SUM = SUM + A(K,MM-K+J1) * B(MM-K+J1)
  105       CONTINUE
            B(J1) = Z(J1) - SUM
  110     CONTINUE
        ENDIF
C
        RETURN
      END
