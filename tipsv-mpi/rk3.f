      SUBROUTINE RK3(NEQ,FUNC,X0,XE,N,Y0,YN,N1,WORK)
**********************************************************************
*     SUBROUTINE RK NUMERICALLY INTEGRATES A SYSTEM OF NEQ           *
*     FIRST ORDER ORDINARY DIFFERENTIAL EQUATIONS OF THE FORM        *
*             DY(I)/DX = F(X, Y(1),..., Y(NEQ)),                     *
*     BY THE CLASSICAL RUNGE-KUTTA FORMULA.                          *
*                                                                    *
*     PARAMETERS                                                     *
*  === INPUT ===                                                     *
*     (1) NEQ: NUMBER OF EQUATIONS TO BE INTEGRATED                  *
*     (2) FUNC: SUBROUTINE FUNC(X,Y,F) TO EVALUATE DERIVATIVES       *
*                F(I)=DY(I)/DX                                       *
*     (3) X0: INITIAL VALUE OF INDEPENDENT VARIABLE                  *
*     (4) XE: OUTPUT POINT AT WHICH THE SOLUTION IS DESIRED          *
*     (5) N: NUMBER OF DIVISIONS                                     *
*        THE INTERVAL (X0, XE) IS DIVIDED INTO N SUBINTERVALS        *
*        WITH THE LENGTH (XE-X0)/N AND IN EACH SUBINTERVAL           *
*        THE CLASSICAL RUNGE-KUTTA FORMULA IS USED.                  *
*     (6) Y0(I) (I=1,..,NEQ): INITIAL VALUE AT X0                    *
*  === OUTPUT ===                                                    *
*     (7) YN(I) (I=1,..,NEQ): APPROXIMATE SOLUTION AT XE             *
*  === OTHER ===                                                     *
*     (8) WORK(): TWO-DIMENTIONAL ARRAY (SIZE=(NEQ,2)) TO BE         *
*                 USED INSIDE RK                                     *
*     COPYRIGHT: M. SUGIHARA, NOVEMBER 15, 1989, V. 1                *
**********************************************************************
       EXTERNAL FUNC
       INTEGER NEQ,N,I,J,N1
       REAL*8 X0,XE
       COMPLEX*16 Y0(NEQ),YN(NEQ),WORK(N1,2)
       REAL*8 H
      H = (XE - X0) / DBLE(N)
      DO 10 I = 1,N
       CALL RKSTEP(NEQ,FUNC,X0,H,Y0,YN,WORK(1,1),WORK(1,2))
       X0 = X0 + H
       DO 20 J = 1,NEQ
        Y0(J) = YN(J)
   20  CONTINUE
   10 CONTINUE
      X0 = XE
      RETURN
      END
C
      SUBROUTINE RKSTEP(NEQ,FUNC,X,H,Y0,YN,AK,W)
       REAL*8 A2,A3,B2,B3,C1,C2,C3,C4
       PARAMETER(A2 = 0.5D0, A3 = A2)
       PARAMETER(B2 = 0.5D0, B3 = B2)
       PARAMETER(C1 = 1.D0/6.D0, C2 = 1.D0/3.D0, C3 = C2, C4 = C1)
       INTEGER NEQ,I
       REAL*8 X,H
       COMPLEX*16 Y0(NEQ),YN(NEQ),AK(NEQ),W(NEQ)
       EXTERNAL FUNC
       CALL FUNC(X,Y0,AK)
      DO 10 I = 1,NEQ
       YN(I) = Y0(I) + DCMPLX( H * C1 ) * AK(I)
   10 CONTINUE
      DO 20 I = 1,NEQ
       W(I) = Y0(I) + DCMPLX( H * B2 ) * AK(I)
   20 CONTINUE
      CALL FUNC(X + A2 * H,W,AK)
      DO 30 I = 1,NEQ
       YN(I) = YN(I) + DCMPLX( H * C2 ) * AK(I)
   30 CONTINUE
      DO 40 I = 1,NEQ
       W(I) = Y0(I) + DCMPLX( H * B3 ) * AK(I)
   40 CONTINUE
      CALL FUNC(X + A3 * H,W,AK)
      DO 50 I = 1,NEQ
       YN(I) = YN(I) + DCMPLX( H * C3 ) * AK(I)
   50 CONTINUE
      DO 60 I = 1,NEQ
       W(I) = Y0(I) + DCMPLX( H ) * AK(I)
   60 CONTINUE
      CALL FUNC(X + H,W,AK)
      DO 70 I = 1,NEQ
       YN(I) = YN(I) + DCMPLX( H * C4 ) * AK(I)
   70 CONTINUE
      RETURN
      END
