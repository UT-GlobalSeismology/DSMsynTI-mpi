*
	SUBROUTINE GLU( A, N, N1, B, EPS, WK, IP, IER )
*
*        GLU, GLUSUB
*               COPYRIGHT : H.HASEGAWA, AUG. 26 1989 V.1
*
*               SOLVES SIMULTANEOUS LINEAR EQUATIONS
*               BY GAUSSIAN ELIMINATION METHOD.
*
*        INPUT - -
*             A(N1,N)  R *8  : 2-DIM. ARRAY CONTAINING THE COEFFICIENTS.
*             N        I *4  : ORDER OF MATRIX.
*             N1       I *4  : SIZE OF ARRAY A.
*             B(N)     R *8  : 1-DIM. ARRAY CONTAINING THE RIGHT-HAND
*                              SIDE VECTOR.
*             EPS      R *8  : PARAMETER TO CHECK SINGULARITY OF THE
*                              MATRIX. ( STANDARD VALUE 3.52D-15 )
*        OUTPUT - -
*             A(N1,N)        : RESULT OF GAUSSIAN ELIMINATION.
*             B(N)           : SOLUTION.
*             IP(N)    I *4  : PIVOT NUMBER.
*             IER      I *4  : = 0,  FOR NORMAL EXECUTION.
*                              = 1,  FOR SINGULAR MATRIX.
*                              = 2,  FOR SINGULAR ORIGINAL MATRIX.
*                              = 3,  FOR INVALID ARGUEMENT.
*        WORKING  -
*             WK(N)    R *8  : 1-DIM. ARRAY.
*
	INTEGER N,N1
	INTEGER IP(N),IER
	REAL*8 EPS
	COMPLEX*16 A(N1,N),B(N),WK(N)
	INTEGER I,J,K,IPK
	COMPLEX*16 AMAX,AIK,W,T
*             LEFT-HAND SIDE
	IF( EPS.LT.0.0D0 )  EPS = 3.52D-15
	IF( ( N1.LT.N ).OR.( N.LE.0 ) )  THEN
	   IER = 3
	   WRITE(*,*) '  (SUBR. GLU)  INVALID ARGUMENT.  N1, N =', N1, N
	   RETURN
	END IF
*             CHECK ORIGINAL MATRIX.
	DO 10 I = 1, N
   10 	  WK(I) = cdabs(A(I,1))
	DO 20 J = 2, N
	  DO 30 I = 1, N
	    WK(I) = DMAX1( cdabs( WK(I) ), cdabs(A(I,J)) )
   30 	  CONTINUE
   20 	CONTINUE
	DO 40 I = 1, N
	  IF( cdabs( WK(I) ).LT.EPS )  THEN
	     IER = 2
	     WRITE(*,*) '  (SUBR. GLU)  ORIGINAL MATRIX IS SINGULAR.'
	     RETURN
	  END IF
   40 	CONTINUE
*
	IER = 0
	DO 100 K = 1, N
*             FIND MAXIMUM ELEMENT IN THE K-TH COLUMN.
	  AMAX = cdabs(A(K,K))
	  IPK = K
	  DO 110 I = K+1, N
	    AIK = cdabs(A(I,K))
	    IF( cdabs( AIK ).GT.cdabs( AMAX ) )  THEN
	       IPK = I
	       AMAX = AIK
	    END IF
  110 	  CONTINUE
	  IP(K) = IPK
*
	  IF( cdabs( AMAX ).GT.EPS )  THEN
	     IF( IPK.NE.K )  THEN
	        W = A(IPK,K)
	        A(IPK,K) = A(K,K)
	        A(K,K) = W
	     END IF
*             COMPUTE ALFA
	     DO 120 I = K+1, N
	       A(I,K) = -A(I,K)/A(K,K)
  120 	       WK(I) = A(I,K)
*
	     DO 130 J = K+1, N
	       IF( IPK.NE.K )  THEN
	          W = A(IPK,J)
	          A(IPK,J) = A(K,J)
	          A(K,J) = W
	       END IF
*             GAUSSIAN ELIMINATION
	       T = A(K,J)
	       DO 140 I = K+1, N
  140 	         A(I,J) = A(I,J) + WK(I)*T
  130 	     CONTINUE
*             MATRIX IS SINGULAR.
	  ELSE
	     IER = 1
	     IP(K) = K
	     DO 150 I = K+1, N
  150 	       A(I,K) = 0.0D0
	     WRITE(*,*) '  (SUBR. GLU)  MATRIX IS SINGULAR AT K =', K
	     RETURN
	  END IF
  100 	CONTINUE
*             RIGHT-HAND SIDE
*	ENTRY GLUSUB( A, B )
	ENTRY GLUSUB( A, N, N1, B, EPS, WK, IP, IER )
*             FORWARD ELIMINATION PROCESS
	DO 200 K = 1, N
	  IF( IP(K).NE.K ) THEN
	     W = B(IP(K))
	     B(IP(K)) = B(K)
	     B(K) = W
	  END IF
*
	  T = B(K)
	  DO 210 I = K+1, N
  210 	    B(I) = B(I) + A(I,K)*T
  200 	CONTINUE
*             BACKWARD SUBSTITUTION PROCESS
	B(N) = B(N)/A(N,N)
	DO 300 K = N-1, 1, -1
	  T = B(K+1)
	  DO 310 I = 1, K
  310 	    B(I) = B(I) - A(I,K+1)*T
	  B(K) = B(K)/A(K,K)
  300 	CONTINUE
	RETURN
	END
