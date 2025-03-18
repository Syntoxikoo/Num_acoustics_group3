	PROGRAM SCATDISC
C
C       =============================================================
C       Purpose: This program computes the complex sound pressure on
C                the generator of a infinitely thin hard disc. The
C                excitation is a plane wave perpendicular to the
C                disc surface.
C                
C       Input :  n    --- Number of computed points on each side
C                ka   --- Wavenumber times radius
C       Output:  "scatout.txt" --- file with the results in DOS text
C
C       ==============================================================
C       VC 3/99
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION EG(1200)
	EPS=1.0D-14

	WRITE(*,*)'Enter number of divisions on each side'
	READ(*,*)NPO
	WRITE(*,*)NPO
C	NPO=40
	WRITE(*,*)'Enter value of ka'
	READ(*,*)AK
	WRITE(*,*)AK
C	AK=2
C	WRITE(*,*)'Enter maximum number of terms in the series (20)'
C	READ(*,*)MAXN
C	WRITE(*,*)MAXN
	MAXN=30

	IOUT=41
	OPEN(IOUT,FILE='scatout.txt')
C	WRITE(IOUT,2000)

C	CORRECTED TO AVOID LACK OF ACCURACY (DO 120 AND IRHO LIMITS)
	DO 120 SIDE=0,1
	DO 110 IRHO=SIDE*2-1,NPO+1-SIDE
C	  CONVERT RHO (CYLINDRICAL) TO NU (OBLATE SPHEROIDAL)
	  GNU=(1-2*SIDE)*DSQRT(1-(1-DABS(DBLE(IRHO+SIDE*NPO-NPO)
     &        /DBLE(NPO)))**2)
	  SPRE=0.0D0
	  SPIM=0.0D0
	  DO 100 IN=0,MAXN
	    CALL SEGV(0,IN,AK,-1,CV,EG)
	    CALL RSWFO(0,IN,AK,0.0D0,CV,3,R1F,R1D,R2F,R2D)
	    CALL ASWFA(0,IN,AK,-1.0D0,-1,CV,S1FA,S1DA)
	    CALL NORP(IN,AK,-1,CV,PNMN)
C	    FIX EDGE POINT, MAY NOT BE NECESSARY:
	    IF (GNU.NE.0D0) THEN
	      CALL ASWFA(0,IN,AK,GNU,-1,CV,S1FB,S1DB)
	    ELSE
	      CALL ASWFABIS(IN,S1FB)
	    ENDIF
	    FACTOR=2*S1FA*S1FB/
     &             (AK*PNMN*(R1D**2+R2D**2))
	    IF (IN.EQ.2*INT(IN/2)) THEN
C	      REAL PART - EVEN N
	      SPRE=SPRE+FACTOR*R1D*(-1)**(IN/2)
C	      IMAGINARY PART - EVEN N
	      SPIM=SPIM+FACTOR*R2D*(-1)**(IN/2+1)
	    ELSE
C	      REAL PART - ODD N
	      SPRE=SPRE+FACTOR*R2D*(-1)**((IN-1)/2)
C	      IMAGINARY PART - ODD N
	      SPIM=SPIM+FACTOR*R1D*(-1)**((IN-1)/2)
	    ENDIF
	    IF (DABS(SPRETMP-SPRE).LT.DABS(SPRE)*EPS.AND.
     &          DABS(SPIMTMP-SPIM).LT.DABS(SPIM)*EPS.AND.IN.GT.0) 
     &          GOTO 105
	    SPRETMP=SPRE
	    SPIMTMP=SPIM
C	    OUTPUT PARTIAL RESULTS TO CHECK ACCURACY
C	    IF ((IRHO.LT.3.OR.IRHO.GT.(NPO-2)).AND.IN.LT.6)
C     &        WRITE(IOUT,2100) INT(IRHO+SIDE*NPO),IN,CV,S1FA,
C     &                         R1F,R2F,R1D,R2D,S1FB,PNMN
100	  CONTINUE
105	  IF (IRHO.GT.(SIDE-1).AND.IRHO.LE.NPO)
     &      WRITE(IOUT,2200) INT(IRHO+SIDE*NPO),IN-1,-SPIM,-SPRE,
     &                       DSQRT(SPRE**2+SPIM**2)
110	CONTINUE
120	CONTINUE


C2000	FORMAT(//5X,'Complex sound pressure on the generator'/
C     &         5X,'of a infinitely thin hard disc.'/
C     &         5X,'The excitation is a plane wave perpendicular'/
C     &         5X,'to the disc surface.'/
C     &         5X,'IN,CV,S1FA,R1F,R2F,R1D,R2D,S1FB,PNMN'/)
C2100	FORMAT(2I5,8D25.16)
2200	FORMAT(2I5,3D25.16)

        CLOSE(IOUT)
C
	END




	SUBROUTINE NORP(N,C,KD,CV,PNMN)
C
C       ===========================================================
C       Purpose: Compute the prolate and oblate spheroidal 
C                normalization parameters. (see Flammer)
C                Implemented only for m=0
C       Input :  n  --- Mode parameter,  n = m,m+1,...
C                c  --- Spheroidal parameter
C                KD --- Function code
C                       KD=1 for prolate;  KD=-1 for oblate
C                cv --- Characteristic value
C       Output:  PNMN-- Normalization parameter Nmn(c)
C
C       Routine called:
C                SDMN for computing expansion coefficients dk
C       ===========================================================
C       VC 3/99


	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION DF(1200)
	EPS=1.0D-14
	IP=1
	IF (N.EQ.2*INT((N)/2)) IP=0
C	N ODD, IP=1; N EVEN, IP=0
	NM=25+INT(N/2+C)
	NM2=2*NM
	CALL SDMN(0,N,C,CV,KD,DF)
	SU1=0.0D0
	DO 10 K=IP,NM2,2
	   MK=INT(K/2)+1
	   SU1=SU1+DF(MK)**2/(2*K+1)
	   IF (DABS(SW-SU1).LT.DABS(SU1)*EPS) GOTO 15
10         SW=SU1
15      PNMN=2*SU1
	RETURN
	END

	SUBROUTINE ASWFABIS(N,S1F)
C
C       ===========================================================
C       Purpose: Compute the prolate and oblate spheroidal angular
C                functions of the first kind and their derivatives
C                when m=0 and x=0
C       ===========================================================
C       VC 3/99
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	EPS=1.0D-14
	IF (N.NE.2*INT((N)/2)) THEN
	  S1F=0
	ELSE
	  R1=1.0D0
	  DO 10 J=2,N
10	    R1=R1*J
	  R2=1.0D0
	  DO 20 J=2,N/2
20	    R2=R2*J
	  S1F=(-1)**(N/2)*R1/(2**N*R2**2)
	ENDIF
	RETURN
	END


	SUBROUTINE ASWFA(M,N,C,X,KD,CV,S1F,S1D)
C
C       ===========================================================
C       Purpose: Compute the prolate and oblate spheroidal angular
C                functions of the first kind and their derivatives
C       Input :  m  --- Mode parameter,  m = 0,1,2,...
C                n  --- Mode parameter,  n = m,m+1,...
C                c  --- Spheroidal parameter
C                x  --- Argument of angular function, |x| < 1.0
C                KD --- Function code
C                       KD=1 for prolate;  KD=-1 for oblate
C                cv --- Characteristic value
C       Output:  S1F --- Angular function of the first kind
C                S1D --- Derivative of the angular function of
C                        the first kind
C       Routine called:
C                SCKB for computing expansion coefficients ck
C       ===========================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION CK(1200),DF(1200)
	EPS=1.0D-14
	X0=X
	X=DABS(X)
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	NM=25+INT((N-M)/2+C)
	NM2=NM/2-2
	CALL SDMN(M,N,C,CV,KD,DF)
	CALL SCKB(M,N,C,DF,CK)
	X1=1.0D0-X*X
	IF (M.EQ.0.AND.X1.EQ.0.0D0) THEN
	   A0=1.0D0
	ELSE
	   A0=X1**(0.5D0*M)
	ENDIF
	SU1=CK(1)
	DO 10 K=1,NM2
	   R=CK(K+1)*X1**K
	   SU1=SU1+R
	   IF (K.GE.10.AND.DABS(R/SU1).LT.EPS) GO TO 15
10         CONTINUE
15      S1F=A0*X**IP*SU1
	IF (X.EQ.1.0D0) THEN
	   IF (M.EQ.0) S1D=IP*CK(1)-2.0D0*CK(2)
	   IF (M.EQ.1) S1D=-1.0D+100
	   IF (M.EQ.2) S1D=-2.0D0*CK(1)
	   IF (M.GE.3) S1D=0.0D0
	ELSE
	   D0=IP-M/X1*X**(IP+1.0D0)
	   D1=-2.0D0*A0*X**(IP+1.0D0)
	   SU2=CK(2)
	   DO 20 K=2,NM2
	      R=K*CK(K+1)*X1**(K-1.0D0)
	      SU2=SU2+R
	      IF (K.GE.10.AND.DABS(R/SU2).LT.EPS) GO TO 25
20            CONTINUE
25         S1D=D0*A0*SU1+D1*SU2
	ENDIF
	IF (X0.LT.0.0D0.AND.IP.EQ.0) S1D=-S1D
	IF (X0.LT.0.0D0.AND.IP.EQ.1) S1F=-S1F
	X=X0
	RETURN
	END


	SUBROUTINE RSWFO(M,N,C,X,CV,KF,R1F,R1D,R2F,R2D)
C
C       ==========================================================
C       Purpose: Compute oblate radial functions of the first
C                and second kinds, and their derivatives
C       Input :  m  --- Mode parameter,  m = 0,1,2,...
C                n  --- Mode parameter,  n = m,m+1,m+2,...
C                c  --- Spheroidal parameter
C                x  --- Argument (x  0)
C                cv --- Characteristic value
C                KF --- Function code
C                       KF=1 for the first kind
C                       KF=2 for the second kind
C                       KF=3 for both the first and second kinds
C       Output:  R1F --- Radial function of the first kind
C                R1D --- Derivative of the radial function of
C                        the first kind
C                R2F --- Radial function of the second kind
C                R2D --- Derivative of the radial function of
C                        the second kind
C       Routines called:
C            (1) SDMN for computing expansion coefficients dk
C            (2) RMN1 for computing prolate or oblate radial
C                function of the first kind
C            (3) RMN2L for computing prolate or oblate radial
C                function of the second kind for a large argument
C            (4) RMN2SO for computing oblate radial functions of
C                the second kind for a small argument
C       ==========================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION DF(1200)
	KD=-1
	CALL SDMN(M,N,C,CV,KD,DF)
	IF (KF.NE.2) THEN
	   CALL RMN1(M,N,C,X,DF,KD,R1F,R1D)
	ENDIF
	IF (KF.GT.1) THEN
	   ID=10
	   IF (X.GT.1.0D-8) THEN
	      CALL RMN2L(M,N,C,X,DF,KD,R2F,R2D,ID)
	   ENDIF
	   IF (ID.GT.-1) THEN
	      CALL RMN2SO(M,N,C,X,CV,DF,KD,R2F,R2D)
	   ENDIF
	ENDIF
	RETURN
	END


	SUBROUTINE SDMN(M,N,C,CV,KD,DF)
C
C       =====================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions, dk
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                cv --- Characteristic value
C                KD --- Function code
C                       KD=1 for prolate; KD=-1 for oblate
C       Output:  DF(k) --- Expansion coefficients dk;
C                          DF(1), DF(2), ... correspond to
C                          d0, d2, ... for even n-m and d1,
C                          d3, ... for odd n-m
C       =====================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION A(1200),D(1200),G(1200),DF(1200)
	NM=25+INT(0.5*(N-M)+C)
	IF (C.LT.1.0D-10) THEN
	   DO 5 I=1,NM
5             DF(I)=0D0
	   DF((N-M)/2+1)=1.0D0
	   RETURN
	ENDIF   
	CS=C*C*KD
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	DO 10 I=1,NM+2
	   IF (IP.EQ.0) K=2*(I-1)
	   IF (IP.EQ.1) K=2*I-1
	   DK0=M+K
	   DK1=M+K+1
	   DK2=2*(M+K)
	   D2K=2*M+K
	   A(I)=(D2K+2.0)*(D2K+1.0)/((DK2+3.0)*(DK2+5.0))*CS
	   D(I)=DK0*DK1+(2.0*DK0*DK1-2.0*M*M-1.0)/((DK2-1.0)
     &          *(DK2+3.0))*CS
	   G(I)=K*(K-1.0)/((DK2-3.0)*(DK2-1.0))*CS
10      CONTINUE
	FS=1.0D0
	F1=0.0D0
	F0=1.0D-100
	KB=0
	DF(NM+1)=0.0D0
	DO 30 K=NM,1,-1
	   F=-((D(K+1)-CV)*F0+A(K+1)*F1)/G(K+1)
	   IF (DABS(F).GT.DABS(DF(K+1))) THEN
	      DF(K)=F
	      F1=F0
	      F0=F
	      IF (DABS(F).GT.1.0D+100) THEN
		 DO 12 K1=K,NM
12                  DF(K1)=DF(K1)*1.0D-100
		 F1=F1*1.0D-100
		 F0=F0*1.0D-100
	      ENDIF  
	   ELSE
	      KB=K
	      FL=DF(K+1)
	      F1=1.0D-100
	      F2=-(D(1)-CV)/A(1)*F1
	      DF(1)=F1
	      IF (KB.EQ.1) THEN
		 FS=F2
	      ELSE IF (KB.EQ.2) THEN
		 DF(2)=F2
		 FS=-((D(2)-CV)*F2+G(2)*F1)/A(2)
	      ELSE 
		 DF(2)=F2
		 DO 20 J=3,KB+1
		    F=-((D(J-1)-CV)*F2+G(J-1)*F1)/A(J-1)
		    IF (J.LE.KB) DF(J)=F
		    IF (DABS(F).GT.1.0D+100) THEN
		       DO 15 K1=1,J
15                        DF(K1)=DF(K1)*1.0D-100
		       F=F*1.0D-100
		       F2=F2*1.0D-100
		    ENDIF  
		    F1=F2
20                  F2=F
		 FS=F
	      ENDIF
	      GO TO 35
	   ENDIF
30      CONTINUE
35      SU1=0.0D0
	R1=1.0D0
	DO 40 J=M+IP+1,2*(M+IP)
40         R1=R1*J
	SU1=DF(1)*R1
	DO 45 K=2,KB
	   R1=-R1*(K+M+IP-1.5D0)/(K-1.0D0)
45           SU1=SU1+R1*DF(K)
	SU2=0.0D0
	DO 50 K=KB+1,NM
	   IF (K.NE.1) R1=-R1*(K+M+IP-1.5D0)/(K-1.0D0)
	   SU2=SU2+R1*DF(K)
	   IF (DABS(SW-SU2).LT.DABS(SU2)*1.0D-14) GOTO 55
50         SW=SU2
55      R3=1.0D0
	DO 60 J=1,(M+N+IP)/2
60         R3=R3*(J+0.5D0*(N+M+IP))
	R4=1.0D0
	DO 65 J=1,(N-M-IP)/2
65         R4=-4.0D0*R4*J
	S0=R3/(FL*(SU1/FS)+SU2)/R4
	DO 70 K=1,KB
70         DF(K)=FL/FS*S0*DF(K)
	DO 75 K=KB+1,NM
75         DF(K)=S0*DF(K)
	RETURN
	END


	SUBROUTINE RMN1(M,N,C,X,DF,KD,R1F,R1D)
C
C       =======================================================
C       Purpose: Compute prolate and oblate spheroidal radial
C                functions of the first kind for given m, n,
C                c and x
C       Routines called:
C            (1) SCKB for computing expansion coefficients c2k
C            (2) SPHJ for computing the spherical Bessel
C                functions of the first kind     
C       =======================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION CK(1200),DF(1200),SJ(0:1520),DJ(0:1520)
	EPS=1.0D-14
	IP=1
	NM1=INT((N-M)/2)
	IF (N-M.EQ.2*NM1) IP=0
	NM=25+NM1+INT(C)
	REG=1.0D0
	IF (M+NM.GT.80) REG=1.0D-200
	R0=REG
	DO 10 J=1,2*M+IP
10         R0=R0*J
	R=R0    
	SUC=R*DF(1)
	DO 15 K=2,NM
	   R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   SUC=SUC+R*DF(K)
	   IF (K.GT.NM1.AND.DABS(SUC-SW).LT.DABS(SUC)*EPS) GO TO 20
15         SW=SUC
20      CONTINUE
	IF (X.EQ.0.0) THEN
	   CALL SCKB(M,N,C,DF,CK)
	   SUM=0.0D0
	   DO 25 J=1,NM
	      SUM=SUM+CK(J)
	      IF (DABS(SUM-SW1).LT.DABS(SUM)*EPS) GO TO 30
25            SW1=SUM
30         R1=1.0D0
	   DO 35 J=1,(N+M+IP)/2
35            R1=R1*(J+0.5D0*(N+M+IP))
	   R2=1.0D0
	   DO 40 J=1,M
40            R2=2.0D0*C*R2*J
	   R3=1.0D0
	   DO 45 J=1,(N-M-IP)/2
45            R3=R3*J
	   SA0=(2.0*(M+IP)+1.0)*R1/(2.0**N*C**IP*R2*R3)
	   IF (IP.EQ.0) THEN
	      R1F=SUM/(SA0*SUC)*DF(1)*REG
	      R1D=0.0D0
	   ELSE IF (IP.EQ.1) THEN
	      R1F=0.0D0
	      R1D=SUM/(SA0*SUC)*DF(1)*REG
	   ENDIF
	   RETURN
	ENDIF
	CX=C*X
	NM2=2*NM+M
	CALL SPHJ(NM2,CX,NM2,SJ,DJ)
	A0=(1.0D0-KD/(X*X))**(0.5D0*M)/SUC  
	R1F=0.0D0
	DO 50 K=1,NM
	   L=2*K+M-N-2+IP
	   IF (L.EQ.4*INT(L/4)) LG=1
	   IF (L.NE.4*INT(L/4)) LG=-1
	   IF (K.EQ.1) THEN
	      R=R0
	   ELSE
	      R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   ENDIF
	   NP=M+2*K-2+IP
	   R1F=R1F+LG*R*DF(K)*SJ(NP)
	   IF (K.GT.NM1.AND.DABS(R1F-SW).LT.DABS(R1F)*EPS) GO TO 55
50         SW=R1F
55      R1F=R1F*A0
	B0=KD*M/X**3.0D0/(1.0-KD/(X*X))*R1F    
	SUD=0.0D0
	DO 60 K=1,NM
	   L=2*K+M-N-2+IP
	   IF (L.EQ.4*INT(L/4)) LG=1
	   IF (L.NE.4*INT(L/4)) LG=-1
	   IF (K.EQ.1) THEN
	      R=R0
	   ELSE
	      R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   ENDIF
	   NP=M+2*K-2+IP
	   SUD=SUD+LG*R*DF(K)*DJ(NP)
	   IF (K.GT.NM1.AND.DABS(SUD-SW).LT.DABS(SUD)*EPS) GO TO 65
60         SW=SUD
65      R1D=B0+A0*C*SUD
	RETURN
	END


	SUBROUTINE SCKB(M,N,C,DF,CK)
C
C       ======================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions, c2k
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                DF(k) --- Expansion coefficients dk
C       Output:  CK(k) --- Expansion coefficients ck;
C                          CK(1), CK(2), ... correspond to
C                          c0, c2, ...
C       ======================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION DF(1200),CK(1200)
	IF (C.LE.1.0D-10) C=1.0D-10
	NM=25+INT(0.5*(N-M)+C)
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	REG=1.0D0
	IF (M+NM.GT.80) REG=1.0D-200
	FAC=-0.5D0**M
	DO 35 K=0,NM-1
	   FAC=-FAC
	   I1=2*K+IP+1
	   R=REG
	   DO 10 I=I1,I1+2*M-1
10            R=R*I
	   I2=K+M+IP
	   DO 15 I=I2,I2+K-1
15            R=R*(I+0.5D0)
	   SUM=R*DF(K+1)
	   DO 20 I=K+1,NM
	      D1=2.0D0*I+IP
	      D2=2.0D0*M+D1
	      D3=I+M+IP-0.5D0
	      R=R*D2*(D2-1.0D0)*I*(D3+K)/(D1*(D1-1.0D0)*(I-K)*D3)
	      SUM=SUM+R*DF(I+1)
	      IF (DABS(SW-SUM).LT.DABS(SUM)*1.0D-14) GOTO 25
20            SW=SUM
25         R1=REG
	   DO 30 I=2,M+K
30            R1=R1*I
35         CK(K+1)=FAC*SUM/R1
	RETURN
	END 


	SUBROUTINE SPHJ(N,X,NM,SJ,DJ)
C
C       =======================================================
C       Purpose: Compute spherical Bessel functions jn(x) and
C                their derivatives
C       Input :  x --- Argument of jn(x)
C                n --- Order of jn(x)  ( n = 0,1,תתת )
C       Output:  SJ(n) --- jn(x)
C                DJ(n) --- jn'(x)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       =======================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION SJ(0:N),DJ(0:N)
	NM=N
	IF (DABS(X).EQ.1.0D-100) THEN
	   DO 10 K=0,N
	      SJ(K)=0.0D0
10            DJ(K)=0.0D0
	   SJ(0)=1.0D0
	   DJ(1)=.3333333333333333D0
	   RETURN
	ENDIF
	SJ(0)=DSIN(X)/X
	SJ(1)=(SJ(0)-DCOS(X))/X
	IF (N.GE.2) THEN
	   SA=SJ(0)
	   SB=SJ(1)
	   M=MSTA1(X,200)
	   IF (M.LT.N) THEN
	      NM=M
	   ELSE
	      M=MSTA2(X,N,15)
	   ENDIF
	   F0=0.0D0
	   F1=1.0D0-100
	   DO 15 K=M,0,-1
	      F=(2.0D0*K+3.0D0)*F1/X-F0
	      IF (K.LE.NM) SJ(K)=F
	      F0=F1
15            F1=F
	   IF (DABS(SA).GT.DABS(SB)) CS=SA/F
	   IF (DABS(SA).LE.DABS(SB)) CS=SB/F0
	   DO 20 K=0,NM
20            SJ(K)=CS*SJ(K)
	ENDIF      
	DJ(0)=(DCOS(X)-DSIN(X)/X)/X
	DO 25 K=1,NM
25         DJ(K)=SJ(K-1)-(K+1.0D0)*SJ(K)/X
	RETURN
	END


	INTEGER FUNCTION MSTA1(X,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward  
C                recurrence such that the magnitude of    
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point   
C       ===================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	A0=DABS(X)
	N0=INT(1.1*A0)+1
	F0=ENVJ(N0,A0)-MP
	N1=N0+5
	F1=ENVJ(N1,A0)-MP
	DO 10 IT=1,20             
	   NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
	   F=ENVJ(NN,A0)-MP
	   IF(ABS(NN-N1).LT.1) GO TO 20
	   N0=N1
	   F0=F1
	   N1=NN
 10        F1=F
 20     MSTA1=NN
	RETURN
	END


	INTEGER FUNCTION MSTA2(X,N,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	A0=DABS(X)
	HMP=0.5D0*MP
	EJN=ENVJ(N,A0)
	IF (EJN.LE.HMP) THEN
	   OBJ=MP
	   N0=INT(1.1*A0)
	ELSE
	   OBJ=HMP+EJN
	   N0=N
	ENDIF
	F0=ENVJ(N0,A0)-OBJ
	N1=N0+5
	F1=ENVJ(N1,A0)-OBJ
	DO 10 IT=1,20
	   NN=N1-(N1-N0)/(1.0D0-F0/F1)
	   F=ENVJ(NN,A0)-OBJ
	   IF (ABS(NN-N1).LT.1) GO TO 20
	   N0=N1
	   F0=F1
	   N1=NN
10         F1=F
20      MSTA2=NN+10
	RETURN
	END

	REAL*8 FUNCTION ENVJ(N,X)
	DOUBLE PRECISION X
	ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
	RETURN
	END


	SUBROUTINE RMN2L(M,N,C,X,DF,KD,R2F,R2D,ID)
C
C       ========================================================
C       Purpose: Compute prolate and oblate spheroidal radial
C                functions of the second kind for given m, n,
C                c and a large cx
C       Routine called:
C                SPHY for computing the spherical Bessel
C                functions of the second kind      
C       ========================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION DF(1200),SY(0:1520),DY(0:1520)
	EPS=1.0D-14
	IP=1
	NM1=INT((N-M)/2)
	IF (N-M.EQ.2*NM1) IP=0
	NM=25+NM1+INT(C)
	REG=1.0D0
	IF (M+NM.GT.80) REG=1.0D-200
	NM2=2*NM+M
	CX=C*X
	CALL SPHY(NM2,CX,NM2,SY,DY)
	R0=REG
	DO 10 J=1,2*M+IP
10         R0=R0*J
	R=R0    
	SUC=R*DF(1)
	DO 15 K=2,NM
	   R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   SUC=SUC+R*DF(K)
	   IF (K.GT.NM1.AND.DABS(SUC-SW).LT.DABS(SUC)*EPS) GO TO 20
15         SW=SUC
20      A0=(1.0D0-KD/(X*X))**(0.5D0*M)/SUC
	R2F=0.0
	DO 50 K=1,NM
	   L=2*K+M-N-2+IP
	   IF (L.EQ.4*INT(L/4)) LG=1
	   IF (L.NE.4*INT(L/4)) LG=-1
	   IF (K.EQ.1) THEN
	      R=R0
	   ELSE
	      R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   ENDIF
	   NP=M+2*K-2+IP
	   R2F=R2F+LG*R*(DF(K)*SY(NP))
	   EPS1=DABS(R2F-SW)
	   IF (K.GT.NM1.AND.EPS1.LT.DABS(R2F)*EPS) GO TO 55
50         SW=R2F
55      ID1=INT(LOG10(EPS1/DABS(R2F)+EPS))
	R2F=R2F*A0
	IF (NP.GE.NM2) THEN
	   ID=10
	   RETURN
	ENDIF
	B0=KD*M/X**3.0D0/(1.0-KD/(X*X))*R2F                
	SUD=0.0D0
	DO 60 K=1,NM
	   L=2*K+M-N-2+IP
	   IF (L.EQ.4*INT(L/4)) LG=1
	   IF (L.NE.4*INT(L/4)) LG=-1
	   IF (K.EQ.1) THEN
	      R=R0
	   ELSE
	      R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   ENDIF
	   NP=M+2*K-2+IP
	   SUD=SUD+LG*R*(DF(K)*DY(NP))
	   EPS2=DABS(SUD-SW)
	   IF (K.GT.NM1.AND.EPS2.LT.DABS(SUD)*EPS) GO TO 65
60         SW=SUD
65      R2D=B0+A0*C*SUD       
	ID2=INT(LOG10(EPS2/DABS(SUD)+EPS))
	ID=MAX(ID1,ID2)
	RETURN
	END


	SUBROUTINE SPHY(N,X,NM,SY,DY)
C
C       ======================================================
C       Purpose: Compute spherical Bessel functions yn(x) and
C                their derivatives
C       Input :  x --- Argument of yn(x) ( x ע 0 )
C                n --- Order of yn(x) ( n = 0,1,תתת )
C       Output:  SY(n) --- yn(x)
C                DY(n) --- yn'(x)
C                NM --- Highest order computed
C       ======================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION SY(0:N),DY(0:N)
	NM=N
	IF (X.LT.1.0D-60) THEN
	   DO 10 K=0,N
	      SY(K)=-1.0D+300
10            DY(K)=1.0D+300
	   RETURN
	ENDIF
	SY(0)=-DCOS(X)/X
	SY(1)=(SY(0)-DSIN(X))/X
	F0=SY(0)
	F1=SY(1)
	DO 15 K=2,N
	   F=(2.0D0*K-1.0D0)*F1/X-F0
	   SY(K)=F
	   IF (DABS(F).GE.1.0D+300) GO TO 20              
	   F0=F1
15         F1=F
20      NM=K-1
	DY(0)=(DSIN(X)+DCOS(X)/X)/X
	DO 25 K=1,NM
25         DY(K)=SY(K-1)-(K+1.0D0)*SY(K)/X
	RETURN
	END


	SUBROUTINE RMN2SO(M,N,C,X,CV,DF,KD,R2F,R2D)
C
C       =============================================================
C       Purpose: Compute oblate radial functions of the second kind
C                with a small argument, Rmn(-ic,ix) & Rmn'(-ic,ix)
C       Routines called:
C            (1) SCKB for computing the expansion coefficients c2k
C            (2) KMN for computing the joining factors
C            (3) QSTAR for computing the factor defined in (15.7.3)
C            (4) CBK for computing the the expansion coefficient
C                defined in (15.7.6)
C            (5) GMN for computing the function defined in (15.7.4)
C            (6) RMN1 for computing the radial function of the first
C                kind
C       =============================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION BK(1200),CK(1200),DF(1200),DN(1200)
	IF (DABS(DF(1)).LE.1.0D-280) THEN
	   R2F=1.0D+300
	   R2D=1.0D+300
	   RETURN
	ENDIF
	EPS=1.0D-14
	PI=3.141592653589793D0
	NM=25+INT((N-M)/2+C)
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	CALL SCKB(M,N,C,DF,CK)
	CALL KMN(M,N,C,CV,KD,DF,DN,CK1,CK2)
	CALL QSTAR(M,N,C,CK,CK1,QS,QT)
	CALL CBK(M,N,C,CV,QT,CK,BK)
	IF (X.EQ.0.0D0) THEN
	   SUM=0.0D0
	   DO 10 J=1,NM
	      SUM=SUM+CK(J)
	      IF (DABS(SUM-SW).LT.DABS(SUM)*EPS) GO TO 15
10            SW=SUM
15         IF (IP.EQ.0) THEN
	      R1F=SUM/CK1
	      R2F=-0.5D0*PI*QS*R1F
	      R2D=QS*R1F+BK(1)
	   ELSE IF (IP.EQ.1) THEN
	      R1D=SUM/CK1
	      R2F=BK(1)
	      R2D=-0.5D0*PI*QS*R1D
	   ENDIF
	   RETURN
	ELSE
	   CALL GMN(M,N,C,X,BK,GF,GD)
	   CALL RMN1(M,N,C,X,DF,KD,R1F,R1D)
	   H0=DATAN(X)-0.5D0*PI
	   R2F=QS*R1F*H0+GF
	   R2D=QS*(R1D*H0+R1F/(1.0D0+X*X))+GD
	ENDIF
	RETURN
	END


	SUBROUTINE QSTAR(M,N,C,CK,CK1,QS,QT)
C
C       =========================================================
C       Purpose: Compute Q*mn(-ic) for oblate radial functions
C                with a small argument
C       =========================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION AP(1200),CK(1200)
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	R=1.0D0/CK(1)**2
	AP(1)=R
	DO 20 I=1,M
	   S=0.0D0
	   DO 15 L=1,I
	      SK=0.0D0
	      DO 10 K=0,L
10               SK=SK+CK(K+1)*CK(L-K+1)
15            S=S+SK*AP(I-L+1)
20      AP(I+1)=-R*S
	QS0=AP(M+1)     
	DO 30 L=1,M
	   R=1.0D0
	   DO 25 K=1,L
25            R=R*(2.0D0*K+IP)*(2.0D0*K-1.0D0+IP)/(2.0D0*K)**2
30         QS0=QS0+AP(M-L+1)*R
	QS=(-1)**IP*CK1*(CK1*QS0)/C
	QT=-2.0D0/CK1*QS
	RETURN
	END


	SUBROUTINE CBK(M,N,C,CV,QT,CK,BK)
C
C       =====================================================
C       Purpose: Compute coefficient Bk's for oblate radial
C                functions with a small argument
C       =====================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION BK(1200),CK(1200),U(1200),V(1200),W(1200)
	EPS=1.0D-14
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	NM=25+INT(0.5*(N-M)+C)
	U(1)=0.0D0
	N2=NM-2
	DO 10 J=2,N2
10         U(J)=C*C
	DO 15 J=1,N2
15         V(J)=(2.0*J-1.0-IP)*(2.0*(J-M)-IP)+M*(M-1.0)-CV
	DO 20 J=1,NM-1
20         W(J)=(2.0*J-IP)*(2.0*J+1.0-IP)
	IF (IP.EQ.0) THEN
	   DO 40 K=0,N2-1
	      S1=0.0D0
	      I1=K-M+1
	      DO 30 I=I1,NM
		 IF (I.LT.0) GO TO 30
		 R1=1.0D0
		 DO 25 J=1,K
25                  R1=R1*(I+M-J)/J
		 S1=S1+CK(I+1)*(2.0*I+M)*R1
		 IF (DABS(S1-SW).LT.DABS(S1)*EPS) GO TO 35
		 SW=S1
30            CONTINUE
35            BK(K+1)=QT*S1
40         CONTINUE
	ELSE IF (IP.EQ.1) THEN
	   DO 60 K=0,N2-1
	      S1=0.0D0
	      I1=K-M+1
	      DO 50 I=I1,NM
		 IF (I.LT.0) GO TO 50
		 R1=1.0D0
		 DO 45 J=1,K
45                  R1=R1*(I+M-J)/J
		 IF (I.GT.0) S1=S1+CK(I)*(2.0*I+M-1)*R1
		 S1=S1-CK(I+1)*(2.0*I+M)*R1
		 IF (DABS(S1-SW).LT.DABS(S1)*EPS) GO TO 55
		 SW=S1
50            CONTINUE
55            BK(K+1)=QT*S1
60         CONTINUE
	ENDIF
	W(1)=W(1)/V(1)
	BK(1)=BK(1)/V(1)
	DO 65 K=2,N2
	   T=V(K)-W(K-1)*U(K)
	   W(K)=W(K)/T
65         BK(K)=(BK(K)-BK(K-1)*U(K))/T
	DO 70 K=N2-1,1,-1
70         BK(K)=BK(K)-W(K)*BK(K+1)
	RETURN
	END


	SUBROUTINE GMN(M,N,C,X,BK,GF,GD)
C
C       ===========================================================
C       Purpose: Compute gmn(-ic,ix) and its derivative for oblate
C                radial functions with a small argument
C       ===========================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION BK(1200)
	EPS=1.0D-14
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	NM=25+INT(0.5*(N-M)+C)
	XM=(1.0D0+X*X)**(-0.5D0*M)
	GF0=0.0D0
	DO 10 K=1,NM
	   GF0=GF0+BK(K)*X**(2.0*K-2.0)
	   IF (DABS((GF0-GW)/GF0).LT.EPS.AND.K.GE.10) GO TO 15
10         GW=GF0
15      GF=XM*GF0*X**(1-IP)
	GD1=-M*X/(1.0D0+X*X)*GF
	GD0=0.0D0
	DO 20 K=1,NM
	   IF (IP.EQ.0) THEN
	      GD0=GD0+(2.0D0*K-1.0)*BK(K)*X**(2.0*K-2.0)
	   ELSE
	      GD0=GD0+2.0D0*K*BK(K+1)*X**(2.0*K-1.0)
	   ENDIF
	   IF (DABS((GD0-GW)/GD0).LT.EPS.AND.K.GE.10) GO TO 25
20         GW=GD0
25      GD=GD1+XM*GD0
	RETURN
	END


	SUBROUTINE KMN(M,N,C,CV,KD,DF,DN,CK1,CK2)
C
C       ===================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions
C                and joining factors
C       ===================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION U(1200),V(1200),W(1200),DF(1200),DN(1200),
     &            TP(1200),RK(1200)
	NM=25+INT(0.5*(N-M)+C)
	NN=NM+M
	CS=C*C*KD
	IP=1
	IF (N-M.EQ.2*INT((N-M)/2)) IP=0
	DO 10 I=1,NN+3     
	   IF (IP.EQ.0) K=-2*(I-1)
	   IF (IP.EQ.1) K=-(2*I-3)
	   GK0=2.0D0*M+K
	   GK1=(M+K)*(M+K+1.0D0)
	   GK2=2.0D0*(M+K)-1.0D0
	   GK3=2.0D0*(M+K)+3.0D0
	   U(I)=GK0*(GK0-1.0D0)*CS/(GK2*(GK2+2.0D0))
	   V(I)=GK1-CV+(2.0D0*(GK1-M*M)-1.0D0)*CS/(GK2*GK3)
10         W(I)=(K+1.0D0)*(K+2.0D0)*CS/((GK2+2.0D0)*GK3)
	DO 20 K=1,M
	   T=V(M+1)
	   DO 15 L=0,M-K-1
15            T=V(M-L)-W(M-L+1)*U(M-L)/T
20         RK(K)=-U(K)/T
	R=1.0D0
	DO 25 K=1,M
	   R=R*RK(K)
25         DN(K)=DF(1)*R
	TP(NN)=V(NN+1)
	DO 30 K=NN-1,M+1,-1
	   TP(K)=V(K+1)-W(K+2)*U(K+1)/TP(K+1)
	   IF (K.GT.M+1) RK(K)=-U(K)/TP(K)
30      CONTINUE
	IF (M.EQ.0) DNP=DF(1)
	IF (M.NE.0) DNP=DN(M)
	DN(M+1)=(-1)**IP*DNP*CS/((2.0*M-1.0)*(2.0*M+1.0-4.0*IP)
     &          *TP(M+1))
	DO 35 K=M+2,NN
35         DN(K)=RK(K)*DN(K-1)
	R1=1.0D0
	DO 40 J=1,(N+M+IP)/2
40         R1=R1*(J+0.5D0*(N+M+IP))
	NM1=(N-M)/2
	R=1.0D0
	DO 45 J=1,2*M+IP
45         R=R*J
	SU0=R*DF(1)
	DO 50 K=2,NM
	   R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
	   SU0=SU0+R*DF(K)
	   IF (K.GT.NM1.AND.DABS((SU0-SW)/SU0).LT.1.0D-14) GO TO 55
50         SW=SU0
55      IF (KD.EQ.1) GOTO 70
	R2=1.0D0
	DO 60 J=1,M
60         R2=2.0D0*C*R2*J
	R3=1.0D0
	DO 65 J=1,(N-M-IP)/2
65         R3=R3*J
	SA0=(2.0*(M+IP)+1.0)*R1/(2.0**N*C**IP*R2*R3*DF(1))
	CK1=SA0*SU0
	IF (KD.EQ.-1) RETURN
70      R4=1.0D0
	DO 75 J=1,(N-M-IP)/2
75         R4=4.0D0*R4*J
	R5=1.0D0
	DO 80 J=1,M
80         R5=R5*(J+M)/C
	G0=DN(M)
	IF (M.EQ.0) G0=DF(1)
	SB0=(IP+1.0)*C**(IP+1)/(2.0*IP*(M-2.0)+1.0)/(2.0*M-1.0)
	CK2=(-1)**IP*SB0*R4*R5*G0/R1*SU0
	RETURN
	END


	SUBROUTINE SEGV(M,N,C,KD,CV,EG)
C
C       =========================================================
C       Purpose: Compute the characteristic values of spheroidal
C                wave functions
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                KD --- Function code
C                       KD=1 for Prolate; KD=-1 for Oblate
C       Output:  CV --- Characteristic value for given m, n and c
C                EG(L) --- Characteristic value for mode m and n'
C                          ( L = n' - m + 1 )
C       =========================================================
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION B(600),H(0:600),D(1800),E(1800),F(1800),CV0(600),
     &            A(1800),G(1800),EG(1200)
	IF (C.LT.1.0D-10) THEN
	   DO 5 I=1,N
5             EG(I)=(I+M)*(I+M-1.0D0)
	   GO TO 70
	ENDIF                                           
	ICM=(N-M+2)/2
	NM=10+INT(0.5*(N-M)+C)
	CS=C*C*KD
	DO 60 L=0,1
	   DO 10 I=1,NM
	      IF (L.EQ.0) K=2*(I-1)
	      IF (L.EQ.1) K=2*I-1
	      DK0=M+K
	      DK1=M+K+1.0D0
	      DK2=2.0D0*(M+K)
	      D2K=2.0D0*M+K
	      A(I)=(D2K+2.0)*(D2K+1.0)/((DK2+3.0)*(DK2+5.0))*CS
	      D(I)=DK0*DK1+(2.0*DK0*DK1-2.0*M*M-1.0)/((DK2-1.0)
     &             *(DK2+3.0))*CS
10            G(I)=K*(K-1.0)/((DK2-3.0)*(DK2-1.0))*CS
	   DO 15 K=2,NM
	      E(K)=DSQRT(A(K-1)*G(K))
15            F(K)=E(K)*E(K)
	   F(1)=0.0D0
	   E(1)=0.0D0
	   XA=D(NM)+DABS(E(NM))
	   XB=D(NM)-DABS(E(NM))
	   NM1=NM-1
	   DO 20 I=1,NM1
	      T=DABS(E(I))+DABS(E(I+1))
	      T1=D(I)+T
	      IF (XA.LT.T1) XA=T1
	      T1=D(I)-T
	      IF (T1.LT.XB) XB=T1
20         CONTINUE
	   DO 25 I=1,ICM
	      B(I)=XA
25            H(I)=XB
	   DO 55 K=1,ICM
	      DO 30 K1=K,ICM
		 IF (B(K1).LT.B(K)) THEN
		    B(K)=B(K1)
		    GO TO 35
		 ENDIF
30            CONTINUE
35            IF (K.NE.1.AND.H(K).LT.H(K-1)) H(K)=H(K-1)
40            X1=(B(K)+H(K))/2.0D0
	      CV0(K)=X1
	      IF (DABS((B(K)-H(K))/X1).LT.1.0D-14) GO TO 50
	      J=0
	      S=1.0D0
	      DO 45 I=1,NM
		 IF (S.EQ.0.0D0) S=S+1.0D-30
		 T=F(I)/S
		 S=D(I)-T-X1
		 IF (S.LT.0.0D0) J=J+1
45            CONTINUE
	      IF (J.LT.K) THEN
		 H(K)=X1
	      ELSE
		 B(K)=X1
		 IF (J.GE.ICM) THEN
		    B(ICM)=X1
		 ELSE
		    IF (H(J+1).LT.X1) H(J+1)=X1
		    IF (X1.LT.B(J)) B(J)=X1
		 ENDIF
	      ENDIF
	      GO TO 40
50            CV0(K)=X1
	      IF (L.EQ.0) EG(2*K-1)=CV0(K)
	      IF (L.EQ.1) EG(2*K)=CV0(K)
55         CONTINUE
60      CONTINUE
70      CV=EG(N-M+1)
	RETURN
	END
