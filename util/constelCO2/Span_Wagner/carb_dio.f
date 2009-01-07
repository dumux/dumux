************************************************************************
*                                                                      *
*                     RUHR-UNIVERSITAET BOCHUM                         *
*                   Fakultaet fuer Maschinenbau                        *
*                   Lehrstuhl fuer Thermodynamik                       *
*                      Universitaetsstr. 150                           *
*                        FRG - 44780 Bochum                            *
*                      Telefax: 0234/7094 163                          *
*                                                                      *
*                     Prof. Dr.-Ing. W. Wagner                         *
*                       Tel.: 0234/700 3033                            *
*                        Dr.-Ing.  R. Span                             *
*                       Tel.: 0234/700 7687                            *
*                                                                      *
************************************************************************
C
C
C
C
C
C
C
C
C*********************************************************************
      SUBROUTINE DLEQ(T,DL)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CSUB2/RG,XMOL,TC,PC,DC
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PB,DLB,DVB
      COMMON/COUT/NIN,NOUT
      IF(T .LT. TTR)THEN
      WRITE(NOUT,*)'<DLEQ> T= ',T,'<TTR. Density will not be',
     *' calculated !'
      RETURN
      END IF
      IF (T.GT.TC) THEN
      WRITE(NOUT,*)'<DLEQ> T= ',T,'>TC. Density will not be',
     *' calculated !'
      RETURN
      END IF
      CALL DLEQN(T,DL)
      RETURN
      END
C*********************************************************************
      SUBROUTINE DVEQ(T,DV)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CSUB2/RG,XMOL,TC,PC,DC
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PB,DLB,DVB
      COMMON/COUT/NIN,NOUT
      DV=0.D0
      IF(T .LT. TTR)THEN
      WRITE(NOUT,*)'<DVEQ> T= ',T,'<TTR. Density will not be',
     *' calculated !'
      RETURN
      END IF
      IF (T .GT. TC) THEN
      WRITE(NOUT,*)'<DVEQ> T= ',T,'>TC. Density will not be',
     *' calculated !'
      RETURN
      END IF
      CALL DVEQN(T,DV)
      RETURN
      END
C*********************************************************************
             SUBROUTINE VPEQN (T,PS)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/XCIS/ISL,ISR
      COMMON/XCNORM/TNORM,YNORM,ISN
      COMMON/XEQN/G(23),TPOT(20),IMAX,EQK(5)
      COMMON/EQNN/GN(23,3),TPOTN(20,3),IMAXN(3),EQKNN(5,3)
      COMMON/CNORMN/TNORMN(3),YNORMN(3),ISNN(3)
      COMMON/CISN/ISLN(3),ISRN(3)
      COMMON/XCWARN/IWARN
      COMMON/CSUB2/RG,XMOL,TC,PC,DC
      ISL=ISLN(1)
      ISR=ISRN(1)
      TNORM=TNORMN(1)
      YNORM=YNORMN(1)
      ISN=ISNN(1)
      IMAX=IMAXN(1)
      DO 10 I=1,IMAX
      G(I)=GN(I,1)
      TPOT(I)=TPOTN(I,1)
  10  CONTINUE
      DO 100 L=1,5
      EQK(L)=EQKNN(L,1)
  100 CONTINUE
      PS=YCALCN(T)
      IF(IWARN .NE. 0)  PS=0.D0
      RETURN
      END
C*********************************************************************
             SUBROUTINE DLEQN (T,ROL)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/XCIS/ISL,ISR
      COMMON/XCNORM/TNORM,YNORM,ISN
      COMMON/XEQN/G(23),TPOT(20),IMAX,EQK(5)
      COMMON/EQNN/GN(23,3),TPOTN(20,3),IMAXN(3),EQKNN(5,3)
      COMMON/CNORMN/TNORMN(3),YNORMN(3),ISNN(3)
      COMMON/CISN/ISLN(3),ISRN(3)
      COMMON/XCWARN/IWARN
      COMMON/CSUB2/RG,XMOL,TC,PC,DC
      ISL=ISLN(2)
      ISR=ISRN(2)
      TNORM=TNORMN(2)
      YNORM=YNORMN(2)
      ISN=ISNN(2)
      IMAX=IMAXN(2)
      DO 10 I=1,IMAX
      G(I)=GN(I,2)
      TPOT(I)=TPOTN(I,2)
 10   CONTINUE
      DO 100 L=1,5
      EQK(L)=EQKNN(L,2)
 100  CONTINUE
      ROL=YCALCN(T)
      IF(IWARN .NE. 0) ROL=0.D0
      RETURN
      END
C*********************************************************************
             SUBROUTINE DVEQN (T,ROV)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/XCIS/ISL,ISR
      COMMON/XCNORM/TNORM,YNORM,ISN
      COMMON/XEQN/G(23),TPOT(20),IMAX,EQK(5)
      COMMON/EQNN/GN(23,3),TPOTN(20,3),IMAXN(3),EQKNN(5,3)
      COMMON/CNORMN/TNORMN(3),YNORMN(3),ISNN(3)
      COMMON/CISN/ISLN(3),ISRN(3)
      COMMON/XCWARN/IWARN
      COMMON/CSUB2/RG,XMOL,TC,PC,DC
      ISL=ISLN(3)
      ISR=ISRN(3)
      TNORM=TNORMN(3)
      YNORM=YNORMN(3)
      ISN=ISNN(3)
      IMAX=IMAXN(3)
      DO 10 I=1,IMAX
      G(I)=GN(I,3)
      TPOT(I)=TPOTN(I,3)
 10   CONTINUE
      DO 100 L=1,5
      EQK(L)=EQKNN(L,3)
 100  CONTINUE
      ROV=YCALCN(T)
      IF(IWARN .NE. 0) ROV=0.D0
      RETURN
      END
C**********************************************************************
      SUBROUTINE FSRN(T,SR,B)
C**********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION B(20)
      COMMON/XCIS/ISL,ISR
      COMMON/XEQN/G(23),TPOT(20),IMAX,EQK(5)
      COMMON/XCNORM/TNORM,YNORM,ISN
      SR=0.D0
      TAU=FTAUN(T)
      DO 50 I=1,IMAX
50    B(I)=TAU**TPOT(I)
      GOTO(100,200,300)ISR
100   B(IMAX+1)=0.D0
      GOTO 1000
200   B(IMAX+1)=DLOG(T/TNORM)
      GOTO 1000
300   B(IMAX+1)=DLOG10(T/TNORM)
1000  CONTINUE
      KEND=IMAX+1
      DO 2000 K=1,KEND
2000  SR=SR+B(K)*G(K)
999   RETURN
      END
C*********************************************************************
      DOUBLE PRECISION FUNCTION YCALCN(T)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DUMMY(20)
      COMMON/XCIS/ISL,ISR
      COMMON/XCNORM/TNORM,YNORM,ISN
      CALL FSRN(T,SR,DUMMY)
      GOTO(10,20,30,40)ISL
10    YCALCN=DEXP(SR)*YNORM
      GOTO 999
20    YCALCN=10.D0**SR*YNORM
      GOTO 999
30    YCALCN=DEXP(SR/T*TNORM)*YNORM
      GOTO 999
40    YCALCN=(SR+1.D0)*YNORM
 999  RETURN
      END
C*********************************************************************
      DOUBLEPRECISION FUNCTION PMELT(T,IST)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/COUT/NIN,NOUT
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PB,DLB,DVB
      IST=0
      IF( T .LT. TTR ) THEN
       PMELT=1.D-10
       IST=1
       WRITE(NOUT,*)'<PMELT> T<TTR= ',T,'. No calculation',
     * ' possible !'
       RETURN
      END IF
      P=PSCHM(T,1)
      RETURN
      END
C*********************************************************************
      DOUBLE PRECISION FUNCTION FTAUN(T)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/XCNORM/TNORM,YNORM,ISN
      COMMON/XCWARN/IWARN
      IWARN=0
      GOTO(10,20,30)ISN
10    FTAUN=1.D0-T/TNORM
      IF(FTAUN .LT. 1.D-10) GOTO 998
      GOTO 999
20    FTAUN=T/TNORM
      GOTO 999
   30 FTAUN=TNORM/T-1.D0
      IF(FTAUN .LT. 1.D-10) GOTO 998
      GOTO 999
  998 FTAUN=1.D-10
      IWARN=1
999   RETURN
      END
C******************************************************************
      SUBROUTINE SATT(T,DV,DL,P)
C******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/COUT/NIN,NOUT
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PBOIL,DLB,DVB
      IF( T .GT. TC ) GOTO 500
      IF( T .LT. TTR ) GOTO 500
      EPS = 1.D-7
      CALL SATCAL(T,DV,DL,P,EPS)
      GOTO 999
 500  CONTINUE
      DL=0.D0
      DV=0.D0
      P=0.D0
      IF( T .LT. TTR )
     *   WRITE(NOUT,*)'<SATT> T=',T,'K is below the triple-point'
      IF(T .GT. TC)WRITE(NOUT,*)'<SATT> T=',T,'K is supercritical'
  999 CONTINUE
      RETURN
      END
C******************************************************************
      SUBROUTINE SATP(T,DV,DL,P)
C******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CEPS/ EPSE
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PBOIL,DLB,DVB
      COMMON/COUT/NIN,NOUT
      IF( P .GT. PC ) GOTO 500
      IF( P .LT. PTR ) GOTO 500
      EPS = 1.D-7
      CALL SATCAP(T,DV,DL,P,EPS)
      GOTO 999
 500  CONTINUE
      DL=0.D0
      DV=0.D0
      T=0.D0
      IF( P .LT. PTR )
     *   WRITE(NOUT,*)'<SATP> p=',P,'MPa is below the triple-point'
      IF(P .GT. PC)WRITE(NOUT,*)'<SATP> p=',P,'MPa is supercritical'
  999 CONTINUE
      RETURN
      END
C******************************************************************************
      DOUBLEPRECISION FUNCTION FNR(T,D)
C******************************************************************************
C
C    THIS SUBROUTINE CALCULATES F'S NORMALIZED REAL PART.
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EQU/G(74),IDPOT(60),ITPOT(60),
     *           IMAXF(0:10),IMAX,DTPOT
     *          ,ALPHA(60),BETA(60),GAMMA(60),DELTA(60),ETA(60),PHI(60)
      COMMON/CNORM/TNORM,DNORM
      COMMON/EQUFIX/TALT(60),DALT(60),DTERM(60),
     *              D2NA(60),ETNA(60), DELTANA(60)
      DIMENSION FNRA(74)
C
      FNR=0.D0
      IMAX1=0
C
      TN=TNORM/T
      TE=TN**(1.D0/DTPOT)
      DN=D/DNORM
C
      DO 30 J=0,6
       IF(IMAX1.EQ.IMAXF(J)) GOTO 30
C
       IF( J .EQ. 0) THEN
C
C            REINES POLYNOM
C
       DO 10 I=IMAX1+1,IMAXF(J)
        FNRA(I)=TE**ITPOT(I)*DN**IDPOT(I)
        FNR=FNR + G(I)*FNRA(I)
  10  CONTINUE
C
       ELSE
C
       CALL EFUNC(J,DN,EX,EXD,EXDD,EXDDD,EXD4,EXD5,0)
C
       DO 20 I=IMAX1+1,IMAXF(J)
        FNRA(I)=TE**ITPOT(I)*DN**IDPOT(I)*EX
        FNR=FNR + G(I)*FNRA(I)
  20  CONTINUE
C
       END IF
       IMAX1=IMAXF(J)
  30  CONTINUE
C
      IF( IMAXF(7) .EQ. IMAXF(6) ) GOTO 50
C
      II=0
C
      DO 40 I=IMAXF(6)+1,IMAXF(7)
C
      II=II+1
      IJ=IMAX+(II-1)*4
C
      DTN=TN-GAMMA(I)
      DDN=DN-DELTA(I)
C
      EX=DDEXP(-BETA(I)*DTN**2-ALPHA(I)*DDN**2)
C
       FNRA(I)=EX*DN**IDPOT(I)*TE**ITPOT(I)
       FHILF=FNRA(I)*G(I)
       FNRA(IJ+1)=FHILF*2.D0*BETA(I)*DTN
       FNRA(IJ+2)=FHILF*2.D0*ALPHA(I)*DDN
       FNRA(IJ+3)=FHILF*(-DTN**2)
       FNRA(IJ+4)=FHILF*(-DDN**2)
       FNR=FNR+FHILF
   40 CONTINUE
C
   50 CONTINUE
C
      IF( IMAXF(8) .EQ. IMAXF(7) ) GOTO 70
C
      TAU = (T-TNORM)/T
      DIF = (D-DNORM)/DNORM
C
      IF (DABS(TAU).LT.1.D-10) TAU=SIGN(1.D-10,TAU)
      IF (DABS(DIF).LT.1.D-10) DIF=SIGN(1.D-10,DIF)
C
      TT  = TAU*TAU
      DD  = DIF*DIF
      DN  = D/DNORM
C
      DO 60 I=(IMAXF(7)+1), IMAXF(8)
C
      II = II + 1
      IJ = IMAX + (IMAXF(7)-IMAXF(6))*4 + 3 * (II-1)
C
      B  = ALPHA(I)
      BET= BETA(I)
      GA = GAMMA(I)
      GC = DELTA(I)
      GB = IDPOT(I)/10.D0
      A  = ITPOT(I)/DTPOT
      GD = ETA(I)
C
      IF ((T.NE.TALT(I)).OR.(D.NE.DALT(I))) THEN
       DALT(I)=D
       TALT(I)=T
       DTERM(I)=(TAU+GA*DD**(.5/BET))
       DELTANA(I)=DTERM(I)**2+GB*DD**A+PHI(I)
       ETNA(I)=DDEXP(-GC*DD-GD*TT)
      END IF
C
      ETERM = ETNA(I)
      DELNA = DELTANA(I)
      DTE   = DTERM(I)
C
      FNRA(I)=(DELNA**B)*DN*ETERM
      FNR = FNR + G(I)*FNRA(I)
C
   60 CONTINUE
C
   70 CONTINUE
C
      RETURN
      END
C******************************************************************************
      DOUBLEPRECISION FUNCTION FNRD(T,D)
C******************************************************************************
C
C    THIS SUBROUTINE CALCULATES THE 1ST DN-DERIVATIVE AND THE 2ND DN,A(I)-
C    DERIVATIVE OF THE F'S NORMALIZED REAL PART.
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EQU/G(74),IDPOT(60),ITPOT(60),
     *           IMAXF(0:10),IMAX,DTPOT
     *          ,ALPHA(60),BETA(60),GAMMA(60),DELTA(60),ETA(60),PHI(60)
      COMMON/CNORM/TNORM,DNORM
      COMMON/EQUFIX/TALT(60),DALT(60),DTERM(60),
     *              D2NA(60),ETNA(60),DELTANA(60)
      DIMENSION FNRDA(74)
C
      FNRD=0.D0
      IMAX1=0
C
      TN=TNORM/T
      TE=TN**(1.D0/DTPOT)
      DN=D/DNORM
C
      DO 30 J=0,6
       IF(IMAX1.EQ.IMAXF(J)) GOTO 30
C
       IF( J .EQ. 0) THEN
C
        DO 10 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         ITPT=ITPOT(I)
         FNRDA(I)=TE**ITPT*IDPT*DN**(IDPT-1)
         FNRD=FNRD+G(I)*FNRDA(I)
  10    CONTINUE
C
       ELSE
C
        CALL EFUNC(J,DN,EX,EXD,EXDD,EXDDD,EXD4,EXD5,1)
C
        DO 20 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         ITPT=ITPOT(I)
         FNRDA(I)= TE**ITPT
     *           *(EXD*DN**IDPT + EX*IDPT*DN**(IDPT-1))
         FNRD=FNRD+G(I)*FNRDA(I)
  20    CONTINUE
C
       END IF
       IMAX1=IMAXF(J)
  30  CONTINUE
C
      IF( IMAXF(7) .EQ. IMAXF(6) ) GOTO 50
C
      II=0
C
      DO 40 I=IMAXF(6)+1,IMAXF(7)
C
      GI=G(I)
      ITPT=ITPOT(I)
      IDPT=IDPOT(I)
C
      II=II+1
      IJ=IMAX+(II-1)*4
C
      DTN=TN-GAMMA(I)
      DDN=DN-DELTA(I)
C
      EX=DDEXP(-BETA(I)*DTN**2-ALPHA(I)*DDN**2)
      EXD=-2.D0*ALPHA(I)*DDN*EX
C
       FNRDA(I)=TE**ITPT
     *          *(EXD*DN**IDPT + EX*IDPT*DN**(IDPT-1))
       FDA=FNRDA(I)
       FNRDA(IJ+1)=GI*FDA*2.D0*BETA(I)*DTN
       FNRDA(IJ+2)=2.D0*GI*ALPHA(I)*(FDA*DDN+TE**ITPT
     *             *EX*DN**IDPT)
       FNRDA(IJ+3)=GI*FDA*(-DTN**2)
       FNRDA(IJ+4)=GI*(FDA*(-DDN**2)+TE**ITPT*EX*2.D0*DN**IDPT
     *             *(-DDN))
       FNRD=FNRD+GI*FDA
   40 CONTINUE
C
   50 CONTINUE
C
C
      IF( IMAXF(8) .EQ. IMAXF(7) ) GOTO 70
C
      TAU = (T-TNORM)/T
      DIF = (D-DNORM)/DNORM
C
      IF (DABS(TAU).LT.1.D-10) TAU=SIGN(1.D-10,TAU)
      IF (DABS(DIF).LT.1.D-10) DIF=SIGN(1.D-10,DIF)
C
      TT  = TAU*TAU
      DD  = DIF*DIF
      DN  = D/DNORM
C
      DO 60 I=(IMAXF(7)+1), IMAXF(8)
C
      II = II + 1
      IJ = IMAX + (IMAXF(7)-IMAXF(6))*4 + 3 * (II-1)
C
      B  = ALPHA(I)
      BET= BETA(I)
      GA = GAMMA(I)
      GC = DELTA(I)
      GB = IDPOT(I)/10.D0
      A  = ITPOT(I)/DTPOT
      GD = ETA(I)
C
      IF ((T.NE.TALT(I)).OR.(D.NE.DALT(I))) THEN
       DALT(I)=D
       TALT(I)=T
       DTERM(I)=(TAU+GA*DD**(.5/BET))
       DELTANA(I)=DTERM(I)**2+GB*DD**A+PHI(I)
       ETNA(I)=DDEXP(-GC*DD-GD*TT)
      END IF
C
      ETERM = ETNA(I)
      DELNA = DELTANA(I)
      DTE   = DTERM(I)
C
      DB  = DELNA**B
      DBDD= (2.*GB*A*DD**(A-1.) + 4.*GA*.5/BET*DD**(.5/BET-1.)*DTE)
     *      *DIF*B*DB/DELNA
      TEX  = ETERM
      TEXDD= -2.*GC*DIF*ETERM
C
      FNRDA(I)= DB*TEX + DN*DBDD*TEX + DN*DB*TEXDD
C
      FNRD = FNRD + G(I)*FNRDA(I)
C
   60 CONTINUE
C
   70 CONTINUE
C
      RETURN
      END
C******************************************************************************
      DOUBLEPRECISION FUNCTION FNRT(T,D)
C******************************************************************************
C
C    THIS SUBROUTINE CALCULATES THE 1ST TN-DERIVATIVE AND THE 2ND TN,A(I)-
C    DERIVATIVE OF THE F'S NORMALIZED REAL PART.
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EQU/G(74),IDPOT(60),ITPOT(60),
     *           IMAXF(0:10),IMAX,DTPOT
     *          ,ALPHA(60),BETA(60),GAMMA(60),DELTA(60),ETA(60),PHI(60)
      COMMON/CNORM/TNORM,DNORM
      COMMON/EQUFIX/TALT(60),DALT(60),DTERM(60),
     *              D2NA(60),ETNA(60),DELTANA(60)
      DIMENSION FNRTA(74)
C
      FNRT=0.D0
      IMAX1=0
C
      TN=TNORM/T
      TE=TN**(1.D0/DTPOT)
      DTETN=(1.D0/DTPOT)*TN**((1.D0/DTPOT)-1.D0)
      DN=D/DNORM
C
      DO 30 J=0,6
       IF(IMAX1.EQ.IMAXF(J)) GOTO 30
C
       IF( J .EQ. 0) THEN
C
C             REINES POLYNOM
C
        DO 10 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         ITPT=ITPOT(I)
         FNRTA(I)=ITPT*TE**(ITPT-1)*DN**IDPT * DTETN
         FNRT=FNRT+G(I)*FNRTA(I)
  10    CONTINUE
C
       ELSE
C
        CALL EFUNC(J,DN,EX,EXD,EXDD,EXDDD,EXD4,EXD5,0)
C
        FAK=EX * DTETN
C
        DO 20 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         ITPT=ITPOT(I)
         FNRTA(I)=ITPT*TE**(ITPT-1)*DN**IDPT*FAK
         FNRT=FNRT+G(I)*FNRTA(I)
  20    CONTINUE
C
       END IF
       IMAX1=IMAXF(J)
  30  CONTINUE
C
      IF( IMAXF(7) .EQ. IMAXF(6) ) GOTO 50
C
      II=0
C
      DO 40 I=IMAXF(6)+1,IMAXF(7)
C
      GI=G(I)
      ITPT=ITPOT(I)
      IDPT=IDPOT(I)
C
      II=II+1
      IJ=IMAX+(II-1)*4
C
      DTN=TN-GAMMA(I)
      DDN=DN-DELTA(I)
C
      EX=DDEXP(-BETA(I)*DTN**2-ALPHA(I)*DDN**2)
      EXT=-2.D0*BETA(I)*DTN*EX
C
      EXF=EX*DTETN
C
       FNRTA(I)= DN**IDPT * ( EXF*ITPT*TE**(ITPT-1)
     *                           +EXT*    TE**ITPT    )
       FTA=FNRTA(I)
       FNRTA(IJ+1)=GI*2.D0*BETA(I)*(FTA*DTN+
     *             DN**IDPT*TE**ITPT*EX)
       FNRTA(IJ+2)=GI*2.D0*ALPHA(I)*DDN*FTA
       FNRTA(IJ+3)=GI*(FTA*(-DTN**2)-DN**IDPT*2.D0*DTN
     *             *TE**ITPT*EX)
       FNRTA(IJ+4)=GI*FTA*(-DDN**2)
       FNRT=FNRT+GI*FTA
   40 CONTINUE
C
   50 CONTINUE
C
      IF( IMAXF(8) .EQ. IMAXF(7) ) GOTO 70
C
      TAU = (T-TNORM)/T
      DIF = (D-DNORM)/DNORM
C
      IF (DABS(TAU).LT.1.D-10) TAU=SIGN(1.D-10,TAU)
      IF (DABS(DIF).LT.1.D-10) DIF=SIGN(1.D-10,DIF)
C
      TT  = TAU*TAU
      DD  = DIF*DIF
      DN  = D/DNORM
C
      DO 60 I=(IMAXF(7)+1), IMAXF(8)
C
      II = II + 1
      IJ = IMAX + (IMAXF(7)-IMAXF(6))*4 + 3 * (II-1)
C
      B  = ALPHA(I)
      BET= BETA(I)
      GA = GAMMA(I)
      GC = DELTA(I)
      GB = IDPOT(I)/10.D0
      A  = ITPOT(I)/DTPOT
      GD = ETA(I)
C
      IF ((T.NE.TALT(I)).OR.(D.NE.DALT(I))) THEN
       DALT(I)=D
       TALT(I)=T
       DTERM(I)=(TAU+GA*DD**(.5/BET))
       DELTANA(I)=DTERM(I)**2+GB*DD**A+PHI(I)
       ETNA(I)=DDEXP(-GC*DD-GD*TT)
      END IF
C
      ETERM = ETNA(I)
      DELNA = DELTANA(I)
      DTE   = DTERM(I)
C
      TEX   = ETERM
      TEXDT = 2.*GD*TAU*ETERM
      DB   = DELNA**B
      DBDT = -2.*B*DTE*DB/DELNA
C
      FNRTA(I) = DN*(DBDT*TEX + DB*TEXDT)
C
      FNRT = FNRT + G(I)*FNRTA(I)
C
   60 CONTINUE
C
   70 CONTINUE
C
      RETURN
      END
C******************************************************************************
      DOUBLEPRECISION FUNCTION FNRDD(T,D)
C******************************************************************************
C
C    THIS SUBROUTINE CALCULATES THE 2ND DN-DERIVATIVE AND THE 3RD DN,DN,A(I)-
C    DERIVATIVE OF THE F'S NORMALIZED REAL PART.
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EQU/G(74),IDPOT(60),ITPOT(60),
     *           IMAXF(0:10),IMAX,DTPOT
     *          ,ALPHA(60),BETA(60),GAMMA(60),DELTA(60),ETA(60),PHI(60)
      COMMON/CNORM/TNORM,DNORM
      COMMON/EQUFIX/TALT(60),DALT(60),DTERM(60),
     *              D2NA(60),ETNA(60),DELTANA(60)
      DIMENSION FNRDDA(74)
C
      FNRDD=0.D0
      IMAX1=0
C
      TN=TNORM/T
      TE=TN**(1.D0/DTPOT)
      DN=D/DNORM
C
      DO 30 J=0,6
       IF(IMAX1.EQ.IMAXF(J)) GOTO 30
C
       IF( J .EQ. 0) THEN
C
        DO 10 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         ITPT=ITPOT(I)
         FNRDDA(I)=TE**ITPT*IDPT*(IDPT-1)*DN**(IDPT-2)
         FNRDD=FNRDD+G(I)*FNRDDA(I)
  10    CONTINUE
C
       ELSE
C
        CALL EFUNC(J,DN,EX,EXD,EXDD,EXDDD,EXD4,EXD5,2)
C
        DO 20 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         FNRDDA(I)= TE**ITPOT(I)
     *       *(EXDD*DN**IDPT +2.D0* EXD*IDPT*DN**(IDPT-1)
     *         + EX*IDPT*(IDPT-1)*DN**(IDPT-2))
         FNRDD=FNRDD + G(I)*FNRDDA(I)
  20    CONTINUE
C
       END IF
       IMAX1=IMAXF(J)
  30  CONTINUE
C
      IF( IMAXF(7) .EQ. IMAXF(6) ) GOTO 50
C
      II=0
C
      DO 40 I=IMAXF(6)+1,IMAXF(7)
C
      GI=G(I)
      ITPT=ITPOT(I)
      IDPT=IDPOT(I)
      II=II+1
      ALI=ALPHA(I)
C
      IJ=IMAX+(II-1)*4
C
      DTN=TN-GAMMA(I)
      DDN=DN-DELTA(I)
C
      EX=DDEXP(-BETA(I)*DTN**2-ALI*DDN**2)
      EXD=-2.D0*ALI*DDN*EX
      EXDD=-2.D0*ALI*EX-2.D0*ALI*DDN*EXD
C
       FNRDDA(I)= TE**ITPT
     *       *(EXDD*DN**IDPT +2.D0* EXD*IDPT*DN**(IDPT-1)
     *         + EX*IDPT*(IDPT-1)*DN**(IDPT-2))
       FDDA=FNRDDA(I)
       FNRDDA(IJ+1)=GI*FDDA*2.D0*BETA(I)*DTN
       FNRDDA(IJ+2)=GI*2.D0*ALI*(FDDA*DDN+
     *              2.D0*TE**ITPT*(EX*IDPT*DN**(IDPT-1)+
     *              DN**IDPT*EXD))
       FNRDDA(IJ+3)=GI*FDDA*(-DTN**2)
       FNRDDA(IJ+4)=GI*(FDDA*(-DDN**2)+EX*TE**ITPT*
     *              (DN**IDPT*(8.D0*ALI*DDN**2-2)
     *              -4.D0*IDPT*DN**(IDPT-1)*DDN))
       FNRDD=FNRDD+GI*FDDA
   40 CONTINUE
C
   50 CONTINUE
C
      IF( IMAXF(8) .EQ. IMAXF(7) ) GOTO 70
C
      TAU = (T-TNORM)/T
      DIF = (D-DNORM)/DNORM
C
C KONZESSION AN DIE PROBLEME EINES FORTRAN-COMPILERS: KEIN COMPILER
C KANN 0 HOCH EINE REAL-ZAHL RECHNEN!
C
      IF (DABS(TAU).LT.1.D-10) TAU=SIGN(1.D-10,TAU)
      IF (DABS(DIF).LT.1.D-10) DIF=SIGN(1.D-10,DIF)
C
      TT  = TAU*TAU
      DD  = DIF*DIF
      DN  = D/DNORM
C
      DO 60 I=(IMAXF(7)+1), IMAXF(8)
C
      II = II + 1
      IJ = IMAX + (IMAXF(7)-IMAXF(6))*4 + 3 * (II-1)
C
      B  = ALPHA(I)
      BET= BETA(I)
      GA = GAMMA(I)
      GC = DELTA(I)
      GB = IDPOT(I)/10.D0
      A  = ITPOT(I)/DTPOT
      GD = ETA(I)
C
      IF ((T.NE.TALT(I)).OR.(D.NE.DALT(I))) THEN
       DALT(I)=D
       TALT(I)=T
       DTERM(I)=(TAU+GA*DD**(.5/BET))
       DELTANA(I)=DTERM(I)**2+GB*DD**A+PHI(I)
       ETNA(I)=DDEXP(-GC*DD-GD*TT)
      END IF
C
      ETERM = ETNA(I)
      DELNA = DELTANA(I)
      DTE   = DTERM(I)
      AT1   = DD**(.5/BET-1.)
      AT2   = DD**(A-1.)
      HT1   = 2.*GA/BET*AT1*DTE + 2.*A*GB*AT2
C
      DB    = DELNA**B
      DBDD  = DIF*B*DB/DELNA * HT1
      DBDDD = DB/DELNA*(HT1*B + HT1*HT1*B*(B-1.)/DELNA*DD
     *       + (4.*GA/BET*DIF*((.5/BET-1.)*AT1/DD*DTE
     *       + GA*.5/BET*AT1*AT1) + 4.*GB*A*(A-1.)*DIF
     *         *AT2/DD)*B*DIF)
      TEX   = ETERM
      TEXDD = -2.*GC*DIF*ETERM
      TEXDDD= (-2.*GC+4.*GC*GC*DD)*ETERM
C
      FNRDDA(I)= 2.*DBDD*TEX+2.*DB*TEXDD+DN*DBDDD*TEX+
     *           2.*DN*DBDD*TEXDD+DN*DB*TEXDDD
C
      FNRDD = FNRDD + G(I)*FNRDDA(I)
C
   60 CONTINUE
C
   70 CONTINUE
C
      RETURN
      END
C******************************************************************************
      DOUBLEPRECISION FUNCTION FNRDT(T,D)
C******************************************************************************
C
C    THIS SUBROUTINE CALCULATES THE 2ND DN,TN-DERIVATIVE AND THE 3RD DN,TN,A(I)-
C    DERIVATIVE OF THE F'S NORMALIZED REAL PART.
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EQU/G(74),IDPOT(60),ITPOT(60),
     *           IMAXF(0:10),IMAX,DTPOT
     *          ,ALPHA(60),BETA(60),GAMMA(60),DELTA(60),ETA(60),PHI(60)
      COMMON/CNORM/TNORM,DNORM
      COMMON/EQUFIX/TALT(60),DALT(60),DTERM(60),
     *              D2NA(60),ETNA(60),DELTANA(60)
      DIMENSION FNRDTA(74)
C
      FNRDT=0.D0
      IMAX1=0
C
      TN=TNORM/T
      TE=TN**(1.D0/DTPOT)
      DTETN=(1.D0/DTPOT)*TN**((1.D0/DTPOT)-1.D0)
      DN=D/DNORM
C
      DO 30 J=0,6
       IF(IMAX1.EQ.IMAXF(J)) GOTO 30
C
       IF( J .EQ. 0) THEN
C
        DO 10 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         ITPT=ITPOT(I)
         FNRDTA(I)=ITPT*TE**(ITPT-1)*IDPT*DN**(IDPT-1) * DTETN
         FNRDT=FNRDT+G(I)*FNRDTA(I)
  10    CONTINUE
C
       ELSE
C
        CALL EFUNC(J,DN,EX,EXD,EXDD,EXDDD,EXD4,EXD5,1)
C
        DO 20 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         TPT=ITPOT(I)/DTPOT
         FNRDTA(I)= TPT*TN**(TPT-1.D0)
     *           *(EXD*DN**IDPT + EX*IDPT*DN**(IDPT-1))
         FNRDT=FNRDT+G(I)*FNRDTA(I)
  20    CONTINUE
C
       END IF
       IMAX1=IMAXF(J)
  30  CONTINUE
C
      IF( IMAXF(7) .EQ. IMAXF(6) ) GOTO 50
C
      II=0
C
      DO 40 I=IMAXF(6)+1,IMAXF(7)
C
      GI=G(I)
      TPOT=ITPOT(I)/DTPOT
      IDPT=IDPOT(I)
C
      II=II+1
      IJ=IMAX+(II-1)*4
C
      DTN=TN-GAMMA(I)
      DDN=DN-DELTA(I)
C
      EX=DDEXP(-BETA(I)*DTN**2-ALPHA(I)*DDN**2)
      EXD=-2.D0*ALPHA(I)*DDN*EX
      EXT=-2.D0*BETA(I)*DTN*EX
      EXDT=4.D0*ALPHA(I)*DDN*BETA(I)*DTN*EX
C
       FNRDTA(I)= TPOT*TN**(TPOT-1.D0)
     *           *(EXD*DN**IDPT + EX*IDPT*DN**(IDPT-1))
     *                    + TN**TPOT
     *           *(EXDT*DN**IDPT + EXT*IDPT*DN**(IDPT-1))
       FDTA=FNRDTA(I)
       FNRDTA(IJ+1)=GI*2.D0*BETA(I)*(FDTA*DTN+TN**TPOT*
     *              (EXD*DN**IDPT+EX*IDPT*DN**(IDPT-1)))
       FNRDTA(IJ+2)=GI*2.D0*ALPHA(I)*(FDTA*DDN+
     *              DN**IDPT*(TPOT*TN**(TPOT-1.D0)*EX+
     *              TN**TPOT*EXT))
       FNRDTA(IJ+3)=GI*(FDTA*(-DTN**2)-TN**TPOT*2.D0*DTN*
     *              (EXD*DN**IDPT+EX*IDPT*DN**(IDPT-1)))
       FNRDTA(IJ+4)=GI*(FDTA*(-DDN**2)-2.D0*DDN*DN**IDPT*
     *             (TPOT*TN**(TPOT-1.D0)*EX+TN**TPOT*EXT))
       FNRDT=FNRDT+GI*FDTA
   40 CONTINUE
C
   50 CONTINUE
C
      IF( IMAXF(8) .EQ. IMAXF(7) ) GOTO 70
C
      TAU = (T-TNORM)/T
      DIF = (D-DNORM)/DNORM
C
      IF (DABS(TAU).LT.1.D-10) TAU=SIGN(1.D-10,TAU)
      IF (DABS(DIF).LT.1.D-10) DIF=SIGN(1.D-10,DIF)
C
      TT  = TAU*TAU
      DD  = DIF*DIF
      DN  = D/DNORM
C
      DO 60 I=(IMAXF(7)+1), IMAXF(8)
C
      II = II + 1
      IJ = IMAX + (IMAXF(7)-IMAXF(6))*4 + 3 * (II-1)
C
      B  = ALPHA(I)
      BET= BETA(I)
      GA = GAMMA(I)
      GC = DELTA(I)
      GB = IDPOT(I)/10.D0
      A  = ITPOT(I)/DTPOT
      GD = ETA(I)
C
      IF ((T.NE.TALT(I)).OR.(D.NE.DALT(I))) THEN
       DALT(I)=D
       TALT(I)=T
       DTERM(I)=(TAU+GA*DD**(.5/BET))
       DELTANA(I)=DTERM(I)**2+GB*DD**A+PHI(I)
       ETNA(I)=DDEXP(-GC*DD-GD*TT)
      END IF
C
      ETERM = ETNA(I)
      DELNA = DELTANA(I)
      DTE   = DTERM(I)
      AT1   = DD**(.5/BET-1.)
      HT1   = 2.*GA/BET*AT1*DTE + 2.*A*GB*DD**(A-1.)
C
      DB   = DELNA**B
      DBDT = -2.*B*DTE*DB/DELNA
      DBDD = DIF*B*DB/DELNA * HT1
      DBDDT= -4.*GA*.5/BET*B*DIF*DB/DELNA*AT1
     *       -2.*DTE*HT1*B*(B-1.)*DIF*DB/DELNA/DELNA
C
      TEX   = ETERM
      TEXDT = 2.*GD*TAU*ETERM
      TEXDD= -2.*GC*DIF*ETERM
      TEXDDT= (-4.*GC*GD*TAU*DIF)*ETERM
C
      FNRDTA(I) = DBDT*TEX+DN*DBDDT*TEX+DN*DBDT*TEXDD+
     *            DB*TEXDT+DN*DBDD*TEXDT+DN*DB*TEXDDT
C
      FNRDT = FNRDT + G(I)*FNRDTA(I)
C
   60 CONTINUE
C
   70 CONTINUE
C
      RETURN
      END
C****************************************************************************
      DOUBLEPRECISION FUNCTION FNRTT(T,D)
C*****************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EQU/G(74),IDPOT(60),ITPOT(60),
     *           IMAXF(0:10),IMAX,DTPOT
     *          ,ALPHA(60),BETA(60),GAMMA(60),DELTA(60),ETA(60),PHI(60)
      COMMON/CNORM/TNORM,DNORM
      COMMON/EQUFIX/TALT(60),DALT(60),DTERM(60),
     *              D2NA(60),ETNA(60),DELTANA(60)
      DIMENSION FNRTTA(74)
C
      FNRTT=0.D0
      IMAX1=0
C
      TN=TNORM/T
      DN=D/DNORM
C
      DO 30 J=0,6
       IF(IMAX1.EQ.IMAXF(J)) GOTO 30
C
       IF( J .EQ. 0) THEN
C
C                  REINES POLYNOM
C
        DO 10 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         TPT=ITPOT(I)/DTPOT
         FNRTTA(I)=TPT*(TPT-1.D0)*TN**(TPT-2.D0)*DN**IDPT
         FNRTT=FNRTT+G(I)*FNRTTA(I)
  10    CONTINUE
C
       ELSE
C
        CALL EFUNC(J,DN,EX,EXD,EXDD,EXDDD,EXD4,EXD5,0)
C
        DO 20 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         TPT=ITPOT(I)/DTPOT
         FNRTTA(I)=TPT*(TPT-1.D0)*TN**(TPT-2.D0)*DN**IDPT*EX
         FNRTT=FNRTT+G(I)*FNRTTA(I)
  20    CONTINUE
C
       END IF
C
       IMAX1=IMAXF(J)
  30  CONTINUE
C
      IF( IMAXF(7) .EQ. IMAXF(6) ) GOTO 50
C
      II=0
C
      DO 40 I=IMAXF(6)+1,IMAXF(7)
C
      GI=G(I)
      TPOT=ITPOT(I)/DTPOT
      IDPT=IDPOT(I)
      II=II+1
      BETI=BETA(I)
C
      IJ=IMAX+(II-1)*4
C
      DTN=TN-GAMMA(I)
      DDN=DN-DELTA(I)
C
      EX=DDEXP(-BETI*DTN**2-ALPHA(I)*DDN**2)
      EXT=-2.D0*BETI*DTN*EX
      EXTT=-2.D0*BETI*EX-2.D0*BETI*DTN*EXT
C
       FNRTTA(I)= DN**IDPT
     *       * (2.D0*EXT*TPOT*TN**(TPOT-1.D0) +EXTT*TN**TPOT
     *       + EX*TPOT*(TPOT-1.D0)*TN**(TPOT-2.D0) )
       FTTA=FNRTTA(I)
       FNRTTA(IJ+1)=GI*2.D0*BETI*(FTTA*DTN+2.D0*DN**IDPT*
     *              (TPOT*TN**(TPOT-1.D0)*EX+TN**TPOT*EXT))
       FNRTTA(IJ+2)=GI*2.D0*ALPHA(I)*DDN*FTTA
       FNRTTA(IJ+3)=GI*(FTTA*(-DTN**2)+EX*DN**IDPT*(TN**TPOT*
     *              (8.D0*BETI*DTN**2-2)-4.D0*DTN*
     *              TPOT*TN**(TPOT-1.D0)))
       FNRTTA(IJ+4)=GI*FTTA*(-DDN**2)
       FNRTT=FNRTT+GI*FTTA
   40 CONTINUE
C
   50 CONTINUE
C
      IF( IMAXF(8) .EQ. IMAXF(7) ) GOTO 70
C
      TAU = (T-TNORM)/T
      DIF = (D-DNORM)/DNORM
C
      IF (DABS(TAU).LT.1.D-10) TAU=SIGN(1.D-10,TAU)
      IF (DABS(DIF).LT.1.D-10) DIF=SIGN(1.D-10,DIF)
C
      TT  = TAU*TAU
      DD  = DIF*DIF
      DN  = D/DNORM
C
      DO 60 I=(IMAXF(7)+1), IMAXF(8)
C
      II = II + 1
      IJ = IMAX + (IMAXF(7)-IMAXF(6))*4 + 3 * (II-1)
C
      B  = ALPHA(I)
      BET= BETA(I)
      GA = GAMMA(I)
      GC = DELTA(I)
      GB = IDPOT(I)/10.D0
      A  = ITPOT(I)/DTPOT
      GD = ETA(I)
C
      IF ((T.NE.TALT(I)).OR.(D.NE.DALT(I))) THEN
       DALT(I)=D
       TALT(I)=T
       DTERM(I)=(TAU+GA*DD**(.5/BET))
       DELTANA(I)=DTERM(I)**2+GB*DD**A+PHI(I)
       ETNA(I)=DDEXP(-GC*DD-GD*TT)
      END IF
C
      ETERM = ETNA(I)
      DELNA = DELTANA(I)
      DTE   = DTERM(I)
C
      DB   = DELNA**B
      DBDT = -2.*B*DTE*DB/DELNA
      DBDTT= DB/DELNA*(2.*B + 4.*B*(B-1.)*DTE*DTE/DELNA)
C
      TEX   = ETERM
      TEXDT = 2.*GD*TAU*ETERM
      TEXDTT= (-2.*GD + 4.*GD*GD*TT)*ETERM
C
      FNRTTA(I) = DN*(DBDTT*TEX+2.*DBDT*TEXDT+DB*TEXDTT)
C
      FNRTT = FNRTT + G(I)*FNRTTA(I)
C
   60 CONTINUE
C
   70 CONTINUE
C
      RETURN
      END
C******************************************************************************
      DOUBLEPRECISION FUNCTION FNRDDD(T,D)
C******************************************************************************
C
C    THIS FUNCTION CALCULATES THE 3RD DN-DERIVATIVE
C    OF F'S NORMALIZED PART
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EQU/G(74),IDPOT(60),ITPOT(60),
     *           IMAXF(0:10),IMAX,DTPOT
     *          ,ALPHA(60),BETA(60),GAMMA(60),DELTA(60),ETA(60),PHI(60)
      COMMON/CNORM/TNORM,DNORM
      COMMON/EQUFIX/TALT(60),DALT(60),DTERM(60),
     *              D2NA(60),ETNA(60),DELTANA(60)
      DIMENSION FRDDDA(74)
C
      FNRDDD=0.D0
      IMAX1=0
C
      TN=TNORM/T
      TE=TN**(1.D0/DTPOT)
      DN=D/DNORM
C
      DO 30 J=0,6
       IF(IMAX1.EQ.IMAXF(J)) GOTO 30
C
       IF( J .EQ. 0) THEN
C
        DO 10 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         ITPT=ITPOT(I)
         FRDDDA(I)=TE**ITPT
     *          *IDPT*(IDPT-1)*(IDPT-2)*DN**(IDPT-3)
         FNRDDD=FNRDDD+FRDDDA(I)*G(I)
  10    CONTINUE
C
       ELSE
C
        CALL EFUNC(J,DN,EX,EXD,EXDD,EXDDD,EXD4,EXD5,3)
C
        DO 20 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         FRDDDA(I)= TE**ITPOT(I)
     *    *(      EXDDD                       *DN**IDPT
     *    +   3.D0*EXDD*IDPT                  *DN**(IDPT-1)
     *    +    3.D0*EXD*IDPT*(IDPT-1)         *DN**(IDPT-2)
     *    +          EX*IDPT*(IDPT-1)*(IDPT-2)*DN**(IDPT-3))
         FNRDDD=FNRDDD+FRDDDA(I)*G(I)
  20    CONTINUE
C
       END IF
       IMAX1=IMAXF(J)
  30  CONTINUE
C
      IF( IMAXF(7) .EQ. IMAXF(6) ) GOTO 50
C
      II=0
C
      DO 40 I=IMAXF(6)+1,IMAXF(7)
C
      GI=G(I)
      ITPT=ITPOT(I)
      IDPT=IDPOT(I)
      II=II+1
      ALI=ALPHA(I)
C
      IJ=IMAX+(II-1)*4
C
      DTN=TN-GAMMA(I)
      DDN=DN-DELTA(I)
C
      EX=DDEXP(-BETA(I)*DTN**2-ALI*DDN**2)
      EXD=-2.D0*ALI*DDN*EX
      EXDD=-2.D0*ALI*EX-2.D0*ALI*DDN*EXD
      EXDDD=-4.D0*ALI*EXD - 2.D0*ALI*DDN*EXDD
C
       FRDDDA(I)= TE**ITPT
     *    *(EXDDD*DN**IDPT +3.D0*EXDD*IDPT*DN**(IDPT-1)
     *    + 3.D0*EXD*IDPT*(IDPT-1)*DN**(IDPT-2)
     *    + EX*IDPT*(IDPT-1)*(IDPT-2)*DN**(IDPT-3))
       FDDDA=FRDDDA(I)
       FRDDDA(IJ+1)=GI*FDDDA*2.D0*BETA(I)*DTN
       FRDDDA(IJ+2)=GI*2.D0*ALI*(FDDDA*DDN+TE**ITPT*
     *               (3.D0*DN**IDPT*EXDD
     *               +6.D0*IDPT*DN**(IDPT-1)*EXD
     *               +3.D0*IDPT*(IDPT-1)*DN**(IDPT-2)*EX))
       FRDDDA(IJ+3)=GI*FDDDA*(-DTN**2)
       FRDDDA(IJ+4)=GI*(FDDDA*(-DDN**2)+TE**ITPT*(12.D0*DN**IDPT*
     *               EXD*(ALI*DDN**2-1)-6.D0*IDPT*DN**(IDPT-1)*
     *               (2.D0*DDN*EXD+EX)-6.D0*DDN*IDPT*(IDPT-1)*
     *               DN**(IDPT-2)*EX))
       FNRDDD=FNRDDD+GI*FDDDA
   40 CONTINUE
C
   50 CONTINUE
C
      IF( IMAXF(8) .EQ. IMAXF(7) ) GOTO 70
C
      TAU = (T-TNORM)/T
      DIF = (D-DNORM)/DNORM
C
      IF (DABS(TAU).LT.1.D-10) TAU=SIGN(1.D-10,TAU)
      IF (DABS(DIF).LT.1.D-10) DIF=SIGN(1.D-10,DIF)
C
      TT  = TAU*TAU
      DD  = DIF*DIF
C
      DO 60 I=IMAXF(7)+1,IMAXF(8)
C
      II = II + 1
      IJ = IMAX + (IMAXF(7)-IMAXF(6))*4 + 3 * (II-1)
C
      B  = ALPHA(I)
      BET= BETA(I)
      GA = GAMMA(I)
      GC = DELTA(I)
      GB = IDPOT(I)/10.D0
      A  = ITPOT(I)/DTPOT
      GD = ETA(I)
C
      IF ((T.NE.TALT(I)).OR.(D.NE.DALT(I))) THEN
       DALT(I)=D
       TALT(I)=T
       DTERM(I)=(TAU+GA*DD**(.5/BET))
       DELTANA(I)=DTERM(I)**2+GB*DD**A+PHI(I)
       ETNA(I)=DDEXP(-GC*DD-GD*TT)
      END IF
C
      ETERM = ETNA(I)
      DELNA = DELTANA(I)
      DTE   = DTERM(I)
      AT1   = DD**(.5/BET-1.)
      AT2   = DD**(A-1.)
      HT1   = 2.*GA/BET*AT1*DTE + 2.*A*GB*AT2
      HT2   = (.5/BET-1.)*AT1/DD*DTE + GA*.5/BET*AT1*AT1
      HT3   = 4.*GA/BET*DIF*HT2 + 4.*GB*A*(A-1.)*DIF*AT2/DD
      HT4   = 4.*GA/BET*(HT2 + 2.*(.5/BET-1.)*((.5/BET-2.)*AT1/DD*DTE
     *        + 3.*GA*.5/BET*AT1*AT1))
     *        + 4.*GB*A*(A-1.)*(2.*A-3.)*AT2/DD
C
      DB    = DELNA**B
      DBDD  = DIF*B*DB/DELNA * HT1
      DBDDD = DB/DELNA*(HT1*B + HT1*HT1*B*(B-1.)/DELNA*DD
     *       + HT3*B*DIF)
      DBDDDD= DB/DELNA*(B*(2.*HT3+DIF*HT4)
     *       +B*(B-1.)/DELNA*(DIF*HT1*HT1+3.*DD*HT1*HT3+2.*DIF
     *                                *HT1*HT1)
     *       +B*(B-1.)*(B-2.)/DELNA/DELNA*HT1*HT1*HT1*DIF*DD)
C
      TEX   = ETERM
      TEXDD = -2.*GC*DIF*ETERM
      TEXDDD= (-2.*GC+4.*GC*GC*DD)*ETERM
      TEXDD3= (12.*GC*GC*DIF-8.*GC*GC*GC*DIF*DD)*ETERM
C
      FRDDDA(I)=DB*(3.*TEXDDD+DN*TEXDD3) + DBDD*(6.*TEXDD+3.*DN*TEXDDD)
     *          +DBDDD*(3.*TEX+3.*DN*TEXDD) + DBDDDD*DN*TEX
C
      FNRDDD = FNRDDD + G(I)*FRDDDA(I)
C
   60 CONTINUE
C
   70 CONTINUE
C
      RETURN
      END
C******************************************************************************
      DOUBLEPRECISION FUNCTION FNRDDT(T,D)
C******************************************************************************
C
C    THIS FUNCTION CALCULATES THE 3RD DN,DN,TN-DERIVATIVE
C    OF F'S NORMALIZED PART
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EQU/G(74),IDPOT(60),ITPOT(60),
     *           IMAXF(0:10),IMAX,DTPOT
     *          ,ALPHA(60),BETA(60),GAMMA(60),DELTA(60),ETA(60),PHI(60)
      COMMON/EQUFIX/TALT(60),DALT(60),DTERM(60),
     *              D2NA(60),ETNA(60),DELTANA(60)
      COMMON/CNORM/TNORM,DNORM
C
      FNRDDT=0.D0
      IMAX1=0
C
      TN=TNORM/T
      TE=TN**(1.D0/DTPOT)
      DTETN=(1.D0/DTPOT)*TN**((1.D0/DTPOT)-1.D0)
      DN=D/DNORM
C
      DO 30 J=0,6
       IF(IMAX1.EQ.IMAXF(J)) GOTO 30
C
       IF( J .EQ. 0) THEN
C
        DO 10 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         ITPT=ITPOT(I)
         FNRDDT=FNRDDT+G(I)*ITPT*TE**(ITPT-1)
     *                     *IDPT*(IDPT-1)*DN**(IDPT-2)
  10    CONTINUE
C
       ELSE
C
        CALL EFUNC(J,DN,EX,EXD,EXDD,EXDDD,EXD4,EXD5,2)
C
        DO 20 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         ITPT=ITPOT(I)
         FNRDDT=FNRDDT+G(I)* ITPT*TE**(ITPT-1)
     *       *(EXDD*DN**IDPT +2.D0* EXD*IDPT*DN**(IDPT-1)
     *         + EX*IDPT*(IDPT-1)*DN**(IDPT-2))
  20    CONTINUE
C
       END IF
       IMAX1=IMAXF(J)
  30  CONTINUE
C
      FNRDDT=FNRDDT * DTETN
C
C
      IF( IMAXF(7) .EQ. IMAXF(6) ) GOTO 50
C
      II=0
C
      DO 40 I=IMAXF(6)+1,IMAXF(7)
C
      II=II+1
C
      TPOT=ITPOT(I)/DTPOT
      IDPT=IDPOT(I)
      BETI=BETA(I)
      ALI=ALPHA(I)
C
      DTN=TN-GAMMA(I)
      DDN=DN-DELTA(I)
C
C
      EX=DDEXP(-BETI*DTN**2-ALI*DDN**2)
      EXD=-2.D0*ALI*DDN*EX
      EXDD=-2.D0*ALI*EX-2.D0*ALI*DDN*EXD
      EXT=-2.D0*BETI*DTN*EX
      EXDT=4.D0*ALI*DDN*BETI*DTN*EX
      EXDDT=(4.D0*ALI*BETI*DTN - 8.D0*ALI**2*BETI*DDN**2*DTN)*EX
C
       FNRDDT=FNRDDT+G(I)* (TPOT*TN**(TPOT-1.D0)
     *       *(EXDD*DN**IDPT +2.D0* EXD*IDPT*DN**(IDPT-1)
     *         + EX*IDPT*(IDPT-1)*DN**(IDPT-2))
     *                      +TN**TPOT
     *       *(EXDDT*DN**IDPT+2.D0*EXDT*IDPT*DN**(IDPT-1)
     *         + EXT*IDPT*(IDPT-1)*DN**(IDPT-2)))
   40 CONTINUE
C
   50 CONTINUE
C
      IF( IMAXF(8) .EQ. IMAXF(7) ) GOTO 70
C
      TAU = (T-TNORM)/T
      DIF = (D-DNORM)/DNORM
C
      IF (DABS(TAU).LT.1.D-10) TAU=SIGN(1.D-10,TAU)
      IF (DABS(DIF).LT.1.D-10) DIF=SIGN(1.D-10,DIF)
C
      TT  = TAU*TAU
      DD  = DIF*DIF
C
      DO 60 I=IMAXF(7)+1,IMAXF(8)
C
      II = II + 1
      IJ = IMAX + (IMAXF(7)-IMAXF(6))*4 + 3 * (II-1)
C
      B  = ALPHA(I)
      BET= BETA(I)
      GA = GAMMA(I)
      GC = DELTA(I)
      GB = IDPOT(I)/10.D0
      A  = ITPOT(I)/DTPOT
      GD = ETA(I)
C
      IF ((T.NE.TALT(I)).OR.(D.NE.DALT(I))) THEN
       DALT(I)=D
       TALT(I)=T
       DTERM(I)=(TAU+GA*DD**(.5/BET))
       DELTANA(I)=DTERM(I)**2+GB*DD**A+PHI(I)
       ETNA(I)=DDEXP(-GC*DD-GD*TT)
      END IF
C
      ETERM = ETNA(I)
      DELNA = DELTANA(I)
      DTE   = DTERM(I)
      AT1   = DD**(.5/BET-1.)
      AT2   = DD**(A-1.)
      HT1   = 2.*GA/BET*AT1*DTE + 2.*A*GB*AT2
      HT2   = (.5/BET-1.)*AT1/DD*DTE + GA*.5/BET*AT1*AT1
      HT3   = 4.*GA/BET*DIF*HT2 + 4.*GB*A*(A-1.)*DIF*AT2/DD
C
      DB    = DELNA**B
      AT3   = DB/DELNA
      DBDD  = DIF*B*AT3 * HT1
      DBDDD = AT3*(HT1*B + HT1*HT1*B*(B-1.)/DELNA*DD + HT3*B*DIF)
      DBDT  = -2.*B*DTE*AT3
      DBDDT = -4.*GA*.5/BET*B*DIF*AT3*AT1
     *       -2.*DTE*HT1*B*(B-1.)*DIF*AT3/DELNA
      DBDDDT= -2.*GA*B/BET*(1./BET-1.)*AT1*AT3
     *        -2.*B*(B-1.)*(AT3/DELNA*(2.*GA/BET*AT1*DD*HT1
     *                                    +HT3*DIF*DTE)
     *                     +HT1*DTE*AT3/DELNA*(1.+(B-2.)*DD*HT1/DELNA) )
      TEX    = ETERM
      TEXDT  = 2.*GD*TAU*ETERM
      TEXDD  = -2.*GC*DIF*ETERM
      TEXDDD = (-2.*GC+4.*GC*GC*DD)*ETERM
      TEXDDT = (-4.*GC*GD*TAU*DIF)*ETERM
      TEXDDDT= -4.*GC*GD*TAU*(1.-2.*GC*DD)*ETERM
C
      FNRDDTA= DB*(2.*TEXDDT+DN*TEXDDDT) + DBDD*(2.*TEXDT+2.*DN*TEXDDT)
     *        +DBDDT*(2.*TEX+2.*DN*TEXDD)+ DBDT*(2.*TEXDD+DN*TEXDDD)
     *        +DBDDD*DN*TEXDT            + DBDDDT*DN*TEX
C
      FNRDDT = FNRDDT + FNRDDTA * G(I)
C
   60 CONTINUE
C
   70 CONTINUE
C
      RETURN
      END
C******************************************************************************
      DOUBLEPRECISION FUNCTION FNRDTT(T,D)
C******************************************************************************
C
C    THIS FUNCTION CALCULATES THE 3RD DN,TN,TN-DERIVATIVE
C    OF F'S NORMALIZED PART
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EQU/G(74),IDPOT(60),ITPOT(60),
     *           IMAXF(0:10),IMAX,DTPOT
     *          ,ALPHA(60),BETA(60),GAMMA(60),DELTA(60),ETA(60),PHI(60)
      COMMON/CNORM/TNORM,DNORM
      COMMON/EQUFIX/TALT(60),DALT(60),DTERM(60),
     *              D2NA(60),ETNA(60),DELTANA(60)
C
      FNRDTT=0.D0
      IMAX1=0
C
      TN=TNORM/T
      DN=D/DNORM
C
      DO 30 J=0,6
       IF(IMAX1.EQ.IMAXF(J)) GOTO 30
C
       IF( J .EQ. 0) THEN
C
        DO 10 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         TPT=ITPOT(I)/DTPOT
         FNRDTT=FNRDTT+G(I)*TPT*(TPT-1.D0)*TN**(TPT-2.D0)
     *                     *IDPT*DN**(IDPT-1)
  10    CONTINUE
C
       ELSE
C
        CALL EFUNC(J,DN,EX,EXD,EXDD,EXDDD,EXD4,EXD5,1)
C
        DO 20 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         TPT=ITPOT(I)/DTPOT
       FNRDTT=FNRDTT+G(I)* TPT*(TPT-1.D0)*TN**(TPT-2.D0)
     *           *(EXD*DN**IDPT + EX*IDPT*DN**(IDPT-1))
  20    CONTINUE
C
       END IF
       IMAX1=IMAXF(J)
  30  CONTINUE
C
      IF( IMAXF(7) .EQ. IMAXF(6) ) GOTO 50
C
      II=0
C
      DO 40 I=IMAXF(6)+1,IMAXF(7)
C
      II=II+1
C
      TPOT=ITPOT(I)/DTPOT
      IDPT=IDPOT(I)
      BETI=BETA(I)
      ALI=ALPHA(I)
C
      DTN=TN-GAMMA(I)
      DDN=DN-DELTA(I)
C
      EX=DDEXP(-BETI*DTN**2-ALI*DDN**2)
      EXD=-2.D0*ALI*DDN*EX
      EXDD=-2.D0*ALI*EX-2.D0*ALI*DDN*EXD
      EXT=-2.D0*BETI*DTN*EX
      EXTT=-2.D0*BETI*EX-2.D0*BETI*DTN*EXT
      EXDT=4.D0*ALI*DDN*BETI*DTN*EX
      EXDTT=-2.D0*BETI*( 1.D0 +DTN )* EXD
C
       FNRDTT= FNRDTT + G(I)* (IDPT*DN**(IDPT-1)
     *       * (2.D0*EXT*TPOT*TN**(TPOT-1.D0) +EXTT*TN**TPOT
     *       + EX*TPOT*(TPOT-1.D0)*TN**(TPOT-2.D0) )
     *                         +DN**IDPT
     *       * (2.D0*EXDT*TPOT*TN**(TPOT-1.D0) +EXDTT*TN**TPOT
     *       + EXD*TPOT*(TPOT-1.D0)*TN**(TPOT-2.D0) ))
   40 CONTINUE
C
   50 CONTINUE
C
      IF( IMAXF(8) .EQ. IMAXF(7) ) GOTO 70
C
      TAU = (T-TNORM)/T
      DIF = (D-DNORM)/DNORM
C
      IF (DABS(TAU).LT.1.D-10) TAU=SIGN(1.D-10,TAU)
      IF (DABS(DIF).LT.1.D-10) DIF=SIGN(1.D-10,DIF)
C
      TT  = TAU*TAU
      DD  = DIF*DIF
C
      DO 60 I=IMAXF(7)+1,IMAXF(8)
C
      II = II + 1
      IJ = IMAX + (IMAXF(7)-IMAXF(6))*4 + 3 * (II-1)
C
      B  = ALPHA(I)
      BET= BETA(I)
      GA = GAMMA(I)
      GC = DELTA(I)
      GB = IDPOT(I)/10.D0
      A  = ITPOT(I)/DTPOT
      GD = ETA(I)
C
      IF ((T.NE.TALT(I)).OR.(D.NE.DALT(I))) THEN
       DALT(I)=D
       TALT(I)=T
       DTERM(I)=(TAU+GA*DD**(.5/BET))
       DELTANA(I)=DTERM(I)**2+GB*DD**A+PHI(I)
       ETNA(I)=DDEXP(-GC*DD-GD*TT)
      END IF
C
      ETERM = ETNA(I)
      DELNA = DELTANA(I)
      DTE   = DTERM(I)
      AT1   = DD**(.5/BET-1.)
      AT2   = DD**(A-1.)
      HT1   = 2.*GA/BET*AT1*DTE + 2.*A*GB*DD**(A-1.)
C
      DB    = DELNA**B
      AT3   = DB/DELNA
      DBDD  = DIF*B*AT3 * HT1
      DBDT  = -2.*B*DTE*AT3
      DBDTT = AT3*(2.*B + 4.*B*(B-1.)*DTE*DTE/DELNA)
      DBDDT = AT3*(-4.*GA*.5/BET*B*DIF*AT1
     *             -2.*DTE*HT1*B*(B-1.)*DIF/DELNA)
      DBDDTT= 2.*B*(B-1.)*DIF*AT3/DELNA*(HT1+4.*GA/BET*AT1*DTE
     *                                   +2.*(B-2.)*DTE*DTE*HT1/DELNA)
C
      TEX    = ETERM
      TEXDT  = 2.*GD*TAU*ETERM
      TEXDD  = -2.*GC*DIF*ETERM
      TEXDTT = (-2.*GD + 4.*GD*GD*TT)*ETERM
      TEXDDT = (-4.*GC*GD*TAU*DIF)*ETERM
      TEXDDTT= 4.*GC*GD*DIF*(1.-2.*GD*TT)*ETERM
C
      FNRDTTA= 2.*DN*DBDDT*TEXDT + DBDTT*(TEX+DN*TEXDD)
     *        +DBDT*(2.*TEXDT+2.*DN*TEXDDT) + DN*DBDDTT*TEX
     *        +DN*DBDD*TEXDTT + DB*(TEXDTT+DN*TEXDDTT)
C
      FNRDTT = FNRDTT + FNRDTTA * G(I)
C
   60 CONTINUE
C
   70 CONTINUE
C
      RETURN
      END
C******************************************************************************
      DOUBLEPRECISION FUNCTION FNRTTT(T,D)
C******************************************************************************
C
C    THIS FUNCTION CALCULATES THE 3RD TN-DERIVATIVE
C    OF F'S NORMALIZED PART
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EQU/G(74),IDPOT(60),ITPOT(60),
     *           IMAXF(0:10),IMAX,DTPOT
     *          ,ALPHA(60),BETA(60),GAMMA(60),DELTA(60),ETA(60),PHI(60)
      COMMON/CNORM/TNORM,DNORM
      COMMON/EQUFIX/TALT(60),DALT(60),DTERM(60),
     *              D2NA(60),ETNA(60),DELTANA(60)
C
      FNRTTT=0.D0
      IMAX1=0
C
      TN=TNORM/T
      DN=D/DNORM
C
      DO 30 J=0,6
       IF(IMAX1.EQ.IMAXF(J)) GOTO 30
C
       IF( J .EQ. 0) THEN
C
C               REINES POLYNOM
C
        DO 10 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         TPT=ITPOT(I)/DTPOT
         FNRTTT=FNRTTT+G(I)*TPT*(TPT-1.D0)*(TPT-2.D0)
     *                     *TN**(TPT-3.D0)*DN**IDPT
  10    CONTINUE
C
       ELSE
C
        CALL EFUNC(J,DN,EX,EXD,EXDD,EXDDD,EXD4,EXD5,0)
C
        DO 20 I=IMAX1+1,IMAXF(J)
         IDPT=IDPOT(I)
         TPT=ITPOT(I)/DTPOT
         FNRTTT=FNRTTT+G(I)*TPT*(TPT-1.D0)*(TPT-2.D0)
     *                     *TN**(TPT-3.D0)*DN**IDPT*EX
  20    CONTINUE
C
       END IF
C
       IMAX1=IMAXF(J)
  30  CONTINUE
C
      IF( IMAXF(7) .EQ. IMAXF(6) ) GOTO 50
C
      II=0
C
      DO 40 I=IMAXF(6)+1,IMAXF(7)
C
      II=II+1
C
      TPOT=ITPOT(I)/DTPOT
      BETI=BETA(I)
C
      DTN=TN-GAMMA(I)
      DDN=DN-DELTA(I)
C
      EX=DDEXP(-BETI*DTN**2-ALPHA(I)*DDN**2)
      EXT=-2.D0*BETI*DTN*EX
      EXTT=-2.D0*BETI*EX-2.D0*BETI*DTN*EXT
      EXTTT=-4.D0*BETI*EXT - 2.D0*BETI*DTN*EXTT
C
       FNRTTT=FNRTTT+G(I)* DN**IDPOT(I)
     *    *(EXTTT*TN**TPOT +3.D0*EXTT*TPOT*TN**(TPOT-1.D0)
     *    + 3.D0*EXT*TPOT*(TPOT-1.D0)*TN**(TPOT-2.D0)
     *    + EX*TPOT*(TPOT-1.D0)*(TPOT-2.D0)*TN**(TPOT-3.D0))
   40 CONTINUE
C
   50 CONTINUE
C
      IF( IMAXF(8) .EQ. IMAXF(7) ) GOTO 70
C
      TAU = (T-TNORM)/T
      DIF = (D-DNORM)/DNORM
C
      IF (DABS(TAU).LT.1.D-10) TAU=SIGN(1.D-10,TAU)
      IF (DABS(DIF).LT.1.D-10) DIF=SIGN(1.D-10,DIF)
C
      TT  = TAU*TAU
      DD  = DIF*DIF
C
      DO 60 I=IMAXF(7)+1,IMAXF(8)
C
      II = II + 1
      IJ = IMAX + (IMAXF(7)-IMAXF(6))*4 + 3 * (II-1)
C
      B  = ALPHA(I)
      BET= BETA(I)
      GA = GAMMA(I)
      GC = DELTA(I)
      GB = IDPOT(I)/10.D0
      A  = ITPOT(I)/DTPOT
C
      IF ((T.NE.TALT(I)).OR.(D.NE.DALT(I))) THEN
       DALT(I)=D
       TALT(I)=T
       DTERM(I)=(TAU+GA*DD**(.5/BET))
       DELTANA(I)=DTERM(I)**2+GB*DD**A+ETA(I)
       ETNA(I)=DDEXP(-GC*(DD+TT))
      END IF
C
      ETERM = ETNA(I)
      DELNA = DELTANA(I)
      DTE   = DTERM(I)
C
      DB    = DELNA**B
      AT3   = DB/DELNA
      DBDT  = -2.*B*DTE*AT3
      DBDTT = AT3*(2.*B + 4.*B*(B-1.)*DTE*DTE/DELNA)
      DBDTTT= B*(B-1.)*DTE*AT3/DELNA*(-12. -8.*(B-2.)*DTE*DTE/DELNA)
C
      TEX    = ETERM
      TEXDT  = 2.*GC*TAU*ETERM
      TEXDTT = (-2.*GC + 4.*GC*GC*TT)*ETERM
      TEXDTTT= (-12.*GC*GC*TAU+8.*GC*GC*GC*TT*TAU)*ETERM
C
      FNRTTTA= DN*(DB*TEXDTTT+3.*DBDT*TEXDTT+3.*DBDTT*TEXDT
     *             +DBDTTT*TEX)
C
      FNRTTT = FNRTTT + FNRTTTA * G(I)
C
   60 CONTINUE
C
   70 CONTINUE
C
      RETURN
      END
C******************************************************************************
      SUBROUTINE EFUNC(J,DN,EX,EXD,EXDD,EXDDD,EXD4,EXD5,K)
C******************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  JMINUS(5), JJMINU(5)
      COMMON/EQU/G(74),IDPOT(60),ITPOT(60),
     *           IMAXF(0:10),IMAX,DTPOT
     *          ,ALPHA(60),BETA(60),GAMMA(60),DELTA(60),ETA(60),PHI(60)
      J2 = J + J
      DO 10 I=1,5
         JMINUS(I) = J-I
         JJMINU(I) = J2-I
10    CONTINUE
        GAM1=G(IMAX+J)
        ARG1=GAM1*DN**J
        EX=DEXP(ARG1)
        IF (K .EQ. 0) GOTO 100
        EXD  = GAM1*J*DN**JMINUS(1)*EX
        IF (K .EQ. 1) GOTO 100
        EXDD = GAM1*J*(JMINUS(1)*DN**JMINUS(2)*EX
     &               +         DN**JMINUS(1)*EXD)
        IF (K .EQ. 2) GOTO 100
        EXDDD= GAM1*J*(JMINUS(1)*(JMINUS(2)*DN**JMINUS(3)*EX
     &                           +   2.D0 *DN**JMINUS(2)*EXD)
     &                +                    DN**JMINUS(1)*EXDD)
        IF (K .EQ. 3) GOTO 100
        EXd4 = gam1*j
     &         *(jminus(1)*(jminus(2)*(jminus(3)*dn**jminus(4) * ex
     &                                 + 3.d0   *dn**jminus(3) * exd)
     &                      +            3.d0   *dn**jminus(2) * exdd)
     &           +                               dn**jminus(1) * exddd)
        IF( K .EQ. 4) GOTO 100
        exd5 = gam1*j
     &         *(jminus(1)*(jminus(2)*(jminus(3)*(jminus(4)
     &                                            *dn**jminus(5)*ex
     &                                     + 4.d0 *dn**jminus(4)*exd)
     &                                 +     6.D0 *DN**JMINUS(3)*EXDD)
     &                      +                4.d0 *dn**jminus(2)*exddd)
     &           +                                 dn**jminus(1)*exd4)
100     GAM2=G(IMAX+7+J)
        IF( GAM2 .LT. -999.9D0 ) RETURN
        IF( DABS(GAM2) .LT. 1.D-10 ) THEN
         IF( DABS(ARG1) .LT. 1.D-5 ) THEN
          EX= ARG1 + .5D0*ARG1**2
         ELSE
          EX=EX -1.D0
         END IF
         RETURN
        END IF
         ARG2=GAM2*DN**J
       IF(       DABS(ARG1) .LT. 1.D-5
     *     .AND. DABS(ARG2) .LT. 1.D-5 ) THEN
         DNJ    =  DN**J
         GAMDIF = GAM1 - GAM2
         GAMSUM = GAM1 + GAM2
         GOTO (200, 201, 202, 203, 204, 205), K+1
205      EXD5  = GAMDIF * J * DN**JMINUS(5)
     &           * ( JMINUS(1) * JMINUS(2) * JMINUS(3) * JMINUS(4)
     &              + GAMSUM * JJMINU(1)*JJMINU(2)*JJMINU(3)*JJMINU(4)
     &                       * DNJ)
204      EXD4  = GAMDIF*J* DN**JMINUS(4)
     &           * (  JMINUS(1) * JMINUS(2) * JMINUS(3)
     &              + GAMSUM * JJMINU(1) * JJMINU(2) * JJMINU(3) * DNJ)
203      EXDDD = GAMDIF*J*DN**JMINUS(3)
     &           * (JMINUS(1)*JMINUS(2) +GAMSUM*JJMINU(1)*JJMINU(2)*DNJ)
202      EXDD  = GAMDIF*J*DN**JMINUS(2)*(JMINUS(1)+GAMSUM*JJMINU(1)*DNJ)
201      EXD   = GAMDIF*DN**JMINUS(1)*J* (1.D0 + GAMSUM*DNJ)
200      EX    = GAMDIF*DNJ * (1.D0 + .5D0*GAMSUM*DNJ)
       ELSE
         E2    = DEXP(ARG2)
         EX    = EX - E2
         IF (K .EQ. 0) GOTO 900
         E2D   = GAM2*J* DN**JMINUS(1) * E2
         EXD   = EXD - E2D
         IF (K .EQ. 1) GOTO 900
         E2DD  = GAM2*J* (JMINUS(1) *DN**JMINUS(2) * E2
     &                      +        DN**JMINUS(1) * E2D)
         EXDD  = EXDD - E2DD
         IF (K .EQ. 2) GOTO 900
         E2DDD = GAM2*J* (JMINUS(1) * (JMINUS(2) * DN**JMINUS(3) * E2
     &                                 +    2.D0 * DN**JMINUS(2) * E2D)
     &                    +                        DN**JMINUS(1) * E2DD)
         EXDDD = EXDDD - E2DDD
         IF (K .EQ. 3) GOTO 900
         E2D4  = GAM2*J* (JMINUS(1)* (JMINUS(2)* (JMINUS(3)
     &                                            *DN**JMINUS(4) *E2
     &                                      +3.D0* DN**JMINUS(3) *E2D)
     &                                +3.D0*       DN**JMINUS(2) *E2DD)
     &                    +                        DN**JMINUS(1) *E2DDD)
         EXD4  = EXD4 - E2D4
         IF (K .EQ. 4) GOTO 900
         E2D5  = GAM2*J* (JMINUS(1)* (JMINUS(2)* (JMINUS(3)* (JMINUS(4)
     &                   * DN**JMINUS(5) * E2
     &              +4.D0* DN**JMINUS(4) * E2D)
     &              +6.D0* DN**JMINUS(3) * E2DD)
     &              +4.D0* DN**JMINUS(2) * E2DDD)
     &              +      DN**JMINUS(1) * E2D4)
         EXD5  = EXD5 - E2D5
       END IF
900   RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION DDEXP(X)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      IF ( DABS(X) .LT. 1.D-13) THEN
       DDEXP=1.D0
       RETURN
      END IF
      IF( X .LT. -35.D0 ) THEN
       DDEXP=0.D0
       RETURN
      END IF
      DDEXP=DEXP(X)
      RETURN
      END
C*********************************************************************
           DOUBLEPRECISION FUNCTION FNI(T,D)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CFNI/B1,B2,B(18),TPID(18),TNOR,DNOR,NPOL,N,RGI
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      DN= D/DNOR
      TN= TNOR/T
      FNI= DLOG(DN) +  (B1*TN + B2) * DLOG(TN)
      DO 10 I= 1,NPOL
       FNI= FNI + B(I) * TN**TPID(I)
   10 CONTINUE
      IF(NPOL .EQ. N) GOTO 999
      DO 20 I= NPOL+1,N
       FNI= FNI + B(I) * DLOG(1.D0 - DEXP(-TPID(I)*TN))
   20 CONTINUE
  999 CONTINUE
      FNI=FNI*RGI/R
      RETURN
      END
C*********************************************************************
      DOUBLEPRECISION FUNCTION FNIT(T,D)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CFNI/B1,B2,B(18),TPID(18),TNOR,DNOR,NPOL,N,RGI
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CNORM/TNORM,DNORM
      DN= D/DNOR
      TN= TNOR/T
      FNIT = B1 + B1*DLOG(TN) + B2/TN
      DO 10 I= 1,NPOL
       FNIT= FNIT + B(I) * TPID(I) * TN**(TPID(I) -1.D0)
   10 CONTINUE
      IF(NPOL .EQ. N) GOTO 999
      DO 20 I= NPOL+1,N
       FNIT= FNIT + TPID(I)*B(I)
     *           *((1.D0-DEXP(-TPID(I)*TN))**(-1) - 1.D0)
   20 CONTINUE
  999 FNIT= FNIT * TNOR/TNORM * RGI/R
      RETURN
      END
C*********************************************************************
           DOUBLEPRECISION FUNCTION FNITT(T,D)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CFNI/B1,B2,B(18),TPID(18),TNOR,DNOR,NPOL,N,RGI
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CNORM/TNORM,DNORM
      DN= D/DNOR
      TN= TNOR/T
      FNITT = B1/TN -  B2/TN**2
      DO 10 I= 1,NPOL
       FNITT= FNITT + B(I)*TPID(I)*(TPID(I) -1.D0)*TN**(TPID(I)-2.D0)
   10 CONTINUE
      IF(NPOL .EQ. N) GOTO 999
      DO 20 I= NPOL+1,N
       XI=TPID(I)*TN
       EX=DEXP(-XI)
       FNITT= FNITT- TPID(I)**2*B(I) * (1.D0- EX)**(-2) *  EX
   20 CONTINUE
  999 FNITT= FNITT * (TNOR/TNORM)**2 * RGI/R
      RETURN
      END
C*********************************************************************
           DOUBLEPRECISION FUNCTION FNITTT(T,D)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CFNI/ B1,B2,B(18),TPID(18),TNOR,DNOR,NPOL,N,RGI
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CNORM/TNORM,DNORM
      DN= D/DNOR
      TN= TNOR/T
      FNITTT = -B1/TN**2 + 2.D0*B2/TN**3
      DO 10 I= 1,NPOL
       FNITTT= FNITTT + B(I)*TPID(I)*(TPID(I) -1.D0)
     *                    *(TPID(I)-2.D0)*TN**(TPID(I)-3.D0)
   10 CONTINUE
      IF(NPOL .EQ. N) GOTO 999
      DO 20 I= NPOL+1,N
       XI=TPID(I)*TN
       EX=DEXP(-XI)
       FNITTT = FNITTT + B(I) * TPID(I)**3 * EX* (1.D0-EX)**(-2)
     *                        * (1.D0+ 2.D0* EX / (1.D0 - EX))
   20 CONTINUE
  999 FNITTT= FNITTT * (TNOR/TNORM)**3 * RGI/R
      RETURN
      END
C*********************************************************************
      DOUBLEPRECISION FUNCTION FNID(T,D)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CFNI/ B1,B2,B(18),TPID(18),TNOR,DNOR,NPOL,N,RGI
      COMMON/CNORM/TNORM,DNORM
      TDUMMY=T
      DN=D/DNOR
      FNID=1.D0/DN
      FNID=FNID*DNORM/DNOR
      RETURN
      END
C*********************************************************************
      DOUBLEPRECISION FUNCTION FNIDD(T,D)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CNORM/TNORM,DNORM
      COMMON /CFNI/ B1,B2,B(18),TPID(18),TNOR,DNOR,NPOL,N,RGI
      DN=D/DNOR
      TDUMMY=T
      FNIDD=-1.D0/DN**2
      FNIDD=FNIDD*(DNORM/DNOR)**2
      RETURN
      END
C*********************************************************************
      DOUBLEPRECISION FUNCTION FNIDDD(T,D)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CNORM/TNORM,DNORM
      COMMON /CFNI/ B1,B2,B(18),TPID(18),TNOR,DNOR,NPOL,N,RGI
      DN=D/DNOR
      TDUMMY=T
      FNIDDD=2.D0/DN**3
      FNIDDD=FNIDDD*(DNORM/DNOR)**3
      RETURN
      END
C*******************************************************************
      SUBROUTINE WNULL3(XA,XB,F,P,T,EPS,X,IX)
C*******************************************************************
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3 ) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       GOTO 10
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      IF(F2*F1 .LE.F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
      X=(X1+X3)/2.D0
      F2=F(X,T,P)
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
  10  IX=1
 999  RETURN
      END
C*********************************************************************
      DOUBLEPRECISION FUNCTION RSMXW(T,DL,DV)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      RSMXW=(FNR(T,DL)-FNR(T,DV)+DLOG(DL/DV))*R*T
      RETURN
      END
C*******************************************************************
      SUBROUTINE WNULL2(XA,XB,F,P,T,EPS,X,IX)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F
      X1=XA
      F1=F(X1,T,P)
      X3=XB
      F3=F(X3,T,P)
      IX=0
      IF((X1-X3)*(F1-F3) .LT. 0.D0) THEN
        IX=2
        RETURN
      END IF
C++ SCHLEIFE ZUR ERMITTLUNG DER NULLSTELLE +++++++++++++++++++++++++++
      DO 100 I=1,40
      IF( F1 .NE. F3) THEN
       X=X1+(X3-X1)*F1/(F1-F3)
      ELSE
       GOTO 10
      END IF
      IF(X .LT. 0.D0) X=(X1+X3)/2.D0
      IF(DABS(X) .LT. 1.D-8) THEN
        IF(DABS(X-X1).LT.EPS) RETURN
      ELSE
        IF(DABS((X-X1)/X).LT.EPS) RETURN
      END IF
      F2=F(X,T,P)
      X2=X1-(X1-X3)/2.D0
      IF(F2*F1 .LE.F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
      IF((X2-X3)*(X2-X1) .GE. 0.D0) GOTO 100
      X=(X1+X3)/2.D0
      F2=F(X,T,P)
      IF(F2*F1 .LE. F2*F3) THEN
        X3=X1
        F3=F1
      END IF
      X1=X
      F1=F2
 100  CONTINUE
C++ ENDE DER SCHLEIFE ++++++++++++++++++++++++++++++++++++++++++++++++
  10  IX=1
 999  RETURN
      END
C*******************************************************************
      SUBROUTINE WNULL(XA,XB,F,P,T,EPS,ITX,Z,N,I,ITEST)
C*******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER N,I
      DIMENSION Z(I)
      EXTERNAL F
      IDUMMY = ITEST
      H=(XB-XA)/ITX
      N=0
      X1=XA
      F1=F(X1,T,P)
      F3=F1
  5   X3=X1+H
      IF (DABS(X3-XB).LT.DABS(H/2.D0)) X3=XB-1.D-8*H
      XCON=(XB-XA)*(XB-X3)
      IF(XCON.LE.0.D0) GOTO 999
      EX=X3
      F1=F3
      F3=F(X3,T,P)
  3   IF(F1*F3.GT.0.D0) GOTO 4
      X=X3
 100  CONTINUE
      X2=(X1+X3)/2.D0
      F2=F(X2,T,P)
      IF(F1*F2.GT.0.D0) GOTO 11
      X3=X2
      F3=F2
      GOTO 12
  11  X1=X2
      F1=F2
  12  X2=X
      X=X1+(X3-X1)*F1/(F1-F3)
      F2=F(X,T,P)
  13  IF(DABS(F2) .LT. EPS) GOTO 16
      IF(DABS(X) .GT. 1.D-8) GOTO 17
      IF(DABS(X-X2).LT. EPS) GOTO 16
      GOTO 18
  17  IF(DABS((X-X2)/X).LT. EPS) GOTO 16
  18  X2=X
      IF(F1*F2.GT.0.D0) GOTO 14
      X3=X2
      F3=F2
      GOTO 15
  14  X1=X2
      F1=F2
  15  GOTO 100
  16  N=N+1
      Z(N)=X
      IF(N.GE.I) GOTO 999
      IF((EX-X).LT.EPS) X3=EX+10.D0*EPS
      IF((EX-X).GE.EPS) X3=EX
      F3=F(X3,T,P)
  4   X1=X3
      GOTO 5
  999 RETURN
      END
C*********************************************************************
      DOUBLEPRECISION FUNCTION NULLF(D,T,P)
C*********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      NULLF=PB(T,D)-P
  999 RETURN
      END
C*********************************************************************
      SUBROUTINE SATCAP(T,DV,DL,P,EPS)
C*********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 NUFU2,NUFU3
      COMMON/CSUB3/ TTR,PTR,DLTR,DVTR,TB,PBX,DLB,DVB
      COMMON/COUT/NIN,NOUT
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CCPEQ/ TCEQ,PCEQ,DCEQ
      EXTERNAL NUFU2
      EXTERNAL NUFU3
      IH=0
      T=0.D0
      DL=0.D0
      DV=0.D0
      IF (P.GT.PCEQ) THEN
       WRITE(NOUT,30)
   30  FORMAT('<SATP> P >PC')
       GOTO 999
      ENDIF
      IF ((DABS(PCEQ-P)).LT.1.D-10) THEN
       T= TCEQ
       DL=DCEQ
       DV=DCEQ
       GOTO 999
      ENDIF
      IF (EPS.LT.1.D-7) EPS=1.D-7
      IF( P .GE. .998D0*PC ) THEN
       PD= PC - (PCEQ - P)
       CALL SATPEQ(TEST,DVEST,DLEST,PD,EPS)
       TEST= TCEQ - (TC - TEST)
      ELSE
       CALL SATPEQ(TEST,DVEST,DLEST,P,EPS)
      END IF
      IF( TEST .LT. 1.D-10) GOTO 999
      T1=TEST - 1.D1*EPS
      T2=TEST
      IF(T2 .GT. TCEQ) THEN
       T2=TCEQ-.05D0
       T1=T2-10.D0*EPS
      END IF
  50  CALL WNULL3(T2,T1,NUFU2,P,EPS,EPS,TEST,IY)
      IF (IY.EQ.0) THEN
       T=TEST
       CALL SATCAL (T,DV,DL,PSE,EPS)
      ELSE  IF(DABS(P-PCEQ).LT.0.01D0) THEN
       T2= TCEQ-.05D0
       T1= T2- 1.D1*EPS
       CALL WNULL3(T2,T1,NUFU3,P,DCEQ,EPS,Z,IZ)
       IF(IZ .EQ. 0 .AND. IH .NE. 1) THEN
        IH=1
        T2=Z
        T1=T2-.05D0
        GOTO 50
       END IF
       IF(IZ .NE. 0) GOTO 999
       T=Z
      END IF
  999 CONTINUE
      RETURN
      END
C*********************************************************************
      DOUBLEPRECISION FUNCTION NUFU2(TSZ,EPS,P)
C*********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CCPEQ/TCEQ,PCEQ,DCEQ
      IF( TSZ .GT. TCEQ) THEN
       PSZ=PCEQ
      ELSE
       CALL SATCAL (TSZ,DV,DL,PSZ,EPS)
      END IF
      NUFU2= PSZ-P
      RETURN
      END
C*********************************************************************
      DOUBLEPRECISION FUNCTION NUFU3(TSZ,RHO,P)
C*********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CCPEQ/TCEQ,PCEQ,DCEQ
      NUFU3= PB(TSZ,RHO)-P
      RETURN
      END
C*********************************************************************
      SUBROUTINE SATCAL(T,DV,DL,P,EPS)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CSUB2/RG,XMOL,TC,PC,DC
      COMMON/CCPEQ/TCEQ,PCEQ,DCEQ
      COMMON/COUT/NIN,NOUT
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PBOIL,DLB,DVB
      SAVE TOLD,DLOLD,DVOLD,PSOLD
      DATA TOLD /1.D9/
C
      IF (DABS(T-TOLD).LT.1.D-6) THEN
       DV = DVOLD
       DL = DLOLD
       P  = PSOLD
       RETURN
      END IF
C
      IF( TCEQ .GT. 1.D0 ) THEN
        TCE=TCEQ
        DCE=DCEQ
        PCE=PCEQ
       ELSE
        TCE=TC
        DCE=DC
        PCE=PC
       END IF
       IF( T .GT. TCE ) GOTO 505
      ICMAX=0
      ICMIN=0
      PMAX=1000.D0
      PMIN=0.D0
      DMIN=DCE
      DMAX=0.D0
      DMAX1=0.D0
      ITX=100
      PX1=2.D0
      PX2=1.D0
      EPSD=EPS/10.D0
      TD=T
      IF(TC-T .LT. 0.5D0) TD=TC-0.5D0
      CALL DVEQN(TD,DVEST)
      CALL DLEQN(TD,DLEST)
      TX=(T-TTR)/(TC-TTR)
      IF(TX .GT. 0.995D0) GOTO 11
      CALL VPEQN(T,PEST)
      GOTO 30
  11  CONTINUE
       PEST=PB(T,DCE)
      GOTO 30
  12  CONTINUE
      IF(ICMAX .EQ. 0) GOTO 13
      IF( DMAX1 .GT. 1.D-2 ) THEN
       DMAX=DMAX1
       DMAX1=0.D0
       GOTO 130
      ELSE
      WRITE(NOUT,1003)
 1003 FORMAT('<SATT> DID NOT CONVERGE!')
      GOTO 500
      END IF
  13  DV1=DVEST*1.D0
      DL2=DLEST*1.D0
      DELD=(DL2-DV1)/ITX
      D=DV1-DELD
      DO 100 I=1,ITX
      D=D+DELD
      P=PB(T,D)
      DPX1=P-PX1
      DPX2=PX1-PX2
      DPDP=DPX1*DPX2
      IF(DPDP.GT.0.D0) GOTO 110
      IF(DPX1.GT.0.D0) GOTO 120
       ICMAX=1
       IF(PX1.LT.PMAX) PMAX=PX1
       IF( D .LT. DMIN) DMIN=D-DELD
      GOTO 110
 120  CONTINUE
       ICMIN=1
       IF(PX1.GT.PMIN) PMIN=PX1
       IF(D .GT. DMAX) THEN
        DMAX1=DMAX
        DMAX=D
       END IF
 110  CONTINUE
      PX2=PX1
      PX1=P
 100  CONTINUE
 130  IF(ICMIN.NE.1.OR.ICMAX.NE.1) GOTO 505
      PEST=(PMAX+PMIN)/2.D0
      DVEST=DMIN - 0.2D0*(DMAX-DMIN)
      DLEST=DMAX + 0.2D0*(DMAX-DMIN)
      GOTO 30
  33  PEST=P
  30  CONTINUE
      P1=PEST
      DL=DENS(P1,T,DLEST,EPSD)
      IF(DL.LT.1.D-20) GOTO 12
      DV=DENS(P1,T,DVEST,EPSD)
      IF(DV.LT.1.D-20) GOTO 12
      IF( DABS((DL-DV)/DV).LE.10.D0*EPSD) GOTO 12
      DELV=1.D0/DV-1.D0/DL
      DFR=RSMXW(T,DL,DV)
      VAR1=P1*DELV-DFR
      P=DFR/DELV
      IF(P.GT.0.D0) GOTO 301
      P1=P1/5.D0
      PEST=P1
      GOTO 30
  301 CONTINUE
      DL0=DL
      DV0=DV
      IF(DL.LE.DCE) DL=DLEST
      IF(DV.GE.DCE) DV=DVEST
      DL=DENS(P,T,DL,EPSD)
      IF(DL.LT.1.D-20) GOTO 12
      DV=DENS(P,T,DV,EPSD)
      IF(DV.LT.1.D-20) GOTO 12
      IF( DABS((DL-DV)/DV).LE.10.D0*EPSD) GOTO 12
      DDL=DABS((DL-DL0)/DL)
      DPS=DABS((P-P1)/P)
      DD=DDL**2+DPS**2
      IF(DD.GT. EPS**2) GOTO 31
      DDV=DABS((DV-DV0)/DV)
      IF(DDV.LT.EPS) GOTO 999
   31 CONTINUE
      DELV=1.D0/DV-1.D0/DL
      DFR=RSMXW(T,DL,DV)
      VAR3=P*DELV-DFR
      P3=P
      DO 300 I=1,30
      P=P1+(P3-P1)*VAR1/(VAR1-VAR3)
      DL0=DL
      DV0=DV
      IF(DL.LE.DCE) DL=DLEST
      IF(DV.GE.DCE) DV=DVEST
      DL=DENS(P,T,DL,EPSD)
      IF(DL.LT.1.D-20) GOTO 12
      DV=DENS(P,T,DV,EPSD)
      IF(DV.LT.1.D-20) GOTO 12
      IF( DABS((DL-DV)/DV).LE.10.D0*EPSD) GOTO 12
      DDL=DABS((DL-DL0)/DL)
      DPS=DABS((P-P1)/P)
      DD=DDL**2+DPS**2
      IF(DD.GT. EPS**2) GOTO 310
      DDV=DABS((DV-DV0)/DV)
      IF(DDV.LT.EPS) GOTO 999
  310 CONTINUE
      DELV=1.D0/DV-1.D0/DL
      DFR=RSMXW(T,DL,DV)
      VAR2=P*DELV-DFR
      IF(VAR1*VAR2.GT.VAR2*VAR3) GOTO 320
      P3=P1
      VAR3=VAR1
 320  P1=P
      VAR1=VAR2
 300  CONTINUE
      XX1=1.D0
      WRITE(NOUT,1004)
 1004  FORMAT('<SATCAL> DID NOT CONVERGE!')
 500  CONTINUE
      DL=0.D0
      DV=0.D0
      P=0.D0
      GOTO 999
 505  CONTINUE
      P=0.D0
      DL=0.D0
      DV=0.D0
      WRITE(NOUT,1005)
 1005 FORMAT('<SATCAL> T>TC')
  999 CONTINUE
      TOLD = T
      DLOLD = DL
      DVOLD = DV
      PSOLD = P
      RETURN
      END
C*********************************************************************
      DOUBLEPRECISION FUNCTION DB(T,P)
C*********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PB,DLB,DVB
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/COUT/NIN,NOUT
      IF( T .LE. TTR.OR. P .LE. 0.D0 ) THEN
       WRITE(NOUT,1001) T,P
 1001  FORMAT('<DB> T=',F8.3,'or P=',D10.4,'out of range of validity')
       RETURN
      END IF
        DEST= BDENS(T,P,0)
        DB= DENS(P,T,DEST,1.D-6)
      RETURN
      END
C*********************************************************************
      DOUBLEPRECISION FUNCTION PB(T,D)
C*********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CNORM/TNORM,DNORM
      PB=0.D0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      DN=D/DNORM
      PB=D*R*T*(1.D0+DN*FNRD(T,D))
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION DPDTB(T,D)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CNORM/TNORM,DNORM
      DPDTB=0.D0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      TN=TNORM/T
      DN=D/DNORM
      DPDTB=D*R*(1.D0+DN*FNRD(T,D)-DN*TN*FNRDT(T,D))
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION DPDDB(T,D)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CNORM/TNORM,DNORM
      DPDDB=0.D0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      DN=D/DNORM
      DPDDB=R*T*(1.D0+2.D0*DN*FNRD(T,D)+DN**2*FNRDD(T,D))
      RETURN
      END
C*********************************************************************
      DOUBLEPRECISION FUNCTION FUGAB(T,D)
C*********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CNORM/TNORM,DNORM
      FUGAB=0.D0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      DN=D/DNORM
      ZB=(1.D0+DN*FNRD(T,D))
      PB=ZB*D*R*T
C
      IF (ZB.LT.0.D0) RETURN
C
      F=ZB-1.D0+FNR(T,D)-DLOG(ZB)
      FUGAB=PB*DEXP(F)
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION CVB(T,D)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CSUB4/FW,FC
      COMMON/CNORM/TNORM,DNORM
      CVB=0.D0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      TN=TNORM/T
      CVB=-FC*R*TN**2*(FNITT(T,D)+FNRTT(T,D))
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION CPB(T,D)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CSUB4/FW,FC
      COMMON/CNORM/TNORM,DNORM
      CPB=0.D0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      TN=TNORM/T
      DN=D/DNORM
      CPB=CVB(T,D)+FC*R*(1.D0+DN*FNRD(T,D)-DN*TN*FNRDT(T,D))**2/
     1    (1.D0+2.D0*DN*FNRD(T,D)+DN**2*FNRDD(T,D))
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION WB(T,D)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CSUB4/FW,FC
      COMMON/CNORM/TNORM,DNORM
      WB=0.D0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      TN=TNORM/T
      DN=D/DNORM
      WB2=CPB(T,D)/CVB(T,D)*(1.D0+2.D0*DN*FNRD(T,D)+DN**2*FNRDD(T,D))
      IF(WB2 .GT. 0.D0) THEN
       WB=DSQRT(WB2*R*FW*T)
      ELSE
       WB=0.D0
      END IF
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION FB(T,D)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CNORM/TNORM,DNORM
      COMMON/CSUB4/FW,FC
      COMMON /CNULLP/T0,D0,F0,H0,U0,S0,G0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      FB=R*FC*T * (FNR(T,D) + FNI(T,D))
      FB=FB+F0
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION HB(T,D)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CNORM/TNORM,DNORM
      COMMON/CSUB4/FW,FC
      COMMON /CNULLP/T0,D0,F0,H0,U0,S0,G0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      TN=TNORM/T
      DN=D/DNORM
      HB=R*FC*T* (1.D0 + DN*FNRD(T,D) + TN*( FNIT(T,D) + FNRT(T,D)))
      HB=HB+H0
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION SB(T,D)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CNORM/TNORM,DNORM
      COMMON/CSUB4/FW,FC
      COMMON /CNULLP/T0,D0,F0,H0,U0,S0,G0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      TN=TNORM/T
      SB=R*FC * ( TN*(FNIT(T,D)+FNRT(T,D)) - (FNI(T,D)+FNR(T,D)) )
      SB=SB+S0
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION UB(T,D)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CNORM/TNORM,DNORM
      COMMON/CSUB4/FW,FC
      COMMON /CNULLP/T0,D0,F0,H0,U0,S0,G0
      UB=0.D0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      TN=TNORM/T
      UB=R*FC*T * TN*(FNIT(T,D)+FNRT(T,D))
      UB=UB+U0
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION GB(T,D)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CNORM/TNORM,DNORM
      COMMON/CSUB4/FW,FC
      COMMON /CNULLP/T0,D0,F0,H0,U0,S0,G0
      GB=0.D0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      TN=TNORM/T
      DN=D/DNORM
      GB=R*FC*T * (1.D0+DN*FNRD(T,D)+FNI(T,D)+FNR(T,D))
      GB=GB+G0
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION BB(T)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CNORM/TNORM,DNORM
      DATA D /1.D-200/
      BB=FNRD(T,D)/DNORM
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION CB(T)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CNORM/TNORM,DNORM
      DATA D /1.D-200/
      CB=FNRDD(T,D)/DNORM**2
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION RJTB(T,D)
C***********************************************************************
      IMPLICIT REAL * 8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CSUB4/FW,FC
      COMMON/CNORM/TNORM,DNORM
      RJTB=0.D0
      IF( T .LT. 1.D0 ) RETURN
      IF( D .LT. 1.D-10) RETURN
      TN=TNORM/T
      DN=D/DNORM
      CVBN=-TN**2*(FNITT(T,D)+FNRTT(T,D))
      DPDDN=1.D0+2.D0*DN*FNRD(T,D)+DN**2*FNRDD(T,D)
      RJTB=-((D*R)**(-1))*(DN*FNRD(T,D)+DN**2*FNRDD(T,D)
     *        +DN*TN*FNRDT(T,D))/
     *     ((1.D0+DN*FNRD(T,D)-DN*TN*FNRDT(T,D))**2+CVBN*DPDDN)
      RETURN
      END
C************************************************************************
      DOUBLEPRECISION FUNCTION DENS(P,T,DEST,EPS)
C************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 NULLF
      DIMENSION Z(3)
      COMMON/CCPEQ/TCEQ,PCEQ,DCEQ
      COMMON/CSUB2/RG,XMOL,TC,PC,DC
      EXTERNAL NULLF
      IF( DCEQ .GT. 1.D0) THEN
       DC1=DCEQ
      ELSE
       DC1=DC
      END IF
      IF(DEST.GT.2.D0*DC1) THEN
       DD=.01D0*DEST
      ELSE
       DD=0.05D0*DEST
      END IF
      D1=DEST
      IF(DEST .GT. DC1) THEN
       D2=D1+DD
      ELSE
       D2=D1-DD
      END IF
      IF(EPS. LT. 1.D-10) EPS=1.D-10
      CALL WNULL2(D1,D2,NULLF,P,T,EPS,X,IX)
      IF( IX .LE. 0 ) THEN
        DENS=X
      ELSE
        D1=DEST-5*DD
        D2=DEST+5*DD
        CALL WNULL(D1,D2,NULLF,P,T,EPS,20,Z,N,3,IX)
        IF( N .GT. 0 ) THEN
          IF( DEST .GT. DC1 ) THEN
            DENS=Z(N)
          ELSE
            DENS=Z(1)
          END IF
        ELSE
          DENS=0.D0
       END IF
      END IF
 999  RETURN
      END
C***************************************************************
      DOUBLEPRECISION FUNCTION BDENS(T,P,IREG)
C***************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CSUB2/RG,XMOL,TC,PC,DC
      COMMON/CCPEQ/TCEQ,PCEQ,DCEQ
      IF( TCEQ .GT. 1.D0 ) THEN
       DC1=DCEQ
       PC1=PCEQ
       TC1=TCEQ
      ELSE
       DC1=DC
       PC1=PC
       TC1=TC
      END IF
      DXC=DSOAVE(TC1,PC1,1)
      DDC=DC1-DXC
      D=DSOAVE(T,P,IREG)
      IF(D+DDC .GT. DC1) D=D+DDC
      BDENS=D
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION DSOAVE(T,P,IREG)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(3)
      COMMON/CSUB2/ RG,XMOL,TC,PC,DC
      COMMON/COUT/NIN,NOUT
      DSOAVE=0.D0
      IF(P.LT.0.D0.OR.T.LT.0.D0) THEN
      WRITE(NOUT,1000)
 1000 FORMAT('CALCULATION OF NOT POSSIBLE, CHECK P AND T')
      GOTO 999
      ENDIF
      IF(P. LT. 1.D-4) P=1.D-4
      PI= DACOS(0.D0)*2.D0
      ATC= .42748025D0*(RG*TC)**2/PC
      B= .08664035D0*RG*TC/PC
      TSOAVE= .7D0*TC
      CALL VPEQN(TSOAVE,PSOAVE)
      PRSOAVE= PSOAVE/PC
      OMEGA= -DLOG10(PRSOAVE) -1.D0
      EM= .47979D0 +1.576D0*OMEGA -.1925D0*OMEGA**2 +.025D0*OMEGA**3
      ALPHA=(1.D0+ EM*(1.D0-(T/TC)**.5D0))**2
      A= ALPHA*ATC
      BETA= -(RG*T)/P
      GAMMA= -(B**2)-(B*RG*T)/P+A/P
      DELTA= -(A*B)/P
      Q= (BETA**3)/27.D0-(BETA*GAMMA)/6.D0+DELTA/2.D0
      PE= (GAMMA-(BETA**2)/3.D0)/3.D0
      DISK= (Q**2)+(PE**3)
      IF(DABS(PE).LT.1.D-20) THEN
       Y(1)= ((DABS(2.D0*Q))**(1.D0/3.D0))
       IF(Q.GT.0.D0) Y(1)=-Y(1)
       Y(2)= 0.D0
       Y(3)= 0.D0
       KK=1
       GOTO 500
      ENDIF
      IF(DISK.LT.0.D0) THEN
       R= ((DABS(PE))**.5D0)
      IF(Q.LT.0.D0) R=-R
       PHI= DACOS(Q/(R**3))
       Y(1)= -(2.D0*R)*DCOS(PHI/3.D0)
       Y(2)=  (2.D0*R)*DCOS((PI/3.D0)-(PHI/3.D0))
       Y(3)=  (2.D0*R)*DCOS((PI/3.D0)+(PHI/3.D0))
       KK= 3
      ELSE
       DISKN=(DISK**.5D0)
       SUM1= -Q+DISKN
        IF(SUM1.LT.0.D0) THEN
         U=-((DABS(SUM1))**(1.D0/3.D0))
        ELSE
         U= (SUM1**(1.D0/3.D0))
        ENDIF
       SUM2= -Q-DISKN
        IF(SUM2.LT.0.D0)THEN
         V=-((DABS(SUM2))**(1.D0/3.D0))
        ELSE
         V= (SUM2**(1.D0/3.D0))
        ENDIF
       Y(1)= U+V
       Y(2)= 0.D0
       Y(3)= 0.D0
       KK=1
      ENDIF
  500 CONTINUE
      DMAX=0.D0
      DMIN=1.D100
      DO 50 I=1,KK
       D=1.D0/(Y(I)-BETA/3.D0)
       IF(D.GT.DMAX) DMAX=D
       IF(D.LT.DMIN) DMIN=D
   50 CONTINUE
      DSOAVE=DMIN
      IF(DMIN.LT.0.D0) DSOAVE=DMAX
      IF(T .LT. TC) THEN
       IREGH=IREG
       IF(IREG .EQ. 0) THEN
        CALL VPEQN(T,PS)
        IREGH=2
        IF(P .GT. PS) IREGH=1
       END IF
       IF(IREGH .EQ. 1) THEN
        DSOAVE=DMAX
        CALL DLEQN(T,DL)
       ELSE
        DSOAVE=DMIN
        CALL DVEQN(T,DV)
       END IF
      END IF
  999 CONTINUE
      RETURN
      END
C*********************************************************************
      SUBROUTINE SATPEQ(TS,DV,DL,P,EPS)
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8   NUFU
      COMMON/COUT/NIN,NOUT
      COMMON/CSUB2/RG,XMOL,TC,PC,DC
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TBB,PBB,DLBB,DVBB
      EXTERNAL NUFU
      TS=0.D0
      DV=0.D0
      DL=0.D0
      IF( DABS(P-PC) .LT. 1.D-6 ) THEN
       DV=DC
       DL=DC
       TS=TC
       P=PC
       GOTO 999
      END IF
      IF(P.GT.PC) THEN
       GOTO 999
      ENDIF
      IF(P.LT.PTR) THEN
       GOTO 999
      ENDIF
      XKONST= DLOG(PC/PTR)/(1.D0/TTR-1.D0/TC)
      T1=(1.D0/TC-1.D0/XKONST*DLOG(P/PC))**(-1)
      T2= T1-2.D0
      CALL WNULL2(T1,T2,NUFU,P,EMPTY,EPS,TS,IX)
      IF(IX.NE.0) THEN
       WRITE(NOUT,50)
   50  FORMAT('<SATPEQ> CANNOT FIND TS(P)')
       TS=0.D0
       GOTO 999
      ENDIF
      CALL DVEQN(TS,DV)
      CALL DLEQN(TS,DL)
  999 CONTINUE
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUFU(T,EMPTY,P)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DUMMY=EMPTY
      CALL VPEQNL(T,PS)
      NUFU= PS-P
      RETURN
      END
C***********************************************************************
      DOUBLEPRECISION FUNCTION NUFUD(T,EMPTY,D)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      DUMMY = EMPTY
      IF(D .GT. DC) THEN
       CALL DLEQN(T,DB)
       NUFUD=D-DB
      ELSE
       CALL DVEQN(T,DB)
       NUFUD=DB-D
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE VPEQNL(T,P)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PBX,DLB,DVB
      IF(T.GT.TC) THEN
       T1=TC-.1D0
       CALL VPEQN(T1,P1)
       P=PC+(PC-P1)/(TC-T1)*(T-TC)
      ELSE
       CALL VPEQN(T,P)
      END IF
      RETURN
      END
C*********************************************************************
      BLOCKDATA FNIDAT
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /CFNI/B1,B2,B(18),TPID(18),TNOR,DNOR,NPOL,N,RGI
      COMMON /CNULLP/T0,D0,F0,H0,U0,S0,G0
      DATA T0/0.298150D+03/D0/0.179885D+01/F0/0.000000D+00
     *    /H0/0.000000D+00/U0/0.000000D+00/S0/0.000000D+00
     *    /G0/0.000000D+00/
      DATA TNOR/0.304128200D+03/DNOR/0.467600000D+03
     *     /RGI/0.188924100D-03/NPOL/ 2/N/ 7/
      DATA B1/0.000000000D+00/B2/0.250000000D+01/
     *     B /0.837304456D+01,-.370454304D+01,0.832767753D-01,
     *        0.199427042D+01,0.104028922D+01,0.621052475D+00,
     *        0.411952928D+00,0.000000000D+00,0.000000000D+00,
     *        0.000000000D+00,0.000000000D+00,0.000000000D+00,
     *        0.000000000D+00,0.000000000D+00,0.000000000D+00,
     *        0.000000000D+00,0.000000000D+00,0.000000000D+00/
      DATA TPID/0.000000000D+00,0.100000000D+01,0.270879200D+02,
     *          0.315163000D+01,0.113238400D+02,0.611190000D+01,
     *          0.677708000D+01,0.000000000D+00,0.000000000D+00,
     *          0.000000000D+00,0.000000000D+00,0.000000000D+00,
     *          0.000000000D+00,0.000000000D+00,0.000000000D+00,
     *          0.000000000D+00,0.000000000D+00,0.000000000D+00/
      END
C*********************************************************************
      BLOCKDATA  SUBST
C*********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CSUB2/RG,XMOL,TC,PC,DC
      COMMON/CSUB21/R1,XMOL1,TC1,PC1,DC1
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PB,DLB,DVB
      COMMON/CSUB4/FW,FC
      COMMON/EQNN/GN(23,3),TPOTN(20,3),IMAXN(3),EQKNN(5,3)
      COMMON/CNORMN/TNORMN(3),YNORMN(3),ISNN(3)
      COMMON/CCPEQN/TCEQN(3),PCEQN(3),DCEQN(3)
      COMMON/CISN/ISLN(3),ISRN(3)
      COMMON/H/HCPN(3)
      DATA RG   /0.188924100D-03/R1   /0.188924100D-03/
     *     XMOL /0.440098000D-01/XMOL1/0.440098000D-01/
     *     TC   /0.304128200D+03/TC1  /0.304128200D+03/
     *     PC   /0.737730000D+01/PC1  /0.737730000D+01/
     *     DC   /0.467600000D+03/DC1  /0.467600000D+03/
      DATA TTR/0.216592000D+03/PTR/0.517950000D+00/DLTR/0.117853000D+04
     *   /DVTR/0.137614000D+02/TB/0.000000000D+00/PB/0.000000000D+00/
      DATA FW/0.100000000D+07/FC/0.100000000D+04/
      DATA ISLN(1)/ 3/ISRN(1)/ 1/ISNN(1)/ 1/IMAXN(1)/ 4/
      DATA TPOTN( 1,1)/0.100000000D+01/GN( 1,1)/-.706020874D+01/
      DATA TPOTN( 2,1)/0.150000000D+01/GN( 2,1)/0.193912184D+01/
      DATA TPOTN( 3,1)/0.200000000D+01/GN( 3,1)/-.164635966D+01/
      DATA TPOTN( 4,1)/0.400000000D+01/GN( 4,1)/-.329956340D+01/
      DATA TCEQN(1)/0.304130000D+03/HCPN(1)/0.000000000D+00
     *    /TNORMN(1)/0.304128200D+03/YNORMN(1)/0.737730000D+01/
      DATA ISLN(2)/ 1/ISRN(2)/ 1/ISNN(2)/ 1/IMAXN(2)/ 4/
      DATA TPOTN( 1,2)/0.340000000D+00/GN( 1,2)/0.192451083D+01/
      DATA TPOTN( 2,2)/0.500000000D+00/GN( 2,2)/-.623855554D+00/
      DATA TPOTN( 3,2)/0.166666667D+01/GN( 3,2)/-.327311265D+00/
      DATA TPOTN( 4,2)/0.183333333D+01/GN( 4,2)/0.392451423D+00/
      DATA TCEQN(2)/0.304130000D+03/HCPN(2)/0.000000000D+00
     *    /TNORMN(2)/0.304128200D+03/YNORMN(2)/0.467600000D+03/
      DATA ISLN(3)/ 1/ISRN(3)/ 1/ISNN(3)/ 1/IMAXN(3)/ 5/
      DATA TPOTN( 1,3)/0.340000000D+00/GN( 1,3)/-.170748787D+01/
      DATA TPOTN( 2,3)/0.500000000D+00/GN( 2,3)/-.822746702D+00/
      DATA TPOTN( 3,3)/0.100000000D+01/GN( 3,3)/-.460085488D+01/
      DATA TPOTN( 4,3)/0.233333333D+01/GN( 4,3)/-.101111780D+02/
      DATA TPOTN( 5,3)/0.466666666D+01/GN( 5,3)/-.297422520D+02/
      DATA TCEQN(3)/0.304130000D+03/HCPN(3)/0.000000000D+00
     *    /TNORMN(3)/0.304128200D+03/YNORMN(3)/0.467600000D+03/
      END
C******************************************************************
      BLOCKDATA  EQS
C******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/EQU/G(74),IDPOT(60),ITPOT(60),
     *           IMAXF(0:10),IMAX,DTPOT
     *          ,ALPHA(60),BETA(60),GAMMA(60),DELTA(60),ETA(60),PHI(60)
      COMMON/CNORM/TNORM,DNORM
      COMMON/CNORM1/TNORM1,DNORM1
      COMMON/CCPEQ/TCEQ,PCEQ,DCEQ
      COMMON/CCPEQ1/TCEQ1,PCEQ1,DCEQ1
      COMMON/CLIM/TLIM,PLIM
      DATA DTPOT/0.80D+01/IMAX/42/
      DATA IMAXF /  7, 16, 23, 26, 32, 33, 34, 39, 42,  0,  0/
      DATA G( 1)/0.388568232032D+00/IDPOT( 1)/  1/ITPOT( 1)/  0/
      DATA G( 2)/0.293854759427D+01/IDPOT( 2)/  1/ITPOT( 2)/  6/
      DATA G( 3)/-.558671885349D+01/IDPOT( 3)/  1/ITPOT( 3)/  8/
      DATA G( 4)/-.767531995925D+00/IDPOT( 4)/  1/ITPOT( 4)/ 16/
      DATA G( 5)/0.317290055804D+00/IDPOT( 5)/  2/ITPOT( 5)/  6/
      DATA G( 6)/0.548033158978D+00/IDPOT( 6)/  2/ITPOT( 6)/ 16/
      DATA G( 7)/0.122794112203D+00/IDPOT( 7)/  3/ITPOT( 7)/  6/
      DATA G( 8)/0.216589615432D+01/IDPOT( 8)/  1/ITPOT( 8)/ 12/
      DATA G( 9)/0.158417351097D+01/IDPOT( 9)/  2/ITPOT( 9)/ 12/
      DATA G(10)/-.231327054055D+00/IDPOT(10)/  4/ITPOT(10)/ 20/
      DATA G(11)/0.581169164314D-01/IDPOT(11)/  5/ITPOT(11)/  0/
      DATA G(12)/-.553691372054D+00/IDPOT(12)/  5/ITPOT(12)/ 12/
      DATA G(13)/0.489466159094D+00/IDPOT(13)/  5/ITPOT(13)/ 16/
      DATA G(14)/-.242757398435D-01/IDPOT(14)/  6/ITPOT(14)/  0/
      DATA G(15)/0.624947905017D-01/IDPOT(15)/  6/ITPOT(15)/  8/
      DATA G(16)/-.121758602252D+00/IDPOT(16)/  6/ITPOT(16)/ 16/
      DATA G(17)/-.370556852701D+00/IDPOT(17)/  1/ITPOT(17)/ 24/
      DATA G(18)/-.167758797004D-01/IDPOT(18)/  1/ITPOT(18)/ 48/
      DATA G(19)/-.119607366380D+00/IDPOT(19)/  4/ITPOT(19)/ 24/
      DATA G(20)/-.456193625088D-01/IDPOT(20)/  4/ITPOT(20)/ 48/
      DATA G(21)/0.356127892703D-01/IDPOT(21)/  4/ITPOT(21)/ 64/
      DATA G(22)/-.744277271321D-02/IDPOT(22)/  7/ITPOT(22)/ 48/
      DATA G(23)/-.173957049024D-02/IDPOT(23)/  8/ITPOT(23)/  0/
      DATA G(24)/-.218101212895D-01/IDPOT(24)/  2/ITPOT(24)/ 56/
      DATA G(25)/0.243321665592D-01/IDPOT(25)/  3/ITPOT(25)/ 96/
      DATA G(26)/-.374401334235D-01/IDPOT(26)/  3/ITPOT(26)/128/
      DATA G(27)/0.143387157569D+00/IDPOT(27)/  5/ITPOT(27)/176/
      DATA G(28)/-.134919690833D+00/IDPOT(28)/  5/ITPOT(28)/192/
      DATA G(29)/-.231512250535D-01/IDPOT(29)/  6/ITPOT(29)/128/
      DATA G(30)/0.123631254929D-01/IDPOT(30)/  7/ITPOT(30)/192/
      DATA G(31)/0.210583219729D-02/IDPOT(31)/  8/ITPOT(31)/ 64/
      DATA G(32)/-.339585190264D-03/IDPOT(32)/ 10/ITPOT(32)/ 16/
      DATA G(33)/0.559936517716D-02/IDPOT(33)/  4/ITPOT(33)/224/
      DATA G(34)/-.303351180556D-03/IDPOT(34)/  8/ITPOT(34)/112/
      DATA G(35)/-.213654886883D+03/IDPOT(35)/  2/ITPOT(35)/  8/
      DATA G(36)/0.266415691493D+05/IDPOT(36)/  2/ITPOT(36)/  0/
      DATA G(37)/-.240272122046D+05/IDPOT(37)/  2/ITPOT(37)/  8/
      DATA G(38)/-.283416034240D+03/IDPOT(38)/  3/ITPOT(38)/ 24/
      DATA G(39)/0.212472844002D+03/IDPOT(39)/  3/ITPOT(39)/ 24/
      DATA G(40)/-.666422765408D+00/IDPOT(40)/  3/ITPOT(40)/ 28/
      DATA G(41)/0.726086323499D+00/IDPOT(41)/  3/ITPOT(41)/ 28/
      DATA G(42)/0.550686686128D-01/IDPOT(42)/ 10/ITPOT(42)/ 24/
      DATA G(43)/-.100000000D+01/G(50)/-.100000000D+05/
      DATA G(44)/-.100000000D+01/G(51)/-.100000000D+05/
      DATA G(45)/-.100000000D+01/G(52)/-.100000000D+05/
      DATA G(46)/-.100000000D+01/G(53)/-.100000000D+05/
      DATA G(47)/-.100000000D+01/G(54)/-.100000000D+05/
      DATA G(48)/-.100000000D+01/G(55)/-.100000000D+05/
      DATA G(49)/0.000000000D+00/G(56)/0.000000000D+00/
      DATA ALPHA(35) /0.250000000D+02/
      DATA  BETA(35) /0.325000000D+03/
      DATA GAMMA(35) /0.116000000D+01/
      DATA DELTA(35) /0.100000000D+01/
      DATA   ETA(35) /0.000000000D+00/
      DATA   PHI(35) /0.000000000D+00/
      DATA ALPHA(36) /0.250000000D+02/
      DATA  BETA(36) /0.300000000D+03/
      DATA GAMMA(36) /0.119000000D+01/
      DATA DELTA(36) /0.100000000D+01/
      DATA   ETA(36) /0.000000000D+00/
      DATA   PHI(36) /0.000000000D+00/
      DATA ALPHA(37) /0.250000000D+02/
      DATA  BETA(37) /0.300000000D+03/
      DATA GAMMA(37) /0.119000000D+01/
      DATA DELTA(37) /0.100000000D+01/
      DATA   ETA(37) /0.000000000D+00/
      DATA   PHI(37) /0.000000000D+00/
      DATA ALPHA(38) /0.150000000D+02/
      DATA  BETA(38) /0.275000000D+03/
      DATA GAMMA(38) /0.125000000D+01/
      DATA DELTA(38) /0.100000000D+01/
      DATA   ETA(38) /0.000000000D+00/
      DATA   PHI(38) /0.000000000D+00/
      DATA ALPHA(39) /0.200000000D+02/
      DATA  BETA(39) /0.275000000D+03/
      DATA GAMMA(39) /0.122000000D+01/
      DATA DELTA(39) /0.100000000D+01/
      DATA   ETA(39) /0.000000000D+00/
      DATA   PHI(39) /0.000000000D+00/
      DATA ALPHA(40) /0.875000000D+00/
      DATA  BETA(40) /0.300000000D+00/
      DATA GAMMA(40) /0.700000000D+00/
      DATA DELTA(40) /0.100000000D+02/
      DATA   ETA(40) /0.275000000D+03/
      DATA   PHI(40) /0.000000000D+00/
      DATA ALPHA(41) /0.925000000D+00/
      DATA  BETA(41) /0.300000000D+00/
      DATA GAMMA(41) /0.700000000D+00/
      DATA DELTA(41) /0.100000000D+02/
      DATA   ETA(41) /0.275000000D+03/
      DATA   PHI(41) /0.000000000D+00/
      DATA ALPHA(42) /0.875000000D+00/
      DATA  BETA(42) /0.300000000D+00/
      DATA GAMMA(42) /0.700000000D+00/
      DATA DELTA(42) /0.125000000D+02/
      DATA   ETA(42) /0.275000000D+03/
      DATA   PHI(42) /0.000000000D+00/
      DATA TNORM /0.304128200D+03/DNORM /0.467600000D+03/
      DATA TLIM  /0.500000000D+04/PLIM/0.500000000D+04/
      DATA TNORM1/0.304128200D+03/DNORM1/0.467600000D+03/
      DATA TCEQ  /0.304128200D+03/TCEQ1 /0.304128200D+03/
      DATA PCEQ  /0.737730000D+01/PCEQ1 /0.737730000D+01/
      DATA DCEQ  /0.467599998D+03/DCEQ1 /0.467599998D+03/
      END
C*******************************************************
       DOUBLEPRECISION FUNCTION  PSCHM(T,ID)
C*******************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      DIMENSION TPOT(2),COEF(2)
      DATA TPOT /1.D0,2.D0/
      DATA COEF /1.95553898D3,2.05545926D3/
      DATA TTR  /216.592D0/
      DATA PTR  /0.51795D0/
      IDUM=ID
      IF(T.LT.TTR) THEN
       WRITE(*,100)
 100   FORMAT('PSCHM  CALCULATION NOT POSSIBLE: T < TTR')
       PSCHM=PTR
       RETURN
      ENDIF
      TN=T/TTR-1.d0
      P=1.D0
      DO 10 I=1,2
  10  P=P+COEF(I)*(TN**TPOT(I))
      PSCHM=P*PTR
      RETURN
      END
C*******************************************************
       DOUBLEPRECISION FUNCTION  PSUB(T,ID)
C*******************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      SAVE
      DIMENSION TPOT(3),COEF(3)
      DATA TPOT /1.D0,1.9D0,2.9D0/
      DATA COEF /-1.47408463D1,2.4327015D0,-5.3061778D0/
      DATA TTR  /216.592D0/
      DATA PTR  /0.51795D0/
      IDUM=ID
      IF(T.GT.TTR) THEN
       WRITE(*,100)
 100   FORMAT('PSUB  CALCULATION NOT POSSIBLE: T > TTR')
       PSUB=PTR
       RETURN
      ENDIF
      TN=1.D0-T/TTR
      P=0.D0
      DO 10 I=1,3
  10  P=P+COEF(I)*(TN**TPOT(I))
      PSUB=PTR*DEXP(P*TTR/T)
      RETURN
      END
C***********************************************************************
      subroutine TBER(P, D, EPST, T, IPHASE)
C *********************************************************************
      implicit real*8 (A-H,O-Z)
      real*8 NULLP
      external NULLP
C
      COMMON/CSUB2/RG,XMOL,TC,PC,DC
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PB,DLB,DVB
      DIMENSION Z(10)
C
      IPHASE = 1
      TS = TTR
      DL = DC
      DV = 0.D0
      if ((P.GT.PTR).AND.(P.LT.PC)) then
       call SATPEQ(TS,DV,DL,P,1.D-6)
       if ((D.GT.(0.99D0*DV)).AND.(D.LT.(1.01D0*DL))) THEN
        call SATCAP(TS, DV, DL, P, 1.D-6)
        if ((D.GT.DV).AND.(D.LT.DL)) THEN
         IPHASE = 0
         T = TS
         goto 999
        end if
       end if
       if (D.GT.DL) IPHASE = 2
      end if
C
      TEST = 0.D0
      call TBVDW(P,D,TEST)
C
      if ((TEST.LT.TTR).AND.(TEST.LT.TS)) TEST = TS
C
      IDKRIT = 0
      IPKRIT = 0
      if ((D.GT.(0.9D0*DC)).AND.(D.LT.(1.5D0*DC))) IDKRIT = 1
      if ((P.GT.(0.95D0*PC)).AND.(P.LT.(1.4D0*PC))) IPKRIT = 1
      if ((IDKRIT*IPKRIT).EQ.1) TEST = 1.06D0*TC
C
      if (TEST.NE.TS) then
       T1 = 0.95D0*TEST
       T2 = 1.05D0*TEST
      else
       T1 = TS
       T2 = 1.5D0*TC
      end if
C
      IF (EPST.LT.1.D-7) EPST = 1.D-7
C
      CALL WNULL2(T1,T2,NULLP,P,D,EPST,X,IX)
C
      IF (IX.LT.1) THEN
       T = X
       IF ((T.LT.TC).AND.(P.GT.PC)) THEN
        CALL DLEQN(T,DL)
        IF (D.LT.DL) IX = 1
       END IF
      END IF
      if (IX.GE.1) then
c      WRITE (*,1001)
       T1 = TS
       T2 = 2.D0*TC
       Z(1) = 0.D0
       IDIM = 10
       IANZ = 20
       CALL WNULL(T1,T2,NULLP,P,D,EPST,IANZ,Z,N,IDIM,IX)
       DO 98 I=1,IDIM
        IF (Z(I).GT.TC) THEN
         T = Z(I)
         GOTO 999
        END IF
        IF (Z(I).GT.TTR) THEN
         T = Z(I)
         CALL SATCAL(T,DV,DL,P,1.D-6)
         IF ((D.LT.DV).OR.(D.GT.DL)) GOTO 999
        END IF
   98  continue
c      write (*,1002) P, D
       T = 0.D0
      end if
C
c1001 format (/,' <TCAL>: WNULL2 failed, try again with WNULL')
c1002 format (/,' <TCAL>: Iteration failed for P =',
c    *        F9.4,' and D =',F9.4,' => T=0.D0',/)
C
 999  continue
      return
      end

C ************************************************************************
      subroutine TBVDW (P,D,TEST)
C ************************************************************************
      implicit real*8(A-H,O-Z)
C
      COMMON/CSUB2/RG,XMOL,TC,PC,DC
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PB,DLB,DVB
C
      A = 27.D0/64.D0 *RG*RG*TC*TC/PC
      B = RG*TC/(8.D0*PC)
      V = 1.D0/D
C
      TEST = (P+A/(V*V))*(V-B)/RG
C
      return
      end
C***********************************************************************

      doubleprecision function NULLP (TAKT,D,P)
      implicit real*8 (A-H,O-Z)
C
      PAKT = PB(TAKT,D)
      NULLP = PAKT - P
      return
      end

C******************************************************************
      SUBROUTINE PSITER (P,S,T,D,EPSS)
C******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL SNULL
      REAL*8 SNULL
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PBOIL,DLB,DVB
      COMMON/PSOLD/POLD,SOLD,TOLD,DOLD
      COMMON/CPHASE/IPHASE
C
      IF ((DABS(P-POLD).LT.1.D-6).AND.(DABS(S-SOLD).LT.1.D-6)) THEN
       T = TOLD
       D = DOLD
       RETURN
      END IF
C
      EPST = EPSS/1.D2
      IPHASE = 0
C
      IF (P.GE.PC) THEN
       CALL TBER(P, DC, EPST, T, IPHASE)
       TOLD = T
       IF (T.LT.TTR) GOTO 900
       SDC = SB(T,DC)
       IF (S.GT.SDC) TFAK = 1.2D0
       IF (S.LE.SDC) TFAK = 0.9D0
       SDOLD = S-SDC
       SOLD = SDC
C
 101   CONTINUE
       T = TOLD*TFAK
       IF (T.LT.TTR) T = TTR+1.D-3
       IF (T.EQ.TOLD) GOTO 900
       D = DB(T,P)
       SNEU = SB(T,D)
       SDNEU = S-SNEU
       IF ((SDOLD*SDNEU).LE.0.D0) GOTO 102
       SOLD = SNEU
       SDOLD = SDNEU
       TOLD = T
       GOTO 101
C
 102   TNEU = T
       CALL WNULL3(TOLD,TNEU,SNULL,S,P,EPST,T,IX)
       IF (IX.EQ.0) THEN
        D = DB(T,P)
       ELSE
        GOTO 900
       END IF
C
      ELSE IF (P.LT.PC) THEN
       IF (P.GT.PTR) THEN
        EPS = 1.D-6
        CALL SATCAP(TS,DV,DL,P,EPS)
        SV=SB(TS,DV)
        SL=SB(TS,DL)
        IF ((SL.LE.S).AND.(SV.GE.S)) THEN
         T = TS
         X = (S-SL)/(SV-SL)
         V = 1.D0/DL + X*(1.D0/DV-1.D0/DL)
         D = 1.D0/V
         RETURN
        END IF
       ELSE
        TS = TTR+0.1D0
        D = DB(TS,P)
        SV = SB(TS,D)
        SL = -1.D6
       END IF
C
       IF (S.GT.SV) TFAK = 1.2D0
       IF (S.GT.SV) SOLD = SV
       IF (S.GT.SV) IPHASE = 2
       IF (S.LT.SL) TFAK = 1.D0-(TS-TTR)/(5.D0*TS)
       IF (S.LT.SL) SOLD = SL
       IF (S.LT.SL) IPHASE = 1
       TOLD = TS
       SDOLD = S-SOLD
C
 103   CONTINUE
       T = TOLD*TFAK
       D = DB(T,P)
       SNEU = SB(T,D)
       SDNEU = S-SNEU
       IF ((SDOLD*SDNEU).LE.0.D0) GOTO 104
       SOLD = SNEU
       SDOLD = SDNEU
       TOLD = T
       GOTO 103
C
 104   TNEU = T
       CALL WNULL3(TOLD,TNEU,SNULL,S,P,EPST,T,IX)
       IF (IX.EQ.0) THEN
        D = DB(T,P)
       ELSE
        GOTO 900
       END IF
      END IF
C
      RETURN
C
 900  T = 0.D0
      D = 0.D0
C
      END
C***********************************************************************
      doubleprecision function SNULL (TAKT,P,S)
      implicit real*8 (A-H,O-Z)
      COMMON/CPHASE/IPHASE
C
      EPSD=1.D-6
      DEST = BDENS(TAKT,P,IPHASE)
      D = DENS(P,TAKT,DEST,EPSD)
      SNULL = SB(TAKT,D) - S
      return
      end
C******************************************************************
      SUBROUTINE PHITER (P,H,T,D,EPSH)
C******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL HNULL
      REAL*8 HNULL
      COMMON/CSUB2/R,XMOL,TC,PC,DC
      COMMON/CSUB3/TTR,PTR,DLTR,DVTR,TB,PBOIL,DLB,DVB
      COMMON/PHOLD/POLD,HOLD,TOLD,DOLD
      COMMON/CPHASE/IPHASE
C
      IF ((DABS(P-POLD).LT.1.D-6).AND.(DABS(H-HOLD).LT.1.D-6)) THEN
       T = TOLD
       D = DOLD
       RETURN
      END IF
C
      EPST = EPSH/1.D2
      IPHASE = 0
C
      IF (P.GE.PC) THEN
       CALL TBER(P, DC, EPST, T, IPHASE)
       TOLD = T
       IF (T.LT.TTR) GOTO 900
       HDC = HB(T,DC)
       IF (H.GT.HDC) TFAK = 1.2D0
       IF (H.LE.HDC) TFAK = 0.9D0
       HDOLD = H-HDC
       HOLD = HDC
C
 101   CONTINUE
       T = TOLD*TFAK
       IF (T.LT.TTR) T = TTR+1.D-3
       IF (T.EQ.TOLD) GOTO 900
       D = DB(T,P)
       HNEU = HB(T,D)
       HDNEU = H-HNEU
       IF ((HDOLD*HDNEU).LE.0.D0) GOTO 102
       HOLD = HNEU
       HDOLD = HDNEU
       TOLD = T
       GOTO 101
C
 102   TNEU = T
       CALL WNULL3(TOLD,TNEU,HNULL,H,P,EPST,T,IX)
       IF (IX.EQ.0) THEN
        D = DB(T,P)
       ELSE
        GOTO 900
       END IF
C
      ELSE IF (P.LT.PC) THEN
       IF (P.GT.PTR) THEN
        EPS = 1.D-6
        CALL SATCAP(TS,DV,DL,P,EPS)
        HV=HB(TS,DV)
        HL=HB(TS,DL)
        IF ((HL.LE.H).AND.(HV.GE.H)) THEN
         T = TS
         X = (H-HL)/(HV-HL)
         V = 1.D0/DL + X*(1.D0/DV-1.D0/DL)
         D = 1.D0/V
         RETURN
        END IF
       ELSE
        TS = TTR+0.1D0
        D = DB(TS,P)
        HV = HB(TS,D)
        HL = -1.D6
       END IF
C
       IF (H.GT.HV) TFAK = 1.2D0
       IF (H.GT.HV) HOLD = HV
       IF (H.GT.HV) IPHASE = 2
       IF (H.LT.HL) TFAK = 1.D0-(TS-TTR)/(5.D0*TS)
       IF (H.LT.HL) HOLD = HL
       IF (H.LT.HL) IPHASE = 1
       TOLD = TS
       HDOLD = H-HOLD
C
 103   CONTINUE
       T = TOLD*TFAK
       D = DB(T,P)
       HNEU = HB(T,D)
       HDNEU = H-HNEU
       IF ((HDOLD*HDNEU).LE.0.D0) GOTO 104
       HOLD = HNEU
       HDOLD = HDNEU
       TOLD = T
       GOTO 103
C
 104   TNEU = T
       CALL WNULL3(TOLD,TNEU,HNULL,H,P,EPST,T,IX)
       IF (IX.EQ.0) THEN
        D = DB(T,P)
       ELSE
        GOTO 900
       END IF
      END IF
C
      RETURN
C
 900  T = 0.D0
      D = 0.D0
C
      END
C***********************************************************************
      doubleprecision function HNULL (TAKT,P,H)
      implicit real*8 (A-H,O-Z)
      COMMON/CPHASE/IPHASE
C
      EPSD=1.D-6
      DEST = BDENS(TAKT,P,IPHASE)
      D = DENS(P,TAKT,DEST,EPSD)
      HNULL = HB(TAKT,D) - H
      return
      end
C***********************************************************************
      doubleprecision function TBPS(P,S)
C***********************************************************************
      implicit real*8 (A-H,O-Z)
C
      EPSS=1.D-5
      CALL PSITER (P,S,T,D,EPSS)
      TBPS = T
      return
      end
C***********************************************************************
      doubleprecision function DBPS(P,S)
C***********************************************************************
      implicit real*8 (A-H,O-Z)
C
      EPSS=1.D-5
      CALL PSITER (P,S,T,D,EPSS)
      DBPS = D
      return
      end
C***********************************************************************
      doubleprecision function HBPS(P,S)
C***********************************************************************
      implicit real*8 (A-H,O-Z)
C
      EPSS=1.D-5
      CALL PSITER (P,S,T,D,EPSS)
      HBPS = HB2(T,D)
      return
      end
C***********************************************************************
      doubleprecision function TBPH(P,H)
C***********************************************************************
      implicit real*8 (A-H,O-Z)
C
      EPSH=1.D-5
      CALL PHITER (P,H,T,D,EPSH)
      TBPH = T
      return
      end
C***********************************************************************
      doubleprecision function DBPH(P,H)
C***********************************************************************
      implicit real*8 (A-H,O-Z)
C
      EPSH=1.D-5
      CALL PHITER (P,H,T,D,EPSH)
      DBPH = D
      return
      end
C***********************************************************************
      doubleprecision function SBPH(P,H)
C***********************************************************************
      implicit real*8 (A-H,O-Z)
C
      EPSH=1.D-5
      CALL PHITER (P,H,T,D,EPSH)
      SBPH = SB2(T,D)
      return
      end
C***********************************************************************
      doubleprecision function HB2(T,D)
C***********************************************************************
      implicit real*8 (A-H,O-Z)
C
      call DLEQN(T,DL)
      call DVEQN(T,DV)
      IF (((1.01D0*DL).GT.D).AND.((0.98*DV).LT.D)) THEN
       CALL SATCAL(T,DV,DL,PS,1.D-6)
       IF ((DL.GT.D).AND.(DV.LT.D)) THEN
        HL = HB(T,DL)
        HV = HB(T,DV)
        VL = 1.D0/DL
        VV = 1.D0/DV
        V = 1.D0/D
        X = (V-VL)/(VV-VL)
        HB2 = HL + X*(HV-HL)
       ELSE
        HB2 = HB(T,D)
       END IF
      ELSE
       HB2 = HB(T,D)
      END IF
      return
      end
C***********************************************************************
      doubleprecision function SB2(T,D)
C***********************************************************************
      implicit real*8 (A-H,O-Z)
C
      call DLEQN(T,DL)
      call DVEQN(T,DV)
      IF (((1.01D0*DL).GT.D).AND.((0.98*DV).LT.D)) THEN
       CALL SATCAL(T,DV,DL,PS,1.D-6)
       IF ((DL.GT.D).AND.(DV.LT.D)) THEN
        SL = SB(T,DL)
        SV = SB(T,DV)
        VL = 1.D0/DL
        VV = 1.D0/DV
        V = 1.D0/D
        X = (V-VL)/(VV-VL)
        SB2 = SL + X*(SV-SL)
       ELSE
        SB2 = SB(T,D)
       END IF
      ELSE
       SB2 = SB(T,D)
      END IF
      return
      end
C***********************************************************************
      doubleprecision function UB2(T,D)
C***********************************************************************
      implicit real*8 (A-H,O-Z)
C
      call DLEQN(T,DL)
      call DVEQN(T,DV)
      IF (((1.01D0*DL).GT.D).AND.((0.98*DV).LT.D)) THEN
       CALL SATCAL(T,DV,DL,PS,1.D-6)
       IF ((DL.GT.D).AND.(DV.LT.D)) THEN
        UL = UB(T,DL)
        UV = UB(T,DV)
        VL = 1.D0/DL
        VV = 1.D0/DV
        V = 1.D0/D
        X = (V-VL)/(VV-VL)
        UB2 = UL + X*(UV-UL)
       ELSE
        UB2 = UB(T,D)
       END IF
      ELSE
       UB2 = UB(T,D)
      END IF
      return
      end
C***********************************************************************
      doubleprecision function PB2(T,D)
C***********************************************************************
      implicit real*8 (A-H,O-Z)
C
      call DLEQN(T,DL)
      call DVEQN(T,DV)
      IF (((1.01D0*DL).GT.D).AND.((0.98*DV).LT.D)) THEN
       CALL SATCAL(T,DV,DL,PS,1.D-6)
       IF ((DL.GT.D).AND.(DV.LT.D)) THEN
        PB2 = PS
       ELSE
        PB2 = PB(T,D)
       END IF
      ELSE
       PB2 = PB(T,D)
      END IF
      return
      end
C***********************************************************************
