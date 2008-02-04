      PROGRAM SIMUL
C  TURNING BAND SIMULATION FOR TWO DIMENSIONS
C
C  BASED ON METHODS DESCRIBED IN:
C        -  BROOKER P.I.: TWO DIMENSIONAL SIMULATION BY TURNING
C           BANDS, MATH.GEOL.,17., 81-91  (1985)
C        -  JOURNEL AND HUIJBREGTS: MINING GEOSTATISTICS,
C           ACADEMIC PRESS, (1978)
C        -  GECKINLI-YAVUZ: DISCRETE FOURIER TRANSFORMATION AND
C           ITS APPLICATIONS TO POWER SPECTRA ESTIMATION, ELSEVIER,(1983)
C
C  PROGRAM  WRITTEN  BY   ANDRAS BARDOSSY
C
C  VERSION 1.2      1.6.1988
C                            INCLUDING BLOCK SIMULATION (1.1)
C                            GENERATING SEVERAL REALIZATIONS
C
C *********************************************************************
C * FIRST MODUL
C * CONTROLS PROGRAM FLOW
C * INPUT DESCRIBED SEPARATELY
C *********************************************************************
      COMMON /FURA/ XF(300),YF(300),ZF(300),ZFS(300),NF,ATZ,SRZ
      COMMON /PARAM/ KR,HSX,HSY,NK,DMAX,NPO,VMAX
      COMMON /VGPAR/ DH(5),ITY(5),C(5),NST,COSFI(5),SINFI(5),ALAM(5)
      COMMON /SIMDAT/ NS,SALAM
      COMMON /FDAT/ PE(257),PS(257),IPREV
      DIMENSION IHOL(100)
      CHARACTER F110*60,FILNAM*18
      DOUBLE PRECISION RSEED
      OPEN(15,FILE='SIMNAM',STATUS='OLD')
      WRITE(*,194)
 194  FORMAT(/,' RANDOM NUMBER SEED  [NOT 0.00] =>')
      READ(*,*) RS
      RSEED=RS
      WRITE(*,193)
 193  FORMAT(/,' ENTER NUMBER OF BANDS ')
      READ(*,*) NLIN
111     READ(15,'(A18)',END=912) FILNAM
      WRITE(*,'(A18)') FILNAM
      OPEN(10,FILE='PARAME',STATUS='OLD')
      OPEN(13,FILE='DATEXP3',STATUS='OLD')
      OPEN(12,FILE='DATGAU2',STATUS='OLD')
      OPEN(20,FILE=FILNAM)
C     RS=0.5
      SL=DRAND(RSEED)
      READ(13,1945)(PE(I),I=1,256)
      READ(12,1945)(PS(I),I=1,256)
      CLOSE(12)
      CLOSE(13)
      OPEN(13,STATUS='SCRATCH',ACCESS='DIRECT',RECL=40)
      IPREV=0
1945  FORMAT(F10.6)
C
C Read control parameters
C
C IKS =1 Conditional simulation
C IVS =1 simulated points from file
C KR  =1 Block simulation
C NST    Number of variogram structures
C NK     Maximum number of points in kriging system
C DMAX = Search radius
C
      READ(10,1000) IKS,IVS,KR,NK,NPO,NST,DMAX,HSX,HSY
      WRITE(*,1000) IKS,IVS,KR,NK,NPO,NST,DMAX,HSX,HSY
      IF(NK.GT.40) NK=40
      IF(NST.GT.5) NST=5
      DO 1 I=1,NST
      READ(10,1010) ITY(I),DH(I),C(I)
      IF((ITY(I).GT.5).OR.(ITY(I).LT.1)) THEN
        WRITE(*,*) 'Error in Parame VG-type', ITY(I)
      ENDIF
      IF(C(I).LE.0.) THEN
        WRITE(*,*) 'Error in Parame Sill', C(I)
      ENDIF
      IF(DH(I).LE.0.) THEN
        WRITE(*,*) 'Error in Parame Range', DH(I)
      ENDIF
      READ(10,*) FI,ALAM(I)
        COSFI(I)=COS(FI)
        SINFI(I)=SIN(FI)
1     CONTINUE
      NF=0
      CILE=-900.
      VMAX=0.
      ZERO=0.
      CALL VARI(DMAX,ZERO,EE)
      VMAX=-EE
      READ(10,*) ATZ,SRZ
C
C Conditioning points
C
      IF(IKS.EQ.1) THEN
          OPEN(11,FILE='FURA',STATUS='OLD')
          NF=1
          READ(11,1190) F110
12          READ(11,F110,END=19) XF(NF),YF(NF),ZF(NF)
          ZFS(NF)=0.
          NF=NF+1
          GO TO 12
19          NF=NF-1
          CLOSE(11)
      ENDIF
      WRITE(*,*) ' Number of points = ' ,NF
C
C Unconditional simulation
C
C     NLIN=1
      CALL UNCON(NLIN,RSEED)
      IF (IKS.EQ.1) THEN
      DO 33 I=1,NS
      READ(13,REC=I) XS,YS,ZS
      X=XS
      Y=YS
      CALL CKRISZ(X,Y,CC1,CC2,SS)
      WRITE(*,2119) I,ZS,CC1,CC2
2119    FORMAT(' SIM,CCF,CCF',I4,3F8.4)
      ZSO=ZS
      ZS=ZS-CC2+CC1
c      IF(ZS.LT.-100.0) THEN
c        WRITE(*,*) XS,YS,ZS,CC1,CC2,ZSO
c        READ(*,*) TETU
c        DO 994 IU=1,NF
c          WRITE(*,*) XF(IU),YF(IU),ZF(IU),ZFS(IU)
c          READ(*,*) LOTE
c994     CONTINUE
c      ENDIF
      WRITE(13,REC=I) XS,YS,ZS
33      CONTINUE
      ENDIF
      JUJ=0
      ZTS=0.0
      XTS=0.0
      YTS=0.0
      DO 4 I=1,NS
      READ(13,REC=I) XS,YS,ZS
      IF (IKS.EQ.0) THEN
        ZS=ZS+ATZ
      ENDIF
C Block transformation
      IF(KR.EQ.1) THEN
      XTS=XTS+XS/4.
      YTS=YTS+YS/4.
      ZTS=ZTS+ZS/4.
      JUJ=JUJ+1
C
      IF(JUJ.EQ.4) THEN
        WRITE(20,2000) XTS,YTS,ZTS
        JUJ=0
        ZTS=0.0
        XTS=0.0
        YTS=0.0
      ENDIF
      ELSE
        IF(XS .LT. 0.0) XS=0.0
        IF(YS .LT. 0.0) YS=0.0
        WRITE(20,2000) XS,YS,ZS
      ENDIF
4     CONTINUE
      CLOSE(10)
      CLOSE(11)
      CLOSE(20)
      GO TO 111
912   STOP
1000  FORMAT(6I2,3F8.2)
1010  FORMAT(I2,2F8.2)
1020  FORMAT(2F8.2)
1110  FORMAT(3F10.2)
1190  FORMAT(A60)
2000  FORMAT(f20.10,f20.10,f20.8)
      END
      
      SUBROUTINE UNCON(NLIN,RSEED)
C
C *********************************************************************
C * UNCON SUBROUTINE: SIMULATION SUBROUTINE                           *
C *                  XS,YS,ZS - POINTS AND SIMULATED VALUES           *
C *                  KR       - IF > 0 BLOCK SIMULATION               *
C *                  HSX,HSY  - BLOCK SIZE                            *
C *                  DMAX     - SEARCH RADIUS                         *
C *                  NK       - MAXIMAL NO. OF POINTS                 *
C *********************************************************************
C
      COMMON /PARAM/ KR,HSX,HSY,NK,DMAX,NPO,VMAX
      COMMON /VGPAR/ DH(5),ITY(5),C(5),NST,COSFI(5),SINFI(5),ALAM(5)
      COMMON /SIMDAT/ NS,SALAM
      COMMON /FURA/ XF(300),YF(300),ZF(300),ZFS(300),NF,ATZ,SRZ
      COMMON /LOCKB/ XLONG(3000)
      DOUBLE PRECISION RSEED
      DATA PI/3.14159265359/
C
C READ POINT COORDINATES FOR SIMULATION
C
      OPEN(11,FILE='SIMKOR',STATUS='OLD')
      READ(11,*,END=9) XS,YS
      NTTS=1
      XMAX=XS+HSX
      YMAX=YS+HSY
      ZS=0.
      IF(KR.EQ.1) THEN
           X1=XS-HSX/4.
           Y1=YS-HSY/4.
           WRITE(13,REC=1) X1,Y1,ZS
           Y1=YS+HSY/4.
           WRITE(13,REC=2) X1,Y1,ZS
           X1=XS+HSX/4.
           Y1=YS-HSY/4.
           WRITE(13,REC=3) X1,Y1,ZS
           Y1=YS+HSY/4.
           WRITE(13,REC=4) X1,Y1,ZS
           NS=4
      ELSE
           WRITE(13,REC=1) XS,YS,ZS
           NS=1
      ENDIF
      XMIN=XMAX-2.*HSX
      YMIN=YMAX-2.*HSY
2       READ(11,*,END=9) XS,YS
      NTTS=1+NTTS
      IF(XS.GT.XMAX) XMAX=XS+HSX
      IF(YS.GT.YMAX) YMAX=YS+HSY
      IF(XS.LT.XMIN) XMIN=XS-HSX
      IF(YS.LT.YMIN) YMIN=YS-HSY
      ZS=0.
      IF(KR.EQ.1) THEN
           X1=XS-HSX/4.
           Y1=YS-HSY/4.
           WRITE(13,REC=1+NS) X1,Y1,ZS
           Y1=YS+HSY/4.
           WRITE(13,REC=2+NS) X1,Y1,ZS
           X1=XS+HSX/4.
           Y1=YS-HSY/4.
           WRITE(13,REC=3+NS) X1,Y1,ZS
           Y1=YS+HSY/4.
           WRITE(13,REC=4+NS) X1,Y1,ZS
           NS=4+NS
      ELSE
           WRITE(13,REC=1+NS) XS,YS,ZS
           NS=1+NS
      ENDIF
      GO TO 2
9       CONTINUE
      CLOSE(11)
C 
C Reset ALAM
C
CShe  ALAM=1.0
      DO 98 IU=1,NF
        IF(XF(IU).LT.XMIN) XMIN=XF(IU)
        IF(YF(IU).LT.YMIN) YMIN=YF(IU)
        IF(XF(IU).GT.XMAX) XMAX=XF(IU)
        IF(YF(IU).GT.YMAX) YMAX=YF(IU)
98    CONTINUE

C
C   INITIALIZE VARIABLES (ANGLE,CENTER)
C
      ANLIN=1./SQRT(FLOAT(NLIN))
      PSZI=PI/FLOAT(NLIN)
      XCENT=0.5*(XMAX+XMIN)
      YCENT=0.5*(YMAX+YMIN)
      XTAV=SQRT((XMAX-XMIN)**2+(YMAX-YMIN)**2)
      NHOS=1
      write(*,*) NS,XTAV,XCENT,YCENT
C
C   REPEAT FOR EACH VARIOGRAM AND LINE
C
      XMA1 = XMAX
      XMI1 = XMIN
      YMA1 = YMAX
      YMI1 = YMIN
      DO 101 IN=1,NST
      CIN=SQRT(C(IN))
      IF(ITY(IN).EQ.1) THEN
           DO 22 I=1,NF
           AB=DRAND(RSEED)
           CALL NDTRI(AB,XA,DA,IE)
           ZFS(I)=ZFS(I)+XA*CIN
 22          CONTINUE
           DO 23 I=1,NS
           AB=DRAND(RSEED)
           CALL NDTRI(AB,XA,DA,IE)
           READ(13,REC=I) XS,YS,ZS
           ZS=ZS+XA*CIN
           WRITE(13,REC=I) XS,YS,ZS
  23         CONTINUE
           GO TO 101
      ENDIF
      CALL CTRANS(HSX,HSY,ALAM(IN),COSFI(IN),SINFI(IN))
        DO 51 I=1,NS
         READ(13,REC=I)XS,YS,ZS
         CALL CTRANS(XS,YS,ALAM(IN),COSFI(IN),SINFI(IN))      
         WRITE(13,REC=I)XS,YS,ZS
 51     CONTINUE
        DO 52 I=1,NF
         CALL CTRANS(XF(I),YF(I),ALAM(IN),COSFI(IN),SINFI(IN))      
 52     CONTINUE
         XMAX = -9999999999.0
         XMIN = +9999999999.0
         YMAX = XMAX
         YMIN = XMIN
         XM1 = XMA1
         XM2 = XMI1
         YM1 = YMA1
         YM2 = YMI1
         CALL CTRANS(XM1,YM1,ALAM(IN),COSFI(IN),SINFI(IN))      
         IF(XMAX .LT. XM1) XMAX = XM1
         IF(XMIN .GT. XM1) XMIN = XM1
         IF(YMAX .LT. YM1) YMAX = YM1
         IF(YMIN .GT. YM1) YMIN = YM1
         XM1 = XMA1
         XM2 = XMI1
         YM1 = YMA1
         YM2 = YMI1
         CALL CTRANS(XM1,YM2,ALAM(IN),COSFI(IN),SINFI(IN))      
         IF(XMAX .LT. XM1) XMAX = XM1
         IF(XMIN .GT. XM1) XMIN = XM1
         IF(YMAX .LT. YM2) YMAX = YM2
         IF(YMIN .GT. YM2) YMIN = YM2
         XM1 = XMA1
         XM2 = XMI1
         YM1 = YMA1
         YM2 = YMI1
         CALL CTRANS(XM2,YM1,ALAM(IN),COSFI(IN),SINFI(IN))      
         IF(XMAX .LT. XM2) XMAX = XM2
         IF(XMIN .GT. XM2) XMIN = XM2
         IF(YMAX .LT. YM1) YMAX = YM1
         IF(YMIN .GT. YM1) YMIN = YM1
         XM1 = XMA1
         XM2 = XMI1
         YM1 = YMA1
         YM2 = YMI1
         CALL CTRANS(XM2,YM2,ALAM(IN),COSFI(IN),SINFI(IN))      
         IF(XMAX .LT. XM2) XMAX = XM2
         IF(XMIN .GT. XM2) XMIN = XM2
         IF(YMAX .LT. YM2) YMAX = YM2
         IF(YMIN .GT. YM2) YMIN = YM2
         XCENT = (XMAX + XMIN)/2.
         YCENT = (YMAX + YMIN)/2.
         XTAV=SQRT((XMAX-XMIN)**2+(YMAX-YMIN)**2)
     
      NTTS=1
      XMAX=XS+HSX
      YMAX=YS+HSY
      FACT=XTAV/DH(IN)
      WRITE(*,*) FACT
      CIN=SQRT(C(IN))
      DO 11 ID=1,NLIN
      UPSZ=PSZI*FLOAT(ID-1)
C
C  PERFORM ONE DIMENSIONAL SIMULATION
C
      CALL SIMLIN(ITY(IN),ALE,RSEED,FACT,NHOS)
C
C     WRITE(*,2121) IN,ID
2121  FORMAT(' AFTER SIMLIN:',2I5)
C
C  SIMULATE CONDITIONING POINTS
C
      DO 12 I=1,NF
      XT=(XF(I)-XCENT)*COS(UPSZ)+(YF(I)-YCENT)*SIN(UPSZ)
      XT=ALE*XT/DH(IN)+0.5
      LO=INT(XT)+(NHOS/2)
      IF(LO.GT.NHOS) THEN
        WRITE(*,*) 'NHOS,LO',NHOS,LO
        READ(*,*) TETU
      ENDIF
      ZFS(I)=ZFS(I)+ANLIN*XLONG(LO)*CIN
 12   CONTINUE
C
C  SIMULATE REQUIRED POINTS
C
      DO 13 I=1,NS
      READ(13,REC=I) XS,YS,ZS
      XT=(XS-XCENT)*COS(UPSZ)+(YS-YCENT)*SIN(UPSZ)
C
      XT=ALE*XT/DH(IN)+0.5
      LO=INT(XT)+(NHOS/2)
      IF(LO.GT.NHOS) THEN
        WRITE(*,*) 'NHOS,LO',NHOS,LO
        READ(*,*) TETU
      ENDIF
      ZS=ZS+ANLIN*XLONG(LO)*CIN
      WRITE(13,REC=I) XS,YS,ZS
  13  CONTINUE
11    CONTINUE
      DO 61 I=1,NS
        READ(13,REC=I)XS,YS,ZS
        CALL CRTRANS(XS,YS,ALAM(IN),COSFI(IN),SINFI(IN))
        WRITE(13,REC=I)XS,YS,ZS
 61   CONTINUE
      DO 62 I=1,NF
        CALL CRTRANS(XF(I),YF(I),ALAM(IN),COSFI(IN),SINFI(IN))
 62   CONTINUE

100    CONTINUE
101   CONTINUE
      RETURN
      END

      SUBROUTINE CKRISZ(X,Y,CCF,CCS,SS)
C
C *********************************************************************
C * CKRISZ SUBRUTINE: PERFORMS POINT KRIGING                          *
C *   USED FOR THE CONDITIONING AND THE UNCONDITIONAL DATA            *
C *                  XF,YF,ZF - CONDITIONING DATA                     *
C *                  ZS       - UNCONDITIONAL RESULT                  *
C *                  X,Y      - POINT TO BE ESTIMATED                 *
C *                  CC1      - KRIGED RESULT (CONDITIONAL)           *
C *                  CC2      - KRIGED RESULT (UNCONDITIONAL)         *
C *                  DMAX     - SEARCH RADIUS                         *
C *                  NK       - MAXIMAL NUMBER OF POINTS              *
C *********************************************************************
C
      COMMON /PARAM/ KR,HSX,HSY,NK,DMAX,NPO,VMAX
      COMMON /VGPAR/ DH(5),ITY(5),C(5),NST,COSFI(5),SINFI(5),ALAM(5)
      COMMON /FURA/ XF(300),YF(300),ZF(300),ZFS(300),NF,ATZ,SRZ
      DIMENSION A(42,42),B(42),IHOL(300),ERE(42),ZM(300)
      DIMENSION DIST(300)
      R2=DMAX*DMAX
      J=0
      DO 9 I=1,NF
      DX=X-XF(I)
      DY=Y-YF(I)
      DIS=DX**2+DY**2
      IF(DIS.GT.R2) GO TO 9
      IF(J.EQ.0) GOTO 201
      IF(J.GE.NK.AND.DIS.GT.DIST(J)) GO TO 9
  201   J=J+1
      DIST(J)=DIS
      IHOL(J)=I
      IF(J.EQ.1) GO TO 9
      CALL DSORT(DIST,IHOL,J)
      IF(J.GT.NK) J=NK
 9      CONTINUE
      IF(J.EQ.0) THEN
          write(*,916) X,Y
916       format(' Simulation error at location ',2f12.2)
          SS=0.
          CCF=0.
          CCS=0.
          read(*,*) tetu
          RETURN
      ENDIF
      N1=J-1
      DO 1 I=1,N1
      JJ=I+1
      DO 1 K=JJ,J
      DX=XF(IHOL(I))-XF(IHOL(K))
      DY=YF(IHOL(I))-YF(IHOL(K))
      CALL VARI(DX,DY,E)
      A(I,K)=E
 1      A(K,I)=E
      N2=J+1
      DO 2 L=1,N2
 2      A(L,L)=VMAX
      A(N2,N2)=0.
      DO 3 L=1,J
      A(N2,L)=1.
      A(L,N2)=1.
      DX=X-XF(IHOL(L))
      DY=Y-YF(IHOL(L))
      CALL VARI(DX,DY,E)
      B(L)=E
 3      A(L,N2+1)=E
      MN=N2+1
      A(N2,MN)=1.
      CALL FINDT(A,ERE,N2)
      SS=0.
      CCF=0.
      CCS=0.
      DO 6 I=1,J
      SS=SS+ERE(I)*B(I)
      CCS=CCS+ERE(I)*ZFS(IHOL(I))
6     CCF=CCF+ERE(I)*ZF(IHOL(I))
      RETURN
      END
      
      SUBROUTINE VARI(X,Y,ERT)
C *********************************************************************
C * THEORETICAL VARIOGRAM SUBROUTINE
C * 1 - NUGGET EFFECT  2 - SPHERICAL 3 - EXPONENCIAL 4 - POLINOMIAL
C * 5 - GAUSSIAN
C *********************************************************************
      COMMON /VGPAR/ D(5),ITY(5),C(5),NST,COSFI(5),SINFI(5),ALAM(5)
      COMMON /PARAM/ KR,HSX,HSY,NK,DMAX,NPO,VMAX
      ERT=VMAX

      DO 1 I=1,NST
      DX=COSFI(I)*X+SINFI(I)*Y
      DY=-SINFI(I)*X+COSFI(I)*Y
      DT=SQRT(DX*DX*ALAM(I)*ALAM(I)+DY*DY)
      GOTO (10,20,30,40,50) ITY(I)
10    ERT=ERT-C(I)
      GO TO 1
20    IF(DT.GT.D(I)) THEN
          ERT=ERT-C(I)
      ELSE
          ERT=ERT-C(I)*(1.5*DT/D(I)-0.5*(DT/D(I))**3)
      ENDIF
      GOTO 1
30    ERT=ERT-C(I)+C(I)*EXP(-DT/D(I))
      GO TO 1
40    ERT=ERT-C(I)*DT**D(I)
      GO TO 1
50    ERT=ERT-(1.-EXP(-DT**2/D(I)**2))*C(I)
      GO TO 1
1     CONTINUE
      RETURN
      END

      SUBROUTINE DSORT(DS,IHOL,N)
C ***********************************************************************
C * DSORT SUBROUTINE: SORTS THE ARRAY DS IN ASCENDING ORDER             *
C ***********************************************************************
      DIMENSION DS(1),IHOL(1)
      DNEW=DS(N)
      NEW=IHOL(N)
      N1=N-1
      DO 20 I=1,N1
      K=I
      IF(DNEW.LT.DS(I)) GO TO 30
 20     CONTINUE
      RETURN
 30     JK=0
      DO 40 I=K,N1
      J=N1-JK
      JK=JK+1
      DS(J+1)=DS(J)
      IHOL(J+1)=IHOL(J)
 40     CONTINUE
      DS(K)=DNEW
      IHOL(K)=NEW
 50     RETURN
      END
      
      SUBROUTINE FINDT(A,X,N)
C *********************************************************************
C * FINDT SZUBRUTIN
C * SOLVES THE KRIGING EQUATIONS - IN THE COVARIANCE FORM
C * THE DIAGONAL DOMINATES THE MATRIX
C *********************************************************************
      DIMENSION A(42,42),X(42)
      MP=N+1
      DO 10 I=1,N
      IP=I+1
      DO 10 J=1,N
      IF(I-J) 6,10,6
6     F=(-A(J,I))/A(I,I)
      DO 9 K=IP,MP
9     A(J,K)=A(J,K)+F*A(I,K)
10    CONTINUE
      DO 20 I=1,N
20    X(I)=A(I,N+1)/A(I,I)
      RETURN
      END

      SUBROUTINE FWEIGHT(X,Y,N,L2N,IVA)
C *********************************************************************
C * WEIGHT CALCULATES THE F(I) VALUES TO BE USED FOR EACH BAND.
C *********************************************************************
      COMPLEX H(1025)
      REAL X(1025),Z(1025),Y(1025),ZDI(5)
      SQ2=SQRT(2.)
      DO 1 I=1,N
1       H(I)=CMPLX(X(I),0.)
C
C CALL FOURIER TRANSFORM FOR 1-DIM VARIOGRAM
C
      CALL FFT(H,N,L2N)
      DO 2 I=1,N
2       Z(I)=REAL(H(I))
      DO 9 I=1,N
9       Y(I)=AIMAG(H(I))
      DO 3 I=1,N
      IF(Z(I).LT.0.) Z(I)=ABS(Z(I))
      TI=1.
      IF(I.GT.(N/2)) TI=-1.
3       H(I)=CMPLX(0.,TI*SQRT(Z(I)))
C
C INVERSE FOURIER TRANSFORM TO OBTAIN WEIGHT VECTOR
C
      CALL FFT(H,N,L2N)
      AN=N
      DO 4 I=1,N
4       Y(I)=SQ2*REAL(H(I))/AN
C
C
C
      IF(IVA.NE.1) RETURN
      N2=N/2
      DO 5 I=1,N2
      ZA=Y(I)
      Y(I)=Y(I+N2)
5       Y(I+N2)=ZA
      RETURN
      END

      SUBROUTINE FFT(X,N,L2N)
C *********************************************************************
C
C   FFT SUBROUTINE TO COMPUTE DFT OF AN N TERM COMPLEX SEQUENCE
C   FROM THE COMPLEX ARRAY X(N), N=2**L2N
C   RESULT REPLACES X(N)
C
C *********************************************************************
      COMPLEX X(N),W,B
      DATA PI/3.14159265359/
      NV2=N/2
      NM1=N-1
      L=1
      DO 3 K=1,NM1
      IF(K.GE.L) GO TO 1
      B=X(L)
      X(L)=X(K)
      X(K)=B
1       M=NV2
2       IF(M.GE.L) GO TO 3
      L=L-M
      M=M/2
      GO TO 2
3       L=L+M
      K=1
      DO 5 L=1,L2N
      M=K
      BM=M
      K=K*2
      W=(1.0,0.0)
      DO 5 J=1,M
      BJ=J
      DO 4 I=J,N,K
      I2=I+M
      B=X(I2)*W
      X(I2)=X(I)-B
4       X(I)=X(I)+B
      ARG=PI*BJ/BM
5       W=CMPLX(COS(ARG),-SIN(ARG))
      RETURN
      END
      FUNCTION DRAND(IX)
C *********************************************************************
C
C PORTABLE RANDOM NUMBER GENERATOR
C  REFERENCE: SCRAGE,L.,"A MORE PORTABLE FORTRAN RANDOM
C                        NUMBER GENERATOR."
C        IN  ACM TRANS. ON MATH. SOFTWARE, V.5.,N.2.,1979
C
C  IX = IX * A MOD P  RECURSION USED
C
C *********************************************************************
      DOUBLE PRECISION A,P,IX,B15,B16,XHI,XALO,LEFTLO,FHI,K
      A    = 16807.D0
      B15  = 32768.D0
      B16  = 65536.D0
      P    = 2147483647.D0
C
C GET 15 HIGH ORDER BITS OF IX
C
      XHI= IX/B16
      XHI= XHI-DMOD(XHI,1.D0)
C
C GET 16 LOWER BITS OF IX AND LOW PRODUCT
C
      XALO = (IX-XHI*B16)*A
C
C GET 15 HIGH ORDER BITS OF LOW PRODUCT
C
      LEFTLO = XALO/B16
      LEFTLO = LEFTLO - DMOD(LEFTLO,1.D0)
C
C FORM THE 31 HIGHEST BITS OF FULL PRODUCT
C
      FHI = XHI*A+LEFTLO
C
C GET OVERFLOW PAST 31ST BIT OF FULL PRODUCT
C
      K = FHI/B15
      K = K-DMOD(K,1.D0)
C
C ASSEMBLE ALL THE PARTS AND PRESUBTRACT P. THE PARANTHESES ARE
C ESSENTIAL
C
      IX = (((XALO-LEFTLO*B16)-P)+(FHI-K*B15)*B16)+K
C
C ADD P BACK IF NECESSARY
C
      IF(IX.LT.0.D0) IX=IX+P
C
C MULTIPLY BY (1/(2**31-1))
C
      DRAND = IX*4.656612875D-10
      RETURN
      END

      SUBROUTINE SIMLIN(IT,ALE,RSEED,FACT,NHOS)
C *********************************************************************
C
C
C *********************************************************************
      COMMON /LOCKB/ XLONG(3000)
      COMMON /LOCKE/ SLONG(3105)
      COMMON /FDAT/ PE(257),PS(257),IPREV
      COMMON /XYBLOCK/ X(1025),Y(1025)
      DOUBLE PRECISION RSEED
      XA=0.
      DA=0.
      IE=0
      NFL=256
      IF(IT.EQ.2) THEN
        ALE=64.
      ENDIF
      IF(IT.EQ.3) THEN
        ALE=20.
      ENDIF
      IF(IT.EQ.5) THEN
        ALE=34.
      ENDIF
      NHOS=INT(ALE*(FACT+0.5))
      NHOS2=NHOS+260
      WRITE(*,*) NHOS,NHOS2
      DO 1 I=1,NHOS2
      ABE=DRAND(RSEED)
      CALL NDTRI(ABE,XA,DA,IE)
      SLONG(I)=XA
1     CONTINUE
      IF(IT.EQ.IPREV) GO TO 79
C
C  NUGGET EFFECT
C
      IF(IT.EQ.1) THEN
          DO 2 I=1,2500
2           XLONG(I)=SLONG(I)
          ALE=-1.
          RETURN
      ENDIF
C
C SPHERICAL VARIOGRAM
C
      IF(IT.EQ.2) THEN
        A=64.
        NFL=256
        L2N=8
        IVA=1
        DO 3 I=1,NFL
        T=FLOAT(I)-1.
        IF(T.GT.A) THEN
             AL=1.-1.5*((T/A)-0.5*(T/A)**3)*ASIN(A/T)
             AL=AL-(0.75*(T/A)**2)*SQRT((1-(A/T)**2))
        ELSE
             AL=1.-(3.*3.14159265/4.)*(T/A-0.5*(T/A)**3)
        ENDIF
3         X(I)=AL
        X(1)=0.5
        ALE=64.
      ENDIF
C
C        EXPONENTIAL VARIOGRAM
C
      IF(IT.EQ.3) THEN
        ALE=20.
        NFL=256
        L2N=8
        IVA=1
        DO 31 I=1,256
31        X(I)=PE(I)
      ENDIF
C
C      GAUSSIAN VARIOGRAM
C
      IF(IT.EQ.5) THEN
        ALE=34.
        NFL=256
        L2N=8
        IVA=1
        DO 32 I=1,256
32        X(I)=PS(I)
      ENDIF
C
C    LINEAR VARIOGRAM NOT INCLUDED - RESULT : STOP
C
      IF(IT.EQ.4) THEN
      STOP
      ENDIF
      CALL FWEIGHT(X,Y,NFL,L2N,IVA)
C     WRITE(*,3121)
3121  FORMAT(' AFTER FWEIGHT')
      IPREV=IT
79    DO 4 I=1,NHOS2-NFL
      XLONG(I)=0.
      DO 5 J=1,NFL
5     XLONG(I)=XLONG(I)+Y(J)*SLONG(I+J-1)
4     CONTINUE
      RETURN
      END

C********************************************************
C
      SUBROUTINE NDTRI(P,X,D,IE)
C
C********************************************************
C       GIVEN THE CUMUALITIVE PROBABILITY P, SUBROUTINE
C       NDTRI CALCULATES THE CORRESPONDING "Z" VALUE OF THE
C       NORMAL DISTRIBUTION N(0,1).
C       THIS SUBROUTINE IS FROM THE IBM SSP VERSION III.
C********************************************************
      BIG=6.0
      IE=0
      X=BIG
      D=X
      IF(P.LT.0.0.OR.P.GT.1.0) GO TO 50
       IF(P.GT.0) GO TO 10
          X=-BIG
            D=0.0
          GO TO 40
10         D=P
       IF(D.LE.0.5) GO TO 20
          D=1.0-D
20         T2=ALOG(1.0/(D*D))
       T=SQRT(T2)
       X=T-(2.515517+0.802853*T+0.010328*T2)/(1.0+1.432788*T+
     1     0.189269*T2+0.001308*T*T2)
       IF(P.GE.0.5) GO TO 30
            X=-X
30         D=0.3989423*EXP(-X*X/2.0)
40         RETURN
50      IE=-1
      RETURN
      END
      
      SUBROUTINE CTRANS(XF,YF,ALAM,COSFI,SINFI)
C ***********************************************************************
C
C ***********************************************************************      
      DX=COSFI*XF+SINFI*YF
      DY=-SINFI*XF+COSFI*YF
      XF=DX*ALAM
      YF=DY
      RETURN
      END
      
      SUBROUTINE CRTRANS(XF,YF,ALAM,COSFI,SINFI)
C ***********************************************************************
C
C ***********************************************************************      
      XF=XF/ALAM
      DX=COSFI*XF-SINFI*YF
      DY=SINFI*XF+COSFI*YF
      XF=DX
      YF=DY
      RETURN
      END
