      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 XI
      PARAMETER (NZ=1024)
      PARAMETER (NSTEP=10000)
      PARAMETER (LANMIN=10)
      PARAMETER (LANMAX=100)
      DIMENSION ZZ(NZ)
      DIMENSION PSIRE(NZ)
      DIMENSION PSIIM(NZ)
      DIMENSION VV(NZ)
      DIMENSION VIM(NZ)
      DIMENSION AKZ(NZ)
      EVEPS=0.9648533822D0
      PI=4.0D0*DATAN(1.0D0)
      HBAR=0.06350781278D0
CCCCCC
      EKIN=0.14D0
      EKIN=EKIN*EVEPS
      WRITE(6,*)'EKIN',EKIN/EVEPS
      AMH=1.008D0
CCCCCC
      ZMIN=-6.0D0
      ZMAX=10.0D0
      DZ=(ZMAX-ZMIN)/DFLOAT(NZ)
      VOLEM=DZ
CCCCCC
      DO I=1,NZ
        ZZ(I)=ZMIN+(I-1)*DZ
      ENDDO
CCCCCC
      SNZ=0.0D0
      DO I=1,NZ
        ZZ1=ZZ(I)
        CALL TRANSZ(EKIN,AMH,ZZ1,XI)
        CALL POT(ZZ1,EN)
        VV(I)=EN
        SNZ=SNZ+XI*CONJG(XI)*DZ
        PSIRE(I)=DBLE(XI)
        PSIIM(I)=DIMAG(XI)
        WRITE(6,'(I4,3F16.8)')I,ZZ(I),ABS(XI),VV(I)
      ENDDO
      WRITE(6,*)'SNZ=',SNZ
      CALL VPOTIM(ZZ,VIM)
      CALL DEFAK(DZ,AKZ)
      DDT=0.01D0
      ISTEP=0
      IST=0
      T=0.0D0
      LAITER=30
      NLAN=LAITER
      EPSL1=1.0D-21
      EPSL2=1.0D-20
      NLDEL=2
      AFLUX=0.0D0
   50 ISTEP=ISTEP+1
      T=T+DDT
      MMM=NLAN
      CALL LANCZ(PSIRE,PSIIM,DDT,VOLEM,RAV,MMM,ISTEP,AKZ,VV)
      RAV=RAV*VOLEM
        IF(RAV.LT.EPSL1) NLAN=NLAN-NLDEL
        IF(RAV.GT.EPSL2) NLAN=NLAN+NLDEL
        IF(NLAN.GT.LANMAX) NLAN=LANMAX
        IF(NLAN.LT.LANMIN) NLAN=LANMIN
      WRITE(1,*)ISTEP,RAV,NLAN
      CALL ENERGY(AMH,AKZ,PSIRE,PSIIM,VV,VAL1,VAL2,VAL)
        ENER1=VAL1*VOLEM/EVEPS
        ENER2=VAL2*VOLEM/EVEPS
        ENER=VAL*VOLEM/EVEPS
      WRITE(2,'(4F12.6)')T,ENER1,ENER2,ENER
      SS=0.0D0
      DO I=1,NZ
        SS=SS+(PSIRE(I)**2+PSIIM(I)**2)*VOLEM
      ENDDO
      IF(MOD(ISTEP,10).EQ.0)THEN
      IST=IST+1
      DO I=1,NZ
        PSI=PSIRE(I)*PSIRE(I)+PSIIM(I)*PSIIM(I)
        WRITE(99+IST,'(1X,3F12.6)')T,ZZ(I),PSI
      ENDDO
      ENDIF
      CALL FLUX(PSIRE,PSIIM,DZ,AMH,VALU)
      AFLUX=AFLUX+VALU*DDT
      WRITE(3,'(F8.4,3F12.6)')T,SS,ENER/SS,AFLUX
C      DO I=1,NZ
C        PSIRE(I)=PSIRE(I)*EXP(-VIM(I)*DDT/HBAR)
C        PSIIM(I)=PSIIM(I)*EXP(-VIM(I)*DDT/HBAR)
C      ENDDO
      IF (ISTEP.EQ.NSTEP) GO TO 40
      GO TO 50
  40  CONTINUE
      STOP
      END
CCCCCC
      SUBROUTINE TRANSZ(EKIN,AMH,ZZ,XI)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 ZI,XI
      ZI=DCMPLX(0.0D0,1.0D0)
      PI=4.0D0*DATAN(1.0D0)
      HBAR=0.06350781278D0
      AK0=-SQRT(2.0D0*AMH*EKIN)/HBAR
      ZZM=5.0D0
      DELTAZ=1.0D0
      AN=1.0D0/(2.0D0*PI*DELTAZ*DELTAZ)
      AN=AN**0.25D0
      QQ=(ZZ-ZZM)*(ZZ-ZZM)/4.0D0/DELTAZ/DELTAZ
      XI=AN*EXP(-ZI*AK0*ZZ-QQ)
      RETURN
      END
!***********************************************************
!                                                          !
!     MODEL POTENTIAL                                     !
!                                                          !
!***********************************************************
      SUBROUTINE POT(ZZ,EN)
      IMPLICIT REAL*8 (A-H,O-Z)
      EVEPS=0.9648533822D0
      V0=0.2D0*EVEPS
      SIGMA=1.0D0
      EN=V0*EXP(-0.5D0*ZZ*ZZ/SIGMA/SIGMA)
      RETURN
      END
!***********************************************************
!                                                          *
!     ABSORBING POTENTIAL                                  *
!                                                          *
!***********************************************************
      SUBROUTINE VPOTIM(ZZ,VIM)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NZ=1024)
      DIMENSION ZZ(NZ)
      DIMENSION VIM(NZ)
      EVEPS=0.9648533822D0
      V0I=0.2D0*EVEPS
      ZZA=-4.0D0
      ZZB=-3.0D0
      DO I=1,NZ
        VIM(I)=0.0D0
      ENDDO
      DO I=1,NZ
        ZZ1=ZZ(I)
        IF(ZZ1.LE.ZZB)THEN
          FAC=(ZZ1-ZZA)/(ZZB-ZZA)
          VIM(I)=V0I*FAC
        ENDIF
C      WRITE(2,*)I,ZZ(I),VIM(I)
      ENDDO
      RETURN
      END
CCCCCC
      SUBROUTINE ENERGY(AMH,AKZ,PSIRE,PSIIM,VV,VAL1,VAL2,VAL)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NZ=1024)
      DIMENSION AKZ(NZ)
      DIMENSION PSIRE(NZ)
      DIMENSION PSIIM(NZ)
      DIMENSION VV(NZ)
      DIMENSION APSI(2*NZ)
      DIMENSION HPSIZ(2*NZ)
      DIMENSION NN(1)
CCCCCC
      HBAR=0.06350781278D0
      RMST=0.5D0*HBAR*HBAR/AMH
      SN=1.0D0/DSQRT(DFLOAT(NZ))
      DO II=1,NZ
         APSI(2*II-1)=PSIRE(II)
         APSI(2*II)=PSIIM(II)
      ENDDO
      NN(1)=NZ
      CALL FOURN(APSI,NN,1,1)
CCCCCC
      DO II=1,NZ
        HPSIZ(2*II-1)=APSI(2*II-1)*AKZ(II)*AKZ(II)*SN
        HPSIZ(2*II)=APSI(2*II)*AKZ(II)*AKZ(II)*SN
      ENDDO
CC
      CALL FOURN(HPSIZ,NN,1,-1)
CC
      VAL1=0.0D0
      VAL2=0.0D0
      DO II=1,NZ
        HPSIRE1=RMST*HPSIZ(2*II-1)*SN
        HPSIRE2=VV(II)*PSIRE(II)
CC
        HPSIIM1=RMST*HPSIZ(2*II)*SN
        HPSIIM2=VV(II)*PSIIM(II)
        VAL1=VAL1+HPSIRE1*PSIRE(II)+HPSIIM1*PSIIM(II)
        VAL2=VAL2+HPSIRE2*PSIRE(II)+HPSIIM2*PSIIM(II)
      ENDDO
      VAL=VAL1+VAL2
      RETURN
      END
CCCCCC
      SUBROUTINE FOURN(DATA,NN,NDIM,ISIGN)
      INTEGER ISIGN,NDIM,NN(NDIM)
C      REAL DATA(*)
      DOUBLE PRECISION DATA(*)
      INTEGER I1,I2,I2REV,I3,I3REV,IBIT,IDIM,IFP1,IFP2,IP1,IP2,IP3,K1,
     *K2,N,NPREV,NREM,NTOT
C      REAL TEMPI,TEMPR
      DOUBLE PRECISION TEMPI,TEMPR
      DOUBLE PRECISION THETA,WI,WPI,WPR,WR,WTEMP
      NTOT=1
      DO IDIM=1,NDIM
         NTOT=NTOT*NN(IDIM)
      ENDDO
      NPREV=1
      DO IDIM=1,NDIM
         N=NN(IDIM)
         NREM=NTOT/(N*NPREV)
         IP1=2*NPREV
         IP2=IP1*N
         IP3=IP2*NREM
         I2REV=1
         DO I2=1,IP2,IP1
            IF (I2.LT.I2REV) THEN
               DO I1=I2,I2+IP1-2,2
                  DO I3=I1,IP3,IP2
                     I3REV=I2REV+I3-I2
                     TEMPR=DATA(I3)
                     TEMPI=DATA(I3+1)
                     DATA(I3)=DATA(I3REV)
                     DATA(I3+1)=DATA(I3REV+1)
                     DATA(I3REV)=TEMPR
                     DATA(I3REV+1)=TEMPI
                  ENDDO
               ENDDO
            ENDIF
            IBIT=IP2/2
   1        IF ((IBIT.GE.IP1).AND.(I2REV.GT.IBIT)) THEN
               I2REV=I2REV-IBIT
               IBIT=IBIT/2
               GOTO 1
            ENDIF
            I2REV=I2REV+IBIT
         ENDDO
         IFP1=IP1
   2     IF (IFP1.LT.IP2) THEN
            IFP2=2*IFP1
            THETA=ISIGN*6.28318530717959D0/(IFP2/IP1)
            WPR=-2.0D0*SIN(0.5D0*THETA)**2
            WPI=SIN(THETA)
            WR=1.0D0
            WI=0.0D0
            DO I3=1,IFP1,IP1
               DO I1=I3,I3+IP1-2,2
                  DO I2=I1,IP3,IFP2
                     K1=I2
                     K2=K1+IFP1
                     TEMPR=SNGL(WR)*DATA(K2)-SNGL(WI)*DATA(K2+1)
                     TEMPI=SNGL(WR)*DATA(K2+1)+SNGL(WI)*DATA(K2)
                     DATA(K2)=DATA(K1)-TEMPR
                     DATA(K2+1)=DATA(K1+1)-TEMPI
                     DATA(K1)=DATA(K1)+TEMPR
                     DATA(K1+1)=DATA(K1+1)+TEMPI
                  ENDDO
               ENDDO
               WTEMP=WR
               WR=WR*WPR-WI*WPI+WR
               WI=WI*WPR+WTEMP*WPI+WI
            ENDDO
            IFP1=IFP2
            GOTO 2
         ENDIF
         NPREV=N*NPREV
      ENDDO
      RETURN
      END
CCCCCC
      SUBROUTINE DEFAK(DZ,AKZ)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NZ=1024)
      DIMENSION AKZ(NZ)
      PI=4.0D0*DATAN(1.0D0)
      PI2=PI*2.0D0
      ACONSZ=PI2/(NZ*DZ)
      NHALFZ = NZ/2+1
      DO I=1,NHALFZ
        AKZ(I)=ACONSZ*(I-1)
      ENDDO
      DO I=NHALFZ+1,NZ
        AKZ(I)=ACONSZ*(I-1-NZ)
      ENDDO
      RETURN
      END
CCCCCC
      SUBROUTINE FLUX(PSIRE,PSIIM,VOLEM,AMH,VAL)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NZ=1024)
      DIMENSION PSIRE(NZ),PSIIM(NZ)
      HBAR=0.06350781278D0
      MD=199
      MDP1=MD-1
      MDP2=MD-2
      MDM1=MD+1
      MDM2=MD+2
      PSIREDR=(-PSIRE(MDM2)+8.0D0*PSIRE(MDM1)
     1       -8.0D0*PSIRE(MDP1)+PSIRE(MDP2))/12.0D0/VOLEM
      PSIIMDR=(-PSIIM(MDM2)+8.0D0*PSIIM(MDM1)
     1       -8.0D0*PSIIM(MDP1)+PSIIM(MDP2))/12.0D0/VOLEM
      VAL=(HBAR/AMH)*(PSIRE(MD)*PSIIMDR-PSIIM(MD)*PSIREDR)
      RETURN
      END
!***********************************************************
!                                                          *
!     HAMILTONIAN OPERATES ON THE WAVE FUNCTION            *
!     BY THE FFT-METHOD                                    *
!                                                          *
!***********************************************************
      SUBROUTINE XPRO(X0,Y0,AKZ,VV)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NZ=1024)
      DIMENSION AKZ(NZ)
      DIMENSION VV(NZ)
      DIMENSION X0(NZ),Y0(NZ)
      DIMENSION APSI(2*NZ)
      DIMENSION HPSIR(2*NZ)
      DIMENSION HPSIRE(NZ)
      DIMENSION HPSIIM(NZ)
      DIMENSION NN(1)
CCCCCC
      AMH=1.008D0
      HBAR=0.06350781278D0
      RMST=0.5D0*HBAR*HBAR/AMH
      RINV=1.0D0/HBAR
      SN=1.0D0/DSQRT(DFLOAT(NZ))
      DO I=1,NZ
         APSI(2*I-1)=X0(I)
         APSI(2*I)=Y0(I)
      ENDDO
      NN(1)=NZ
      CALL FOURN(APSI,NN,1,1)
      DO I=1,NZ
        HPSIR(2*I-1)=APSI(2*I-1)*AKZ(I)*AKZ(I)*SN
        HPSIR(2*I)=APSI(2*I)*AKZ(I)*AKZ(I)*SN
      ENDDO
      CALL FOURN(HPSIR,NN,1,-1)
      DO I=1,NZ
        HPSIRE(I)=RMST*HPSIR(2*I-1)*SN
     1           +VV(I)*X0(I)
        HPSIIM(I)=RMST*HPSIR(2*I)*SN
     1           +VV(I)*Y0(I)
      ENDDO
      DO I=1,NZ
        X0(I)=HPSIRE(I)*RINV
        Y0(I)=HPSIIM(I)*RINV
      ENDDO
      RETURN
      END
CCCCCC
      SUBROUTINE LANCZ(X1,Y1,DT,VOLEM,RAV,M,ISTEP,AKZ,VV)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NZ=1024)
      PARAMETER (LAD=100)
      DIMENSION X1(NZ),Y1(NZ)
      DIMENSION U0A(NZ),U0B(NZ)
      DIMENSION U1A(NZ),U1B(NZ)
      DIMENSION WA(NZ,LAD),WB(NZ,LAD)
      DIMENSION TT(LAD,2),D(LAD),E(LAD),Z(LAD,LAD),RA(LAD),RB(LAD)
      DIMENSION AKZ(NZ)
      DIMENSION VV(NZ)
      HBAR=0.06350781278D0
      SCALEA=DSQRT(VOLEM)
      SCALEI=1.0D0/SCALEA
      DO I=1,NZ
         U0A(I)=X1(I)*SCALEA
         U0B(I)=Y1(I)*SCALEA
         U1A(I)=0.0D0
         U1B(I)=0.0D0
      ENDDO
      CALL CHAIN(U0A,U0B,U1A,U1B,WA,WB,TT,M,AKZ,VV)
C     IF (M.EQ.1) RETURN
      IF (M.EQ.1) GO TO 1
      DO I=1,NZ
         X1(I)=0.0D0
         Y1(I)=0.0D0
      ENDDO
      DO I=2,M
         D(I)=TT(I,1)
         E(I)=TT(I-1,2)
      ENDDO
      D(1)=TT(1,1)
C     WRITE(6,1000) (D(I),I=1,M)
C     WRITE(6,1000) (E(I),I=1,M)
      DO I=1,M
         DO J=1,M
            Z(I,J)=0.0D0
            IF (I.EQ.J) Z(I,J)=1.0D0
         ENDDO
      ENDDO
C     DO 1002 I=1,M
C     PRINT 1000,I,D(I),E(I)
C1002  CONTINUE
1000  FORMAT(2X,2I5,2X,F15.5)
      CALL TQLI(D,E,M,Z)
      DO I=1,M
         SUA=0.0D0
         SUB=0.0D0
         DO J=1,M
            SUA=SUA+Z(I,J)*Z(1,J)*DCOS(D(J)*DT)
            SUB=SUB-Z(I,J)*Z(1,J)*DSIN(D(J)*DT)
         ENDDO
         RA(I)=SUA*SCALEI
         RB(I)=SUB*SCALEI
CCCCCCC  T E S T  CCCCCCCCC
C           SQRVEC=RA(I)*RA(I)+RB(I)*RB(I)
C           PRINT 100,I,RA(I),SQRVEC
C100        FORMAT(' ',I6,2E15.6)
CCCCCCCCCCCCCCCCCCCCCCCCCCC
      ENDDO
CCCCCCC  AVERAGE OF THE FIVE LAST COMPONENTS  CCCCCCC
      RAV=0.0D0
      DO I=M-4,M
         RAV=RAV+RA(I)*RA(I)+RB(I)*RB(I)
      ENDDO
      RAV=RAV*0.2D0
C     WRITE(6,1000) TA,(BR(I),I=1,N)
      DO J=1,M
         DO I=1,NZ
            X1(I)=X1(I)+WA(I,J)*RA(J)-WB(I,J)*RB(J)
            Y1(I)=Y1(I)+WA(I,J)*RB(J)+WB(I,J)*RA(J)
         ENDDO
      ENDDO
C
    1 CONTINUE
      RETURN
      END
CCCCCC
      SUBROUTINE CHAIN(U0A,U0B,U1A,U1B,WA,WB,T,M,AKZ,VV)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NZ=1024)
      PARAMETER (LAD=100)
      DIMENSION U0A(NZ),U0B(NZ)
      DIMENSION U1A(NZ),U1B(NZ)
      DIMENSION X0(NZ),Y0(NZ)
      DIMENSION WA(NZ,LAD),WB(NZ,LAD)
      DIMENSION T(LAD,2)
      DIMENSION AKZ(NZ)
      DIMENSION VV(NZ)
      BN=0.0D0
      IT=0
   10 IT=IT+1
      DO I=1,NZ
         WA(I,IT)=U0A(I)
         WB(I,IT)=U0B(I)
         X0(I)=U0A(I)
         Y0(I)=U0B(I)
      ENDDO
      CALL XPRO(X0,Y0,AKZ,VV)
C     IF (M.EQ.1) RETURN
      IF (M.EQ.1) GO TO 1
      AN=0.0D0
      DO 3 I=1,NZ
    3 AN=AN+U0A(I)*X0(I)+U0B(I)*Y0(I)
      T(IT,1)=AN
      DO I=1,NZ
         X0(I)=X0(I)-AN*U0A(I)-BN*U1A(I)
         Y0(I)=Y0(I)-AN*U0B(I)-BN*U1B(I)
      ENDDO
      BN=0.0D0
      DO I=1,NZ
         BN=BN+X0(I)*X0(I)+Y0(I)*Y0(I)
      ENDDO
      BN=DSQRT(BN)
      BN1=1.0D0/BN
      DO I=1,NZ
         U1A(I)=U0A(I)
         U1B(I)=U0B(I)
         U0A(I)=X0(I)*BN1
         U0B(I)=Y0(I)*BN1
      ENDDO
      T(IT,2)=BN
      IF (DABS(BN).LT.1.D-6) WRITE(6,289) IT,BN
  289 FORMAT(' BN= ',I4,E10.3)
      IF (IT.LT.M) GO TO 10
    1 CONTINUE
      RETURN
      END
CCCCCC
      SUBROUTINE TQLI(D,E,N,Z)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (LAD=100)
      DIMENSION D(LAD),E(LAD),Z(LAD,LAD)
      IF (N.GT.1) THEN
         DO I=2,N
            E(I-1)=E(I)
         ENDDO
         E(N)=0.D0
         DO L=1,N
            ITER=0
 1          DO M=L,N-1
               DD=DABS(D(M))+DABS(D(M+1))
               IF (DABS(E(M))+DD.EQ.DD) GO TO 2
            ENDDO
            M=N
 2          IF (M.NE.L) THEN
               IF (ITER.EQ.30) THEN
                  PRINT *,' TOO MANY ITERATIONS'
                  STOP
               ENDIF
               ITER=ITER+1
               G=(D(L+1)-D(L))/(2.D0*E(L))
               R=DSQRT(G**2+1.D0)
               G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
               S=1.0D0
               C=1.0D0
               P=0.0D0
               DO I=M-1,L,-1
                  F=S*E(I)
                  B=C*E(I)
                  IF (DABS(F).GE.DABS(G)) THEN
                     C=G/F
                     R=DSQRT(C**2+1.D0)
                     E(I+1)=F*R
                     S=1.0D0/R
                     C=C*S
                  ELSE
                     S=F/G
                     R=DSQRT(S**2+1.D0)
                     E(I+1)=G*R
                     C=1.0D0/R
                     S=S*C
                  ENDIF
                  G=D(I+1)-P
                  R=(D(I)-G)*S+2.D0*C*B
                  P=S*R
                  D(I+1)=G+P
                  G=C*R-B
                  DO K=1,N
                     F=Z(K,I+1)
                     Z(K,I+1)=S*Z(K,I)+C*F
                     Z(K,I)=C*Z(K,I)-S*F
                  ENDDO
               ENDDO
               D(L)=D(L)-P
               E(L)=G
               E(M)=0.0D0
               GO TO 1
            ENDIF
         ENDDO
      ENDIF
      RETURN
      END
CCCCCC
      SUBROUTINE ODR(N,A)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N)
      DO  40  L=1,N
      ATEST=1.0E+12
      DO 41 J=L,N
      IF(A(J)-ATEST)42,41,41
   42 ATEST=A(J)
      JTEST=J
   41 CONTINUE
      A(JTEST)=A(L)
      A(L)=ATEST
   40 CONTINUE
      RETURN
      END
