      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 PSI1
      PARAMETER (N_R=512)
      PARAMETER (NSTEP=20000)
      PARAMETER (LANMIN=10)
      PARAMETER (LANMAX=30)
      DIMENSION AKR(N_R)
      DIMENSION RR(N_R)
      DIMENSION VPOT(N_R)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      DIMENSION NN(1)
      DIMENSION VIM(N_R)
!***********************************************************
!                                                          *
!     FOR FFTW ROUTINES                                    *
!                                                          *
!***********************************************************
      INTEGER*8 PLAN1,PLAN2
      DOUBLE COMPLEX ARR1,ARR
      DIMENSION ARR1(N_R)
      DIMENSION ARR(N_R)
      INTEGER FFTW_R2HC
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_TIMELIMIT
      PARAMETER (FFTW_TIMELIMIT=1073741824)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_NO_DFT_R2HC
      PARAMETER (FFTW_NO_DFT_R2HC=512)
      INTEGER FFTW_NO_NONTHREADED
      PARAMETER (FFTW_NO_NONTHREADED=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)
      INTEGER FFTW_NO_SLOW
      PARAMETER (FFTW_NO_SLOW=262144)
      INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
      PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
      INTEGER FFTW_ALLOW_PRUNING
      PARAMETER (FFTW_ALLOW_PRUNING=1048576)
CCCCCC
      EVEPS=0.9648533822D0
      PI=4.0D0*DATAN(1.0D0)
      HBAR=0.06350781278D0
      AM=1.008D0
      RMST=HBAR*HBAR/AM
      RINV=1.0D0/HBAR
      SNG=1.0D0/DSQRT(DFLOAT(N_R))
      LAITER=30
      NLAN=LAITER
      EPSL1=1.0D-08
      EPSL2=1.0D-07
      NLDEL=2
      ISTEP=0
      EKIN=0.20D0
      EKIN=EKIN*EVEPS
      AK0=DSQRT(2.0D0*AM*EKIN)/HBAR
      WRITE(6,'(2F16.8)')EKIN,AK0
CCCCCC
      RMIN=-4.0D0
      RMAX=8.0D0
      DR=(RMAX-RMIN)/DFLOAT(N_R)
      VOLEM=DR
      DO I=1,N_R
        RR(I)=RMIN+(I-1)*DR
      ENDDO
      CALL INITIALIZATION(RR,DR,AK0,AM,PSIRE,PSIIM)
      DO I=1,N_R
        RR1=RR(I)
        CALL POTENTIAL(RR1,EN)
        VPOT(I)=EN
        WRITE(8,'(2F12.6)')RR1,EN
      ENDDO
CCCCCC
      AFLUX=0.0D0
      TOUT=0.0D0
      DDT=0.010D0
CCCCCC
      CALL VPOTIM(RR,VIM)
      CALL DEFAK(DR,AKR)
      CALL DFFTW_PLAN_DFT_1D(PLAN1,N_R,ARR1,ARR1,
     &                       FFTW_FORWARD,FFTW_MEASURE)
      CALL DFFTW_PLAN_DFT_1D(PLAN2,N_R,ARR,ARR,
     &                       FFTW_BACKWARD,FFTW_MEASURE)
!***********************************************************
!                                                          *
!         Time Propagation by Lanczos iteration            *
!                                                          *
!***********************************************************
  50  ISTEP=ISTEP+1
      TOUT=TOUT+DDT
      MMM=NLAN
      CALL LANCZ(PSIRE,PSIIM,DDT,VOLEM,RAV,MMM,ISTEP,AKR,
     1RMST,RINV,SNG,VPOT,PLAN1,PLAN2)
      RAV=RAV*VOLEM
      IF (EPSL1.GT.0.0D0) THEN
        IF(RAV.LT.EPSL1) NLAN=NLAN-NLDEL
        IF(RAV.GT.EPSL2) NLAN=NLAN+NLDEL
        IF(NLAN.GT.LANMAX) NLAN=LANMAX
        IF(NLAN.LT.LANMIN) NLAN=LANMIN
      ENDIF
      WRITE(10,*)ISTEP,RAV,NLAN
      CALL ENERGY(PSIRE,PSIIM,AKR,RMST,SNG,VPOT,
     1VAL1,VAL2,VAL,VOLEM,PLAN1,PLAN2)
      EER1=VAL1*VOLEM
      EER2=VAL2*VOLEM
      EER=VAL*VOLEM
      WRITE(11,'(1X,4F12.6)')TOUT,EER1,EER2,EER
      SS=0.0D0
      DO I=1,N_R
        SS=SS+(PSIRE(I)**2+PSIIM(I)**2)*VOLEM
      ENDDO
      WRITE(12,'(1X,4F12.6)')TOUT,SS,EER/SS
      IF(MOD(ISTEP,500).EQ.0)THEN
        DO I=1,N_R
          PSI1=DCMPLX(PSIRE(I),PSIIM(I))
          PSI2=PSIRE(I)*PSIRE(I)+PSIIM(I)*PSIIM(I)
C          WRITE(30+ISTEP/500,'(1X,5F12.6)')TOUT,RR(I),ABS(PSI1),PSI2
        ENDDO
      ENDIF
      MD=48
      CALL FLUX(PSIRE,PSIIM,VOLEM,AM,MD,VAL)
      AFLUX=AFLUX+VAL*DDT
      WRITE(13,'(1X,2F12.6)')TOUT,AFLUX
CCCCCC
      DO I=1,N_R
        PSIRE(I)=PSIRE(I)*EXP(-VIM(I)*DDT/HBAR)
        PSIIM(I)=PSIIM(I)*EXP(-VIM(I)*DDT/HBAR)
      ENDDO
CCCCCC
      IF (ISTEP.GE.NSTEP) GO TO 40
      GO TO 50
  40  CONTINUE
      CALL DFFTW_DESTROY_PLAN(PLAN1)
      CALL DFFTW_DESTROY_PLAN(PLAN2)
      STOP 121
      END
!***********************************************************
!                                                          *
!     NORMALIZATION CHECK OF THE TRANSLATIONAL             *
!     WAVEFUNCTION (GWP)                                   *
!                                                          *
!***********************************************************
      SUBROUTINE INITIALIZATION(RR,DR,AK0,AM,PSIRE,PSIIM)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,XI,SS
      PARAMETER (N_R=512)
      DIMENSION RR(N_R)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      ZI=DCMPLX(0.0D0,1.0D0)
      PI=4.0D0*DATAN(1.0D0)
      HBAR=0.06350781278D0
      R0=5.0D0
      DELTAR=1.0D0
      AN=1.0D0/(2.0D0*PI*DELTAR*DELTAR)
      AN=AN**0.25D0
      SS=DCMPLX(0.0D0,0.0D0)
      DO I=1,N_R
        RR1=RR(I)
        QQ=(RR1-R0)*(RR1-R0)/4.0D0/DELTAR/DELTAR
        XI=AN*EXP(-ZI*AK0*RR1-QQ)
        PSIRE(I)=DBLE(XI)
        PSIIM(I)=DIMAG(XI)
        SS=SS+CONJG(XI)*XI*DR
        PSI2=PSIRE(I)*PSIRE(I)+PSIIM(I)*PSIIM(I)
        WRITE(7,'(1X,5F12.6)')RR1,DBLE(XI),DIMAG(XI),ABS(XI),PSI2
      ENDDO
      WRITE(6,*)'NORMALIZATION'
      WRITE(6,*)
      WRITE(6,'(1X,2F12.6)')DBLE(SS),DIMAG(SS)
      WRITE(6,*)
      RETURN
      END
!***********************************************************
!                                                          *
!     MODEL POTENTIAL                                      *
!                                                          *
!***********************************************************
      SUBROUTINE POTENTIAL(RR1,EN)
      IMPLICIT REAL*8 (A-H,O-Z)
      EVEPS=0.9648533822D0
      SIGMA=1.0D0
      V0=0.2D0*EVEPS
      EN=V0*EXP(-0.5D0*RR1*RR1/SIGMA/SIGMA)
      RETURN
      END
!***********************************************************
!                                                          *
!     ABSORBING POTENTIAL                                  *
!                                                          *
!***********************************************************
      SUBROUTINE VPOTIM(RR,VIM)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (N_R=512)
      DIMENSION RR(N_R),VIM(N_R)
      EVEPS=0.9648533822D0
      VMAX=0.2D0*EVEPS
      DO I=1,N_R
      VIM(I)=0.0D0
      ENDDO
      RRA=RR(1)
      RRB=RR(44)
      DO I=1,N_R
        IF(RR(I).LE.RRB)THEN
          VIM(I)=VMAX*(RR(I)-RRB)/(RRA-RRB)
        ENDIF
        WRITE(9,*)RR(I),VIM(I)
      ENDDO
      RETURN
      END
!***********************************************************
!                                                          *
!     CALCULATION OF FREQUENCIES                           *
!                                                          *
!***********************************************************
      SUBROUTINE DEFAK(DR,AKR)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N_R=512)
      DIMENSION AKR(N_R)
      PI=4.0D0*DATAN(1.0D0)
      PI2=PI*2.0D0
      ACONSR=PI2/(N_R*DR)
      NHALFR=N_R/2+1
      DO I=1,NHALFR
        AKR(I)=ACONSR*(I-1)
      ENDDO
      DO I=NHALFR+1,N_R
        AKR(I)=ACONSR*(I-1-N_R)
      ENDDO
      RETURN
      END
CCCCCC
      SUBROUTINE LANCZ(PSIRE,PSIIM,DDT,VOLEM,RAV,MMM,ISTEP,AKR,
     1RMST,RINV,SNG,VPOT,PLAN1,PLAN2)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*8 PLAN1,PLAN2
      PARAMETER (N_R=512)
      PARAMETER (LAD=30)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      DIMENSION U0A(N_R),U0B(N_R)
      DIMENSION U1A(N_R),U1B(N_R)
      DIMENSION WA(N_R,LAD),WB(N_R,LAD)
      DIMENSION TT(LAD,2),D(LAD),E(LAD),Z(LAD,LAD),RA(LAD),RB(LAD)
      DIMENSION AKR(N_R)
      DIMENSION VPOT(N_R)
      HBAR=0.06350781278D0
      SCALEA=DSQRT(VOLEM)
      SCALEI=1.0D0/SCALEA
      DO I=1,N_R
        U0A(I)=PSIRE(I)*SCALEA
        U0B(I)=PSIIM(I)*SCALEA
        U1A(I)=0.0D0
        U1B(I)=0.0D0
      ENDDO
      CALL CHAIN(U0A,U0B,U1A,U1B,WA,WB,TT,MMM,AKR,
     1RMST,RINV,SNG,VPOT,PLAN1,PLAN2)
C     IF (MMM.EQ.1) RETURN
      IF (MMM.EQ.1) GO TO 1
      DO I=1,N_R
        PSIRE(I)=0.0D0
        PSIIM(I)=0.0D0
      ENDDO
      DO I=2,MMM
         D(I)=TT(I,1)
         E(I)=TT(I-1,2)
      ENDDO
      D(1)=TT(1,1)
C     WRITE(6,1000) (D(I),I=1,MMM)
C     WRITE(6,1000) (E(I),I=1,MMM)
      DO I=1,MMM
        DO J=1,MMM
          Z(I,J)=0.0D0
          IF (I.EQ.J) Z(I,J)=1.0D0
        ENDDO
      ENDDO
C     DO 1002 I=1,MMM
C     PRINT 1000,I,D(I),E(I)
C1002  CONTINUE
1000  FORMAT(2X,2I5,2X,F15.5)
      CALL TQLI(D,E,MMM,Z)
      DO I=1,MMM
        SUA=0.0D0
        SUB=0.0D0
        DO J=1,MMM
          SUA=SUA+Z(I,J)*Z(1,J)*DCOS(D(J)*DDT)
          SUB=SUB-Z(I,J)*Z(1,J)*DSIN(D(J)*DDT)
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
      DO I=MMM-4,MMM
        RAV=RAV+RA(I)*RA(I)+RB(I)*RB(I)
      ENDDO
      RAV=RAV*0.2D0
C     WRITE(6,1000) TA,(BR(I),I=1,N)
      DO J=1,MMM
        DO I=1,N_R
          PSIRE(I)=PSIRE(I)+WA(I,J)*RA(J)-WB(I,J)*RB(J)
          PSIIM(I)=PSIIM(I)+WA(I,J)*RB(J)+WB(I,J)*RA(J)
        ENDDO
      ENDDO
C
    1 CONTINUE
      RETURN
      END
CCCCCC
      SUBROUTINE CHAIN(U0A,U0B,U1A,U1B,WA,WB,T,MMM,AKR,
     1RMST,RINV,SNG,VPOT,PLAN1,PLAN2)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (N_R=512)
      PARAMETER (LAD=30)
      DIMENSION U0A(N_R),U0B(N_R)
      DIMENSION U1A(N_R),U1B(N_R)
      DIMENSION X0(N_R),Y0(N_R)
      DIMENSION WA(N_R,LAD),WB(N_R,LAD)
      DIMENSION T(LAD,2)
      DIMENSION AKR(N_R)
      DIMENSION VPOT(N_R)
      INTEGER*8 PLAN1,PLAN2
      BN=0.0D0
      IT=0
   10 IT=IT+1
      DO I=1,N_R
        WA(I,IT)=U0A(I)
        WB(I,IT)=U0B(I)
        X0(I)=U0A(I)
        Y0(I)=U0B(I)
      ENDDO
      CALL XPRO(X0,Y0,AKR,RMST,RINV,SNG,VPOT,PLAN1,PLAN2)
C     IF (MMM.EQ.1) RETURN
      IF (MMM.EQ.1) GO TO 1
      AN=0.0D0
      DO 3 I=1,N_R
    3 AN=AN+U0A(I)*X0(I)+U0B(I)*Y0(I)
      T(IT,1)=AN
      DO I=1,N_R
        X0(I)=X0(I)-AN*U0A(I)-BN*U1A(I)
        Y0(I)=Y0(I)-AN*U0B(I)-BN*U1B(I)
      ENDDO
      BN=0.0D0
      DO I=1,N_R
        BN=BN+X0(I)*X0(I)+Y0(I)*Y0(I)
      ENDDO
      BN=DSQRT(BN)
      BN1=1.0D0/BN
      DO I=1,N_R
        U1A(I)=U0A(I)
        U1B(I)=U0B(I)
        U0A(I)=X0(I)*BN1
        U0B(I)=Y0(I)*BN1
      ENDDO
      T(IT,2)=BN
      IF (DABS(BN).LT.1.D-6) WRITE(6,289) IT,BN
  289 FORMAT(' BN= ',I4,E10.3)
      IF (IT.LT.MMM) GO TO 10
    1 CONTINUE
      RETURN
      END
!***********************************************************
!                                                          *
!     HAMILTONIAN OPERATES ON THE WAVE FUNCTION            *
!     BY THE FFT-METHOD                                    *
!                                                          *
!***********************************************************
      SUBROUTINE XPRO(X0,Y0,AKR,RMST,RINV,SNG,VPOT,
     1PLAN1,PLAN2)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE COMPLEX ARR1,ARR,ART,ARP,ARP1
      INTEGER*8 PLAN1,PLAN2
      PARAMETER (N_R=512)
C
      DIMENSION AKR(N_R)
      DIMENSION VPOT(N_R)
      DIMENSION X0(N_R),Y0(N_R)
      DIMENSION ARR1(N_R)
      DIMENSION ARR(N_R)
!***********************************************************
!                                                          *
!   This program is distributed with permission            *
!                                                          *
!***********************************************************
      INTEGER FFTW_R2HC
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_TIMELIMIT
      PARAMETER (FFTW_TIMELIMIT=1073741824)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_NO_DFT_R2HC
      PARAMETER (FFTW_NO_DFT_R2HC=512)
      INTEGER FFTW_NO_NONTHREADED
      PARAMETER (FFTW_NO_NONTHREADED=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)
      INTEGER FFTW_NO_SLOW
      PARAMETER (FFTW_NO_SLOW=262144)
      INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
      PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
      INTEGER FFTW_ALLOW_PRUNING
      PARAMETER (FFTW_ALLOW_PRUNING=1048576)
CCCCCC
      DO I=1,N_R
         ARR1(I)=DCMPLX(X0(I),Y0(I))
      ENDDO
      CALL DFFTW_EXECUTE_DFT(PLAN1,ARR1,ARR1)
      DO I=1,N_R
         ARR(I)=ARR1(I)*AKR(I)*AKR(I)*SNG
      ENDDO
      CALL DFFTW_EXECUTE_DFT(PLAN2,ARR,ARR)
      DO I=1,N_R
         HPSIRE1=0.5D0*RMST*DBLE(ARR(I))*SNG
         HPSIRE2=VPOT(I)*X0(I)
         HPSIIM1=0.5D0*RMST*DIMAG(ARR(I))*SNG
         HPSIIM2=VPOT(I)*Y0(I)
         X0(I)=HPSIRE1+HPSIRE2
         X0(I)=X0(I)*RINV
         Y0(I)=HPSIIM1+HPSIIM2
         Y0(I)=Y0(I)*RINV
      ENDDO
      RETURN
      END
CCCCCC
      SUBROUTINE TQLI(D,E,N,Z)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (LAD=30)
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
!***********************************************************
!                                                          *
!     HAMILTONIAN OPERATES ON THE WAVE FUNCTION            *
!     BY THE FFT-METHOD TO CALCULATE ENERGY                *
!                                                          *
!***********************************************************
      SUBROUTINE ENERGY(PSIRE,PSIIM,AKR,RMST,SNG,VPOT,
     1VAL1,VAL2,VAL,VOLEM,PLAN1,PLAN2)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ARR1,ARR
      PARAMETER (N_R=512)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      DIMENSION AKR(N_R)
      DIMENSION VPOT(N_R)
      DIMENSION ARR1(N_R),ARR(N_R)
      INTEGER*8 PLAN1,PLAN2
!**********************************************************
!                                                         *
!     This program is distributed with permission         *
!                                                         *
!**********************************************************
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_TIMELIMIT
      PARAMETER (FFTW_TIMELIMIT=1073741824)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_NO_DFT_R2HC
      PARAMETER (FFTW_NO_DFT_R2HC=512)
      INTEGER FFTW_NO_NONTHREADED
      PARAMETER (FFTW_NO_NONTHREADED=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)
      INTEGER FFTW_NO_SLOW
      PARAMETER (FFTW_NO_SLOW=262144)
      INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
      PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
      INTEGER FFTW_ALLOW_PRUNING
      PARAMETER (FFTW_ALLOW_PRUNING=1048576)
CCCCCC
      DO I=1,N_R
        ARR1(I)=DCMPLX(PSIRE(I),PSIIM(I))
      ENDDO
      CALL DFFTW_EXECUTE_DFT(PLAN1,ARR1,ARR1)
      DO I=1,N_R
        ARR(I)=ARR1(I)*AKR(I)*AKR(I)*SNG
      ENDDO
      CALL DFFTW_EXECUTE_DFT(PLAN2,ARR,ARR)
      VAL1=0.0D0
      VAL2=0.0D0
      DO I=1,N_R
        HPSIRE1=0.5D0*RMST*DBLE(ARR(I))*SNG
        HPSIRE2=VPOT(I)*PSIRE(I)
        HPSIIM1=0.5D0*RMST*DIMAG(ARR(I))*SNG
        HPSIIM2=VPOT(I)*PSIIM(I)
        VAL1=VAL1+HPSIRE1*PSIRE(I)+HPSIIM1*PSIIM(I)
        VAL2=VAL2+HPSIRE2*PSIRE(I)+HPSIIM2*PSIIM(I)
      ENDDO
      VAL=VAL1+VAL2
      RETURN
      END
CCCCCC
      SUBROUTINE FLUX(PSIRE,PSIIM,VOLEM,AM,MD,VAL)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (N_R=512)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      HBAR=0.06350781278D0
      MDP1=MD+1
      MDP2=MD+2
      MDP3=MD+3
      MDP4=MD+4
      MDM1=MD-1
      MDM2=MD-2
      MDM3=MD-3
      MDM4=MD-4
      PSIREDR=(3.0D0*PSIRE(MDP4)-32.0D0*PSIRE(MDP3)
     1       +168.0D0*PSIRE(MDP2)-672.0D0*PSIRE(MDP1)
     2       +672.0D0*PSIRE(MDM1)-168.0D0*PSIRE(MDM2)
     3       +32.0D0*PSIRE(MDM3)-3.0D0*PSIRE(MDM4))
     4       /840.0D0/VOLEM
      PSIIMDR=(3.0D0*PSIIM(MDP4)-32.0D0*PSIIM(MDP3)
     1       +168.0D0*PSIIM(MDP2)-672.0D0*PSIIM(MDP1)
     2       +672.0D0*PSIIM(MDM1)-168.0D0*PSIIM(MDM2)
     3       +32.0D0*PSIIM(MDM3)-3.0D0*PSIIM(MDM4))
     4       /840.0D0/VOLEM
      VAL=(HBAR/AM)*(PSIRE(MD)*PSIIMDR-PSIIM(MD)*PSIREDR)
      RETURN
      END
