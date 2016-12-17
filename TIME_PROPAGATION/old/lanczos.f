      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 PS
      PARAMETER (N_R=4096)
      PARAMETER (NSTEP=4500)
      PARAMETER (NSTP=4096)
      PARAMETER (NE=140)
      PARAMETER (LANMIN=10)
      PARAMETER (LANMAX=100)
      PARAMETER (NVIM=2000)
      DIMENSION AKR(N_R)
      DIMENSION R(N_R)
      DIMENSION V0(N_R)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      DIMENSION PSRE(N_R),PSIM(N_R)
      DIMENSION X1(N_R),Y1(N_R)
      DIMENSION PS(N_R),PS1(N_R)
      DIMENSION AT(NSTP),ET(NE)
      DIMENSION AKI(NE),CKI(NE)
      DIMENSION NN(1)
      DIMENSION VIM(N_R)
      DIMENSION URE(NSTEP),UIM(NSTEP)
      DIMENSION URT(2*NSTP)
      DIMENSION URRE(NSTP)
      DIMENSION URIM(NSTP)
      DIMENSION URR(NE)
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
      INA=0
      EVEPS=0.9648533822D0
      EPSKJ=100.0D0
      PI=4.0D0*DATAN(1.0D0)
      HBAR=0.06350781278D0
      RMIN=-20.0D0
      RMAX=20.0D0
      AM=1.0D0
      DR=(RMAX-RMIN)/DFLOAT(N_R)
      VOLEM=DR
      RMST=HBAR*HBAR/AM
      RINV=1.0D0/HBAR
      SN_G=1.0D0/DSQRT(DFLOAT(N_R))
      LAITER=30
      NLAN=LAITER
      EPSL1=1.0D-08
      EPSL2=1.0D-07
      NLDEL=2
      AK0=20.0D0
      EKIN=HBAR*HBAR*AK0*AK0/2.0D0/AM
      WRITE(6,*)EKIN
CCCCCC
      ISTEP=0
      CALL TRANSF(PSIRE,PSIIM)
      DO I=1,N_R
        R(I)=RMIN+(I-1)*DR
        PS(I)=DCMPLX(PSIRE(I),PSIIM(I))
        PS1(I)=PSIRE(I)*PSIRE(I)+PSIIM(I)*PSIIM(I)
        WRITE(30,'(1X,I8,4F12.6)')ISTEP,R(I),ABS(PS(I)),PS1(I)
      ENDDO
      DO I=1,N_R
        RH=R(I)
        CALL POT(RH,EN)
        V0(I)=EN
        WRITE(7,'(2F12.6)')RH,EN
      ENDDO
CCCCCC
      MD=500
      IDX=20
      AFLUX=0.0D0
      XG=0.0D0
      T=0.0D0
      DDT=0.010D0
CCCCCC
      CALL DEFAK(DR,AKR)
      CALL INITIAL_WEIGHT(EKIN,AM,IDX,DDT,ET,AKI,CKI)
      SUM=0.0D0
      DO I=1,NE
        SUM=SUM+CKI(I)
      ENDDO
      WRITE(8,*)SUM
      IF(INA.EQ.1)GOTO 100
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
      T=T+DDT
      CALL VPOTIM(R,MD,VIM)
      DO I=1,N_R
        X1(I)=PSIRE(I)
        Y1(I)=PSIIM(I)
      ENDDO
      CALL ENERGY(X1,Y1,AKR,RMST,SN_G,V0,
     1VAL1,VAL2,VAL,VOLEM,PLAN1,PLAN2)
        EER1=VAL1*VOLEM
        EER2=VAL2*VOLEM
        EER=VAL*VOLEM
        MMM=NLAN
      CALL LANCZ(X1,Y1,DDT,VOLEM,RAV,MMM,ISTEP,AKR,
     1RMST,RINV,SN_G,V0,PLAN1,PLAN2)
      RAV=RAV*VOLEM
      IF (EPSL1.GT.0.0D0) THEN
        IF(RAV.LT.EPSL1) NLAN=NLAN-NLDEL
        IF(RAV.GT.EPSL2) NLAN=NLAN+NLDEL
        IF(NLAN.GT.LANMAX) NLAN=LANMAX
        IF(NLAN.LT.LANMIN) NLAN=LANMIN
      ENDIF
      WRITE(14,*)ISTEP,RAV,NLAN
      DO I=1,N_R
        PSRE(I)=X1(I)
        PSIM(I)=Y1(I)
      ENDDO
      WRITE(22,'(1X,4F12.6)')T,EER1,EER2,EER
      DO I=1,N_R
        PSIRE(I)=PSRE(I)
        PSIIM(I)=PSIM(I)
      ENDDO
      SS=0.0D0
      DO I=1,N_R
        SS=SS+(PSRE(I)**2+PSIM(I)**2)*VOLEM
      ENDDO
      WRITE(21,'(1X,4F12.6)')T,SS,EER/SS
      IF(MOD(ISTEP,100).EQ.0)THEN
      DO I=1,N_R
        PS(I)=DCMPLX(PSRE(I),PSIM(I))
        PS1(I)=PSRE(I)*PSRE(I)+PSIM(I)*PSIM(I)
        WRITE(30+ISTEP/100,'(1X,5F12.6)')T,R(I),ABS(PS(I)),PS1(I)
      ENDDO
      ENDIF
CCCCCC
      WRITE(10,'(1X,I8,2F12.6)')ISTEP,PSRE(N_R-MD),PSIM(N_R-MD)
      DO I=1,N_R
        PSIRE(I)=PSIRE(I)*EXP(-VIM(I)*DDT/HBAR)
        PSIIM(I)=PSIIM(I)*EXP(-VIM(I)*DDT/HBAR)
      ENDDO
CCCCCC
      SUM=0.0D0
      DO I=N_R/2+1,N_R
        SUM=SUM+(PSRE(I)*PSRE(I)
     1     +PSIM(I)*PSIM(I))*VOLEM
      ENDDO
      WRITE(11,'(1X,2F12.6)')T,SUM
      CALL FLUX(PSRE,PSIM,VOLEM,AM,MD,VAL)
      AFLUX=AFLUX+VAL*DDT
      WRITE(12,'(1X,2F12.6)')T,AFLUX
CCCCCC
      IF (ISTEP.EQ.NSTEP) GO TO 40
      GO TO 50
  40  CONTINUE
      CALL DFFTW_DESTROY_PLAN(PLAN1)
      CALL DFFTW_DESTROY_PLAN(PLAN2)
  100 CONTINUE
CCCCCC
      IF(INA.EQ.1)THEN
      DO I=1,NSTEP
        READ(10,'(1X,I8,2F12.6)')ISTEP,URE(I),UIM(I)
      ENDDO
CCCCCC
      SN_T=1.0D0/DSQRT(DFLOAT(NSTP))
C     NSTART=200
      NSTART=0
      ISS=0
      DO IS=NSTART+1,NSTP+NSTART
        ISS=ISS+1
        IF(ISS.LE.NSTEP) THEN
          URT(2*ISS-1)=URE(IS)
          URT(2*ISS)=UIM(IS)
        ELSE
          URT(2*ISS-1)=0.0D0
          URT(2*ISS)=0.0D0
        ENDIF
      ENDDO

      NN(1)=NSTP

      CALL FOURN(URT,NN,1,1)

      ISS=0
      DO IS=1,NSTP
        ISS=ISS+1
        URRE(ISS)=URT(2*IS-1)*SN_T
        URIM(ISS)=URT(2*IS)*SN_T
      ENDDO
CCCCCC
      IS=IDX-1
      IEE=1
      DO IE=1,NE
      IS=IS+1
         ARE=URRE(IS)
         AIM=URIM(IS)
         URR(IEE)=AKI(IE)*AKI(IE)*HBAR*HBAR/CKI(IE)*(ARE**2+AIM**2)
      IEE=IEE+1
      ENDDO
      DO IE=1,NE
        WRITE(13,*)ET(IE)*EPSKJ,AKI(IE),URR(IE)
      ENDDO
CCCCCC
      SUM=0.0D0
      DK=AKI(2)-AKI(1)
      DO I=1,NE
        SUM=SUM+CKI(I)*URR(I)*DK
      ENDDO
      WRITE(15,*)SUM
CCCCCC      
      AKMIN=15.0D0
      AKMAX=24.0D0
      A1=1.0D0
      A2=2.0D0
      DO I=1,10
      AAK=AKMIN+(I-1)*1.0D0
      FAC1=2.0D0*PI*AAK/A2
      FAC2=PI*DSQRT(8.0D0*AM*A1/A2/A2/HBAR/HBAR-1.0D0)
      PTRANS=(DCOSH(FAC1)-1.0D0)/(DCOSH(FAC1)+DCOSH(FAC2))
      WRITE(16,*)AAK,PTRANS
      ENDDO
      ENDIF
CCCCCC
      STOP
      END
!***********************************************************
!                                                          *
!     HAMILTONIAN OPERATES ON THE WAVE FUNCTION            *
!     BY THE FFT-METHOD                                    *
!                                                          *
!***********************************************************
      SUBROUTINE XPRO(X0,Y0,AKR,RMST,RINV,SN,V0,
     1PLAN1,PLAN2)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE COMPLEX ARR1,ARR,ART,ARP,ARP1
      INTEGER*8 PLAN1,PLAN2
      PARAMETER (N_R=4096)
C
      DIMENSION AKR(N_R)
      DIMENSION V0(N_R)
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
         ARR(I)=ARR1(I)*AKR(I)*AKR(I)*SN
      ENDDO
      CALL DFFTW_EXECUTE_DFT(PLAN2,ARR,ARR)
      DO I=1,N_R
         HPSIRE1=0.5D0*RMST*DBLE(ARR(I))*SN
         HPSIRE2=V0(I)*X0(I)
         HPSIIM1=0.5D0*RMST*DIMAG(ARR(I))*SN
         HPSIIM2=V0(I)*Y0(I)
         X0(I)=HPSIRE1+HPSIRE2
         X0(I)=X0(I)*RINV
         Y0(I)=HPSIIM1+HPSIIM2
         Y0(I)=Y0(I)*RINV
      ENDDO
      RETURN
      END
CCCCCC
      SUBROUTINE LANCZ(X1,Y1,DT,VOLEM,RAV,M,ISTEP,AKR,
     1RMST,RINV,SN_G,V0,PLAN1,PLAN2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N_R=4096)
      PARAMETER (LAD=100)
      DIMENSION X1(N_R),Y1(N_R)
      DIMENSION U0A(N_R),U0B(N_R)
      DIMENSION U1A(N_R),U1B(N_R)
      DIMENSION WA(N_R,LAD),WB(N_R,LAD)
      DIMENSION TT(LAD,2),D(LAD),E(LAD),Z(LAD,LAD),RA(LAD),RB(LAD)
      DIMENSION AKR(N_R)
      DIMENSION V0(N_R)
      INTEGER*8 PLAN1,PLAN2
      HBAR=0.06350781278D0
      SCALEA=DSQRT(VOLEM)
      SCALEI=1.0D0/SCALEA
      DO I=1,N_R
         U0A(I)=X1(I)*SCALEA
         U0B(I)=Y1(I)*SCALEA
         U1A(I)=0.0D0
         U1B(I)=0.0D0
      ENDDO
      CALL CHAIN(U0A,U0B,U1A,U1B,
     1WA,WB,TT,M,AKR,RMST,RINV,SN_G,V0,
     3PLAN1,PLAN2)
C     IF (M.EQ.1) RETURN
      IF (M.EQ.1) GO TO 1
      DO I=1,N_R
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
         DO I=1,N_R
            X1(I)=X1(I)+WA(I,J)*RA(J)-WB(I,J)*RB(J)
            Y1(I)=Y1(I)+WA(I,J)*RB(J)+WB(I,J)*RA(J)
         ENDDO
      ENDDO
C
    1 CONTINUE
      RETURN
      END
CCCCCC
      SUBROUTINE CHAIN(U0A,U0B,U1A,U1B,
     1WA,WB,T,M,AKR,RMST,RINV,SN_G,V0,
     2PLAN1,PLAN2)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (N_R=4096)
      PARAMETER (LAD=100)
      DIMENSION U0A(N_R),U0B(N_R)
      DIMENSION U1A(N_R),U1B(N_R)
      DIMENSION X0(N_R),Y0(N_R)
      DIMENSION WA(N_R,LAD),WB(N_R,LAD)
      DIMENSION T(LAD,2)
      DIMENSION AKR(N_R)
      DIMENSION V0(N_R)
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
      CALL XPRO(X0,Y0,
     1AKR,RMST,RINV,SN_G,V0,
     2PLAN1,PLAN2)
C     IF (M.EQ.1) RETURN
      IF (M.EQ.1) GO TO 1
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
C----------------------------------------------------------------------
C----------------------------------------------------------------------  
      SUBROUTINE DEFAK(DR,AKR)
*
*     CALCULATION OF FREQUENCIES
*
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N_R=4096)
      DIMENSION AKR(N_R)
*
      PI=2.0D0*ASIN(1.0D0)
      PI2=PI*2.0D0
      ACONSR=PI2/(N_R*DR)
      NHALFR = N_R/2+1
      DO I=1,NHALFR
        AKR(I)=ACONSR*(I-1)
      ENDDO
      DO I=NHALFR+1,N_R
        AKR(I)=ACONSR*(I-1-N_R)
      ENDDO
      RETURN
      END
!***********************************************************
!                                                          *
!     HAMILTONIAN OPERATES ON THE WAVE FUNCTION            *
!     BY THE FFT-METHOD TO CALCULATE ENERGY                *
!                                                          *
!***********************************************************
      SUBROUTINE ENERGY(X1,Y1,AKR,RMST,SN,V0,
     1VAL1,VAL2,VAL,VOLEM,PLAN1,PLAN2)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ARR1,ARR
      PARAMETER (N_R=4096)
      DIMENSION X1(N_R),Y1(N_R)
      DIMENSION AKR(N_R)
      DIMENSION V0(N_R)
      DIMENSION ARR1(N_R),ARR(N_R)
      INTEGER*8 PLAN1,PLAN2
!*******************************************************************************
!
!   This program is distributed with permission
!
!*******************************************************************************
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
         ARR1(I)=DCMPLX(X1(I),Y1(I))
      ENDDO
      CALL DFFTW_EXECUTE_DFT(PLAN1,ARR1,ARR1)
      DO I=1,N_R
         ARR(I)=ARR1(I)*AKR(I)*AKR(I)*SN
      ENDDO
      CALL DFFTW_EXECUTE_DFT(PLAN2,ARR,ARR)
      VAL1=0.0D0
      VAL2=0.0D0
      DO I=1,N_R
         HPSIRE1=0.5D0*RMST*DBLE(ARR(I))*SN
         HPSIRE2=V0(I)*X1(I)
         HPSIIM1=0.5D0*RMST*DIMAG(ARR(I))*SN
         HPSIIM2=V0(I)*Y1(I)
         VAL1=VAL1+HPSIRE1*X1(I)+HPSIIM1*Y1(I)
         VAL2=VAL2+HPSIRE2*X1(I)+HPSIIM2*Y1(I)
      ENDDO
      VAL=VAL1+VAL2
      RETURN
      END
!***********************************************************
!                                                          !
!     ECKART POTENTIAL                                     !
!                                                          !
!***********************************************************
      SUBROUTINE POT(R,EN)
      IMPLICIT REAL*8 (A-H,O-Z)
      A1=1.0D0
      A2=2.0D0
      EN=A1/DCOSH(A2*R)/DCOSH(A2*R)
      RETURN
      END
!***********************************************************
!                                                          *
!     ABSORBING POTENTIAL                                  *
!                                                          *
!***********************************************************
      SUBROUTINE VPOTIM(R,MD,VIM)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (N_R=4096)
      DIMENSION R(N_R)
      DIMENSION VIM(N_R)
      HBAR=0.06350781278D0
      AM=1.0D0
      AK0=20.0D0
      DELTAR=0.2D0
      DELTAK=1.0D0/2.0D0/DELTAR
      AKMIN=AK0-DELTAK
      AKMAX=AK0+DELTAK
      VMAX=HBAR*HBAR*DSQRT(AKMAX*AKMIN*AKMIN*AKMIN)/2.0D0/AM
      DO I=1,N_R
      VIM(I)=0.0D0
      ENDDO
      DO I=1,N_R
      IF(I.LT.(N_R-MD).AND.I.GT.MD) GO TO 10
      IF(I.GT.(N_R-MD))THEN
      FAC=(R(I)-R(N_R-MD))/(R(N_R)-R(N_R-MD))
      VIM(I)=VMAX*FAC
      ENDIF
      IF(I.LT.MD)THEN
      FAC=(R(I)-R(MD))/(R(1)-R(MD))
      VIM(I)=VMAX*FAC
      ENDIF
   10 CONTINUE
C     WRITE(102,*)R(I),VIM(I)
      ENDDO
      RETURN
      END
!***********************************************************
!                                                          *
!     NORMALIZATION CHECK OF THE TRANSLATIONAL             *
!     WAVEFUNCTION (GWP)                                   *
!                                                          *
!***********************************************************
      SUBROUTINE TRANSF(PSIRE,PSIIM)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,XI,SS
      PARAMETER (N_R=4096)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      ZI=DCMPLX(0.0D0,1.0D0)
      PI=4.0D0*DATAN(1.0D0)
      HBAR=0.06350781278D0
      AM=1.0D0
      AK0=20.0D0
      AK0=-AK0
      RMIN=-20.0D0
      RMAX=20.0D0
      DR=(RMAX-RMIN)/DFLOAT(N_R)
      R0=-4.0D0
      DELTAR=0.2D0
      AN=1.0D0/(2.0D0*PI*DELTAR*DELTAR)
      AN=AN**0.25D0
      SS=DCMPLX(0.0D0,0.0D0)
      DO I=1,N_R
      R2=RMIN+(I-1)*DR
      QQ=(R2-R0)*(R2-R0)/4.0D0/DELTAR/DELTAR
      XI=AN*EXP(-ZI*AK0*R2-QQ)
      PSIRE(I)=DBLE(XI)
      PSIIM(I)=DIMAG(XI)
      SS=SS+CONJG(XI)*XI*DR
      WRITE(100,'(1X,4F12.6)')R2,DBLE(XI),DIMAG(XI),ABS(XI)
      ENDDO
      WRITE(6,*)'NORMALIZATION'
      WRITE(6,*)
      WRITE(6,'(1X,2F12.6)')DBLE(SS),DIMAG(SS)
      WRITE(6,*)
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
C----------------------------------------------------------------------
      SUBROUTINE INITIAL_WEIGHT(EKIN,AM,IDX,DDT,ET,AKI,CKI)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NSTP=4096)
      PARAMETER (NE=140)
      DIMENSION AT(NSTP),ET(NE)
      DIMENSION AKI(NE),CKI(NE)
CCCC
      PI=4.0D0*DATAN(1.0D0)
      DELTAR=0.2D0
      HBAR=0.06350781278D0
      EK=EKIN
CC
      AK_0=SQRT(2.0D0*AM*EK)/HBAR
      WRITE(6,*)
      WRITE(6,*)'AK_0',AK_0
      WRITE(6,*)
C
      WRITE(6,*)
      WRITE(6,*)'WEIGHT FACTOR FOR INCOMING WAVE FUNCTION'
      WRITE(6,*)
      CALL DEFAT(DDT,AT)
      IE=0
      DO I=IDX,NE+IDX-1
      IE=IE+1
      ET(IE)=HBAR*AT(I)
      AKI(IE)=SQRT(2.0D0*AM*ET(IE))/HBAR
      CKI(IE)=SQRT(2.0D0/PI)*DELTAR
     1       *EXP(-2.0D0*DELTAR*DELTAR*(AKI(IE)-AK_0)**2)
      WRITE(6,44)IE,ET(IE),AKI(IE),CKI(IE)
      ENDDO
  44  FORMAT(1X,I4,4X,4F16.6)
      RETURN
      END
CCCCCC
      SUBROUTINE FLUX(PSIRE,PSIIM,VOLEM,AM,MD,VAL)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (N_R=4096)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      HBAR=0.06350781278D0
      MD=N_R-MD-5
      MDP1=MD+1
      MDP2=MD+2
      MDM1=MD-1
      MDM2=MD-2
      MDP1=MD+1
      MDP2=MD+2
      MDP3=MD+3
      MDP4=MD+4
      MDM1=MD-1
      MDM2=MD-2
      MDM3=MD-3
      MDM4=MD-4

      PSIREDR=(-PSIRE(MDM2)+8.0D0*PSIRE(MDM1)
     1       -8.0D0*PSIRE(MDP1)+PSIRE(MDP2))
     2       /12.0D0/VOLEM
      PSIIMDR=(-PSIIM(MDM2)+8.0D0*PSIIM(MDM1)
     1       -8.0D0*PSIIM(MDP1)+PSIIM(MDP2))
     2       /12.0D0/VOLEM
C     PSIREDR=(PSIRE(N_R-MDM1)-PSIRE(N_R-MDP1))/2.0D0/VOLEM
C     PSIIMDR=(PSIIM(N_R-MDM1)-PSIIM(N_R-MDP1))/2.0D0/VOLEM
C      PSIREDR=(-PSIRE(N_R-MDM2)+8.0D0*PSIRE(N_R-MDM1)
C     1       -8.0D0*PSIRE(N_R-MDP1)+PSIRE(N_R-MDP2))/12.0D0/VOLEM
C      PSIIMDR=(-PSIIM(N_R-MDM2)+8.0D0*PSIIM(N_R-MDM1)
C     1       -8.0D0*PSIIM(N_R-MDP1)+PSIIM(N_R-MDP2))/12.0D0/VOLEM
C      PSIREDR=(PSIRE(MDP2)-8.0D0*PSIRE(MDP1)
C     1       +8.0D0*PSIRE(MDM1)-PSIRE(MDM2))
C     2       /12.0D0/VOLEM
C      PSIIMDR=(PSIIM(MDM2)-8.0D0*PSIIM(MDM1)
C     1       +8.0D0*PSIIM(MDP1)-PSIIM(MDP2))
C     2       /12.0D0/VOLEM
      VAL=(HBAR/AM)*(PSIRE(MD)*PSIIMDR-PSIIM(MD)*PSIREDR)
      RETURN
      END
CCCCCC
      SUBROUTINE DEFAT(DT,AT)
*
*     CALCULATION OF FREQUENCIES
*
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NSTP=4096)
      DIMENSION AT(NSTP)
*
      PI=2.0D0*ASIN(1.0D0)
      PI2=PI*2.0D0
      ACONST=PI2/(NSTP*DT)
      NHALFT = NSTP/2+1
      DO I=1,NHALFT
        AT(I)=ACONST*(I-1)
      ENDDO
      DO I=NHALFT+1,NSTP
        AT(I)=ACONST*(I-1-NSTP)
      ENDDO
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
