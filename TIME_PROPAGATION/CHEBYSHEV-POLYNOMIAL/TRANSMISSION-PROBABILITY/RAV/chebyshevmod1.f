      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 PSI
      PARAMETER (N_R=512)
      PARAMETER (NSTEP=2500)
      PARAMETER (NSTP=1024)
      PARAMETER (NE=140)
      PARAMETER (NCHMIN=10)
      PARAMETER (NCHMAX=100)
      COMMON/VOLLL/VOLEM,ISTEP
      DIMENSION BSL(0:NCHMAX)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      DIMENSION R(N_R)
      DIMENSION PS(N_R)
      DIMENSION V0(N_R)
      DIMENSION AKR(N_R)
      DIMENSION VIM(N_R)
      DIMENSION AT(NSTP),ET(NE)
      DIMENSION AKI(NE),CKI(NE)
      DIMENSION NN(1)
      DIMENSION URE(NSTEP),UIM(NSTEP)
      DIMENSION URT(2*NSTP)
      DIMENSION URRE(NSTP)
      DIMENSION URIM(NSTP)
      DIMENSION URR(NE)
CCCCCC
      EPSL1=1.0D-08
      EPSL2=1.0D-07
      NCDEL=2
CCCCCC
      EVEPS=0.9648533822D0
      EPSKJ=100.0D0
      PI=4.0D0*DATAN(1.0D0)
      HBAR=0.06350781278D0
      RMIN=-20.0D0
      RMAX=20.0D0
      AM=1.0D0
      DR=(RMAX-RMIN)/DFLOAT(N_R)
      VOLEM=DR
      AK0=20.0D0
      EKIN0=HBAR*HBAR*AK0*AK0/2.0D0/AM
      WRITE(6,*)EKIN0
      VMIN=0.0D0
      VMAX=1.00D0
      EKIN=PI*PI*HBAR*HBAR/2.0D0/AM/VOLEM/VOLEM
      EMIN=VMIN
      EMAX=EKIN+VMAX
CCCCCC
      MD=800
      AFLUX=0.0D0
      T=0.0D0
      DDT=0.01D0
      EGRIDP=0.50D0*(EMAX+EMIN)
      EGRIDM=0.50D0*(EMAX-EMIN)
      NCHEB=INT(DDT*EGRIDM/HBAR+0.5D0)
      NCHEB=NCHEB+10
      WRITE(6,*)'CHEBYSHEV'
      WRITE(6,'(1X,2F12.6,I8)')EKIN,DDT*EGRIDM/HBAR,NCHEB
CCCCCC
      EGRID1=DDT*EGRIDM/HBAR
      DO J=0,NCHMAX
      IF(J.EQ.0)THEN
        BSL(J)=BESSJ0(EGRID1)
      ELSE
        IF(J.EQ.1)THEN
          BSL(J)=BESSJ1(EGRID1)
        ELSE
          BSL(J)=BESSJ(J,EGRID1)
        ENDIF
      ENDIF
      ENDDO
CCCCCC
      CALL TRANSF(PSIRE,PSIIM)
      SS=0.0D0
      DO I=1,N_R
        R(I)=RMIN+(I-1)*DR+0.5D0*DR
        PS(I)=PSIRE(I)*PSIRE(I)+PSIIM(I)*PSIIM(I)
        SS=SS+PS(I)*VOLEM
        WRITE(30,'(1X,2F12.6)')R(I),PS(I)
      ENDDO
      STOP 121
      WRITE(6,*)'SS=',SS
      DO I=1,N_R
        RH=R(I)
        CALL POT(RH,EN)
        V0(I)=EN
        WRITE(7,'(2F12.6)')RH,EN
      ENDDO
CCCCCC
      CALL DEFAK(DR,AKR)
      CALL VPOTIM(R,MD,VIM)
!***********************************************************
!                                                          *
!     Time Propagation by Chebyshev polynomial expansion   *
!     scheme                                               *
!                                                          *
!***********************************************************
      MMM=NCHEB
      ISTEP=0
  50  ISTEP=ISTEP+1
      CALL ENERGY(PSIRE,PSIIM,AKR,V0,VAL1,VAL2,VAL)
        EER1=VAL1*VOLEM
        EER2=VAL2*VOLEM
        EER=VAL*VOLEM
      WRITE(22,'(1X,4F12.6)')T,EER1,EER2,EER
CCCCCC
      T=T+DDT
      CALL CHEBYSHEV(BSL,MMM,DDT,EGRIDP,EGRIDM,
     1AKR,V0,PSIRE,PSIIM)
CCCCCC
      SS=0.0D0
      DO I=1,N_R
        SS=SS+(PSIRE(I)**2+PSIIM(I)**2)*VOLEM
      ENDDO
      WRITE(21,88)T,SS,EER/SS
      IF(MOD(ISTEP,100).EQ.0)THEN
      DO I=1,N_R
        PS(I)=PSIRE(I)*PSIRE(I)+PSIIM(I)*PSIIM(I)
        WRITE(30+ISTEP/100,'(1X,3F12.6)')T,R(I),PS(I)
      ENDDO
      ENDIF
      WRITE(10,'(1X,I8,2F12.6)')ISTEP,
     1          PSIRE(N_R-MD-200),PSIIM(N_R-MD-200)
CCCCCC
      SUM=0.0D0
      DO I=N_R/2+1,N_R
        SUM=SUM+(PSIRE(I)*PSIRE(I)
     1     +PSIIM(I)*PSIIM(I))*VOLEM
      ENDDO
      WRITE(11,'(1X,2F12.6)')T,SUM
      CALL FLUX(PSIRE,PSIIM,VOLEM,AM,MD,VAL)
      AFLUX=AFLUX+VAL*DDT
      WRITE(12,'(1X,2F12.6)')T,AFLUX
CCCCCC
      DO I=1,N_R
        PSIRE(I)=PSIRE(I)*EXP(-VIM(I)*DDT/HBAR)
        PSIIM(I)=PSIIM(I)*EXP(-VIM(I)*DDT/HBAR)
      ENDDO
CCCCCC
      IF (ISTEP.EQ.NSTEP) GO TO 40
      GO TO 50
  40  CONTINUE
  88  FORMAT(1X,4F22.6)
      STOP
      END
!***********************************************************
!                                                          *
!     NORMALIZATION CHECK OF THE TRANSLATIONAL             *
!     WAVEFUNCTION (GWP)                                   *
!                                                          *
!***********************************************************
      SUBROUTINE TRANSF(PSIRE,PSIIM)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,XI
      PARAMETER (N_R=512)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      ZI=DCMPLX(0.0D0,1.0D0)
      PI=4.0D0*DATAN(1.0D0)
      HBAR=0.06350781278D0
      AM=1.0D0
      AK0=20.0D0
      RMIN=-20.0D0
      RMAX=20.0D0
      DR=(RMAX-RMIN)/DFLOAT(N_R)
      R0=-3.0D0
      DELTAR=0.2D0
      AN=1.0D0/(2.0D0*PI*DELTAR*DELTAR)
      AN=AN**0.25D0
      DO I=1,N_R
        R2=RMIN+(I-1)*DR+0.5D0*DR
        QQ=(R2-R0)*(R2-R0)/4.0D0/DELTAR/DELTAR
        XI=AN*EXP(ZI*AK0*R2-QQ)
        PSIRE(I)=DBLE(XI)
        PSIIM(I)=DIMAG(XI)
      ENDDO
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
      PARAMETER (N_R=512)
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
      ENDDO
      RETURN
      END
!***********************************************************
!                                                          *
!     HAMILTONIAN OPERATES ON THE WAVE FUNCTION            *
!     BY THE FFT-METHOD TO CALCULATE ENERGY                *
!                                                          *
!***********************************************************
      SUBROUTINE ENERGY(X0,Y0,AKR,V0,VAL1,VAL2,VAL)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (N_R=512)
      DIMENSION X0(N_R),Y0(N_R)
      DIMENSION AKR(N_R)
      DIMENSION V0(N_R)
      DIMENSION APSI(2*N_R)
      DIMENSION HPSIR(2*N_R)
      DIMENSION NN(1)
CCCCCC
      AM=1.0D0
      HBAR=0.06350781278D0
      RMST=HBAR*HBAR/AM
      SN=1.0D0/DSQRT(DFLOAT(N_R))
      DO I=1,N_R
         APSI(2*I-1)=X0(I)
         APSI(2*I)=Y0(I)
      ENDDO
      NN(1)=N_R
      CALL FOURN(APSI,NN,1,1)
      DO I=1,N_R
        HPSIR(2*I-1)=APSI(2*I-1)*AKR(I)*AKR(I)*SN
        HPSIR(2*I)=APSI(2*I)*AKR(I)*AKR(I)*SN
      ENDDO
      CALL FOURN(HPSIR,NN,1,-1)
      VAL1=0.0D0
      VAL2=0.0D0
      DO I=1,N_R
         HPSIRE1=0.5D0*RMST*HPSIR(2*I-1)*SN
         HPSIRE2=V0(I)*X0(I)
         HPSIIM1=0.5D0*RMST*HPSIR(2*I)*SN
         HPSIIM2=V0(I)*Y0(I)
         VAL1=VAL1+HPSIRE1*X0(I)+HPSIIM1*Y0(I)
         VAL2=VAL2+HPSIRE2*X0(I)+HPSIIM2*Y0(I)
      ENDDO
      VAL=VAL1+VAL2
      RETURN
      END
CCCCCC
      SUBROUTINE CHEBYSHEV(BSL,NN,T,EGRIDP,EGRIDM,
     1AKR,V0,PSIRE,PSIIM)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,FACT,PHICHE,SUM,PSI
      PARAMETER (N_R=512)
      PARAMETER (NCHMAX=100)
      DIMENSION BSL(0:NCHMAX)
      DIMENSION AKR(N_R)
      DIMENSION V0(N_R)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      DIMENSION PHICHE(N_R,0:NCHMAX)
      ZI=DCMPLX(0.0D0,1.0D0)
      HBAR=0.06350781278D0
      FACT=EXP(-ZI*T*EGRIDP/HBAR) 
      CALL RECURSION(T,BSL,AKR,V0,NN,EGRIDP,EGRIDM,PSIRE,PSIIM,PHICHE)
      DO I=1,N_R
        SUM=DCMPLX(0.0D0,0.0D0)
        DO J=0,NN
          IF(J.EQ.0)THEN
            EPSIL=1.0D0
          ELSE
            EPSIL=2.0D0
          ENDIF
        SUM=SUM+EPSIL*BSL(J)*PHICHE(I,J)
        ENDDO
        PSIRE(I)=DBLE(SUM*FACT)
        PSIIM(I)=DIMAG(SUM*FACT)
      ENDDO
      RETURN
      END
CCCCCC
      SUBROUTINE RECURSION(T,BSL,AKR,V0,NN,EGRIDP,EGRIDM,PSIRE,PSIIM,
     1PHICHE)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,ARR,PHICHE,FACT,SUM
      PARAMETER (N_R=512)
      PARAMETER (NCHMIN=10)
      PARAMETER (NCHMAX=100)
      COMMON/VOLLL/VOLEM,ISTEP
      DIMENSION AKR(N_R)
      DIMENSION V0(N_R)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      DIMENSION PHICHE(N_R,0:NCHMAX)
      DIMENSION BSL(0:NCHMAX)
CCCCCC
      EPSL1=1.0D-21
      EPSL2=1.0D-20
      NCDEL=2
CCCCCC
      NNI=2
      NNF=NN
      ZI=DCMPLX(0.0D0,1.0D0)
      HBAR=0.06350781278D0
      FACT=EXP(-ZI*T*EGRIDP/HBAR) 
      DO I=1,N_R
        PHICHE(I,0)=DCMPLX(PSIRE(I),PSIIM(I))
      ENDDO
      CALL XPRO(PSIRE,PSIIM,AKR,V0,EGRIDP,EGRIDM)
      DO I=1,N_R
        PHICHE(I,1)=-ZI*DCMPLX(PSIRE(I),PSIIM(I))
      ENDDO
  90  CONTINUE
      DO I=NNI,NNF
        DO J=1,N_R
          PSIRE(J)=DBLE(PHICHE(J,I-1))
          PSIIM(J)=DIMAG(PHICHE(J,I-1))
        ENDDO
        CALL XPRO(PSIRE,PSIIM,AKR,V0,EGRIDP,EGRIDM)
        DO J=1,N_R
          PHICHE(J,I)=-2.0D0*ZI*DCMPLX(PSIRE(J),PSIIM(J))
     1               +PHICHE(J,I-2)
        ENDDO
      ENDDO 
CCCCCC
      RAV=0.0D0
      DO I=1,N_R
        SUM=DCMPLX(0.0D0,0.0D0)
        DO J=NNF-4,NNF
          IF(J.EQ.0)THEN
            EPSIL=1.0D0
          ELSE
            EPSIL=2.0D0
          ENDIF
        SUM=SUM+EPSIL*BSL(J)*PHICHE(I,J)
        ENDDO
        RAV=RAV+CONJG(SUM*FACT)*SUM*FACT
      ENDDO
      RAV=RAV*VOLEM
      WRITE(14,*)ISTEP,RAV,NNI,NNF
CCCCCC
      IF(RAV.GE.EPSL1.AND.RAV.LE.EPSL2)GOTO 100
        IF(RAV.LT.EPSL1)THEN
          NN=NNF-NCDEL
        ENDIF
        IF(RAV.GT.EPSL2)THEN
        NNI=NNF+1
        NNF=NNF+NCDEL
        NN=NNF
        ENDIF
        IF(RAV.GT.EPSL2)GOTO 90
CCCCCC
  100 CONTINUE
      RETURN
      END
!***********************************************************
!                                                          *
!     HAMILTONIAN OPERATES ON THE WAVE FUNCTION            *
!     BY THE FFT-METHOD                                    *
!                                                          *
!***********************************************************
      SUBROUTINE XPRO(X0,Y0,AKR,V0,EGRIDP,EGRIDM)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (N_R=512)
      DIMENSION X0(N_R),Y0(N_R)
      DIMENSION AKR(N_R)
      DIMENSION V0(N_R)
      DIMENSION APSI(2*N_R)
      DIMENSION HPSIR(2*N_R)
      DIMENSION HPSIRE(N_R)
      DIMENSION HPSIIM(N_R)
      DIMENSION NN(1)
CCCCCC
      AM=1.0D0
      HBAR=0.06350781278D0
      RMST=HBAR*HBAR/AM
      SN=1.0D0/DSQRT(DFLOAT(N_R))
      DO I=1,N_R
        APSI(2*I-1)=X0(I)
        APSI(2*I)=Y0(I)
      ENDDO
      NN(1)=N_R
      CALL FOURN(APSI,NN,1,1)
      DO I=1,N_R
        HPSIR(2*I-1)=APSI(2*I-1)*AKR(I)*AKR(I)*SN
        HPSIR(2*I)=APSI(2*I)*AKR(I)*AKR(I)*SN
      ENDDO
      CALL FOURN(HPSIR,NN,1,-1)
      VAL1=0.0D0
      VAL2=0.0D0
      DO I=1,N_R
        HPSIRE(I)=0.5D0*RMST*HPSIR(2*I-1)*SN
     1           +V0(I)*X0(I)
        HPSIIM(I)=0.5D0*RMST*HPSIR(2*I)*SN
     1           +V0(I)*Y0(I)
        HPSIRE(I)=(HPSIRE(I)-EGRIDP*X0(I))/EGRIDM
        HPSIIM(I)=(HPSIIM(I)-EGRIDP*Y0(I))/EGRIDM
      ENDDO
      DO I=1,N_R
        X0(I)=HPSIRE(I)
        Y0(I)=HPSIIM(I)
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
CCCCCC
      SUBROUTINE DEFAK(DR,AKR)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N_R=512)
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
CCCCCC
      FUNCTION BESSJ(n,x)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (IACC=40,BIGNO=1.e10,BIGNI=1.e-10)
C     USES bessj0,bessj1
      if(n.lt.2)stop 'bad argument n in bessj'
      ax=abs(x)
      if(ax.eq.0.)then
        bessj=0.
      else if(ax.gt.float(n))then
        tox=2./ax
        bjm=bessj0(ax)
        bj=bessj1(ax)
        do 11 j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
11      continue
        bessj=bj
      else
        tox=2./ax
        m=2*((n+int(sqrt(float(IACC*n))))/2)
        bessj=0.
        jsum=0
        sum=0.
        bjp=0.
        bj=1.
        do 12 j=m,1,-1
          bjm=j*tox*bj-bjp
          bjp=bj
          bj=bjm
          if(abs(bj).gt.BIGNO)then
            bj=bj*BIGNI
            bjp=bjp*BIGNI
            bessj=bessj*BIGNI
            sum=sum*BIGNI
          endif
          if(jsum.ne.0)sum=sum+bj
          jsum=1-jsum
          if(j.eq.n)bessj=bjp
12      continue
        sum=2.*sum-bj
        bessj=bessj/sum
      endif
      if(x.lt.0..and.mod(n,2).eq.1)bessj=-bessj
      return
      END
CCCCCC
      FUNCTION BESSJ0(x)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     1s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     1s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     1-.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,
     1.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     1651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,
     1s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,
     159272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*
     1(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
        bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     1p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END
CCCCCC
      FUNCTION BESSJ1(x)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     1s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     1s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     1242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,
     1s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,
     199447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     1.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     1-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     1y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-2.356194491
        bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     1p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.,x)
      endif
      return
      END
CCCCCC
      SUBROUTINE FLUX(PSIRE,PSIIM,VOLEM,AM,MD,VAL)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (N_R=512)
      DIMENSION PSIRE(N_R),PSIIM(N_R)
      HBAR=0.06350781278D0
      MDP1=MD+1
      MDP2=MD+2
      MDM1=MD-1
      MDM2=MD-2
      PSIREDR=(-PSIRE(N_R-MDM2)+8.0D0*PSIRE(N_R-MDM1)
     1       -8.0D0*PSIRE(N_R-MDP1)+PSIRE(N_R-MDP2))/12.0D0/VOLEM
      PSIIMDR=(-PSIIM(N_R-MDM2)+8.0D0*PSIIM(N_R-MDM1)
     1       -8.0D0*PSIIM(N_R-MDP1)+PSIIM(N_R-MDP2))/12.0D0/VOLEM
      VAL=(HBAR/AM)*(PSIRE(N_R-MD)*PSIIMDR-PSIIM(N_R-MD)*PSIREDR)
      RETURN
      END
