C----------------------------------------------------------------------
C  September 2, 2010
C----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,A1,XI,AA,FA
      COMPLEX*16 SS,SS1,PSIIT,PSINIT 
      COMPLEX*16 S3
      PARAMETER (N_R=128,N_T=64,N_P=128)
      PARAMETER (NRTP=N_R*N_T*N_P)
      PARAMETER (JTOT=4)
      PARAMETER (JROT=1)
      COMMON/AMAT/AA(2*JTOT+1,2*JROT+1)
      DIMENSION RHO(N_R),THE(N_T),PHI(N_P)
      DIMENSION CLGD(-JROT:JROT),CLGD1(2*JROT+1)
      DIMENSION PSIIT(NRTP,2*JTOT+1)
      DIMENSION PSINIT(NRTP,2*JTOT+1)
      DIMENSION PLNITA(30,30)
      LTOT=JTOT
      ZI=DCMPLX(0.0D0,1.0D0)
      PI=4.0D0*DATAN(1.0D0)
      KMORSE=0
      HBAR=0.06350781278D0
      DE=4.5810D0 
      BE=1.942D0 
      RE=0.7417D0
      AM1=1.0075D0
      AM2=1.0075D0
      AM3=2.015D0
      EKIN=1.0D0
      AM=AM1+AM2+AM3
      AMR=AM1*AM2/(AM1+AM2)
      TMR=DSQRT(AM1*AM2*AM3/AM)
      AK_0=DSQRT(2.0D0*TMR*EKIN)/HBAR
      S=DSQRT(2.0D0*AMR*DE)*2.0D0/(HBAR*BE)
      R2_0=3.5D0
      R2_F=1.68D0
      SIGMA=0.21D0
      QQ=2.0D0*(R2_0-R2_F)/AK_0
      A1=1.0D0/(4.0D0*SIGMA*SIGMA-ZI*QQ)
      A2=8.0D0*SIGMA*SIGMA/(16.0D0*SIGMA**4+QQ*QQ)
      AN=(A2/PI)**0.25D0
CC
      CALL AKMU
CC
C  CLEBSCH – GORDAN COEFFICIENTS
CC
      WRITE(6,*)
      WRITE(6,*)' ','JROT',' ','M',' ','LTOT',' ','0'
     1,' ','JTOT',' ','M',' ','VAL',' ','CG COFF.'
      WRITE(6,*)
      DO M=0,JROT,1
         CALL CLD(JROT,M,LTOT,0,JTOT,M,VAL)
         WRITE(6,*)JROT,M,LTOT,0,JTOT,M,VAL
         WRITE(6,*)
         IF (M.EQ.0.0D0) THEN
            CLGD(M)=VAL
         ELSE
            CLGD(M)=VAL
            CLGD(-M)=VAL*(-1)**(-JTOT+JROT+LTOT)
         ENDIF
      ENDDO
      WRITE(6,*)
      WRITE(6,*) 'CLEBSCH – GORDAN COEFFICIENTS'
      WRITE(6,*)
      DO I=-JROT,JROT,1
      WRITE(6,*)CLGD(I)
      ENDDO
CC
C NORMALIZATION OF THE CLEBSCH – GORDAN COEFFICIENTS
CC
      S1=0.0D0
      II=1.0D0
      DO M=-JROT,JROT,1
         CLGD1(II)=CLGD(M)
         S1=S1+(2.0*LTOT+1.0D0)*CLGD1(II)*CLGD1(II)
         II=II+1
      ENDDO
      WRITE(6,*) 
      WRITE(6,*)'NORM. OF THE CLEBSCH – GORDAN COEFFICIENTS'
      WRITE(6,*) 
      WRITE(6,*)'S1=',S1
      WRITE(6,*) 
      S2=0.0D0
      DO M=1,2*JROT+1
         DO L=1,2*JTOT+1
            S2=S2+CONJG(AA(L,M))*AA(L,M)*(2*LTOT+1)*CLGD1(M)*CLGD1(M)
         ENDDO
      ENDDO
      WRITE(6,*) 
      WRITE(6,*)'NORM. OF CLEBSCH – GORDAN WITH A_k,mu'
      WRITE(6,*) 
      WRITE(6,*)'S2=',S2
      WRITE(6,*) 
CC
      WRITE(6,*)
      WRITE(6,*)'S3'
      WRITE(6,*)
      DO L=1,2*JTOT+1
         S3=DCMPLX(0.0D0,0.0D0)
         DO M=1,2*JROT+1
            S3=S3+AA(L,M)*CLGD1(M)
         ENDDO
         WRITE(6,*)DBLE(S3),DIMAG(S3)
         WRITE(6,*)
      ENDDO
      WRITE(6,*) 
      WRITE(6,*)'NORM. OF CLEBSCH – GORDAN WITH A_k,mu'
      WRITE(6,*) 
      WRITE(6,*)'S2=',S2
      WRITE(6,*) 

      RMIN=0.2D0
      RMAX=9.0D0
      DR=(RMAX-RMIN)/DFLOAT(N_R) 
      PI=4.0D0*DATAN(1.0D0)
      DT=PI/DFLOAT(N_T)
      DP=2.0D0*PI/DFLOAT(N_P)
      DO I=1,N_R
         RHO(I)=RMIN+DFLOAT(I-1)*DR+0.5D0*DR
      ENDDO
      DO I=1,N_T
         THE(I)=DFLOAT(I-1)*DT+0.5D0*DT
      ENDDO
      DO I=1,N_P
         PHI(I)=DFLOAT(I-1)*DP+0.5D0*DP
      ENDDO
      DO I=1,N_R
         RH=RHO(I)
         DO J=1,N_T
            TH=THE(J)
            JJ=J+(I-1)*N_T
            DO K=1,N_P
               PH=PHI(K)
               KK=K+(JJ-1)*N_P
               R1=DSQRT(0.5D0*RH*RH*(1.0D0+DSIN(TH)*DCOS(PH)))
               R2=DSQRT(0.5D0*RH*RH*(1.0D0-DSIN(TH)*DCOS(PH)))
               DN=1.0D0-DSIN(TH)*DSIN(TH)*DCOS(PH)*DCOS(PH)
               DNO=DSQRT(DN) 
               SNT=DCOS(TH)/DNO
               CNT=-DSIN(TH)*DSIN(PH)/DNO
               AJOC=RH*DSIN(TH)/(2.0D0*DNO)
               ARG=DEXP(-BE*(R1-RE))
               PSIR=DSQRT(BE)*AFI(KMORSE,S,ARG)
CC   
               XI=AN*CDEXP(-ZI*AK_0*R2-A1*(R2-R2_0)**2)
CC
               FA=DSQRT(2.0D0*PI*DABS(SNT)*AJOC*(2.0D0*LTOT+1)*0.5D0)
     1            *PSIR*XI*(-1)**(JROT-LTOT)
      
               CALL BLEG(CNT,PLNITA,JROT+1,1) 
CC
               DO L=1,2*JTOT+1
                  SS=DCMPLX(0.0D0,0.0D0)
                  DO M=1,2*JROT+1
			NP=ABS(M-JROT-1)
                     SS=SS+CLGD1(M)*AA(L,M)*PLNITA(JROT+1,NP+1) 
                  ENDDO
                  PSIIT(KK,L)=FA*SS
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      IF (MOD(JROT,2).EQ.0) THEN
         DO L=1,2*JTOT+1
            DO KK=1,NRTP
               PSINIT(KK,L)=PSIIT(KK,L)
            ENDDO
         ENDDO
      ENDIF 
      IF (MOD(JROT,2).EQ.1) THEN
         DO L=1,JTOT+1
            IF (L.LT.JTOT+1) THEN
               L1=L
               L2=2*JTOT+2-L
               DO KK=1,NRTP
                  PSINIT(KK,L1)=0.5D0*(PSIIT(KK,L1)+PSIIT(KK,L2))
                  PSINIT(KK,L2)=PSINIT(KK,L1)
               ENDDO
            ELSE
               DO KK=1,NRTP
                  PSINIT(KK,L)=PSIIT(KK,L)
               ENDDO
            ENDIF
         ENDDO 
      ENDIF
      SS2=0.0D0
      WRITE(6,*)
      WRITE(6,*)'WEIGHT FACTORS BEFORE NORMALIZATION'
      WRITE(6,*)
      DO L=1,2*JTOT+1
         SS1=DCMPLX(0.0D0,0.0D0)
         DO KK=1,NRTP
            SS1=SS1+CONJG(PSINIT(KK,L))*PSINIT(KK,L)*DR*DT*DP
         ENDDO
         WRITE(6,44)L-JTOT-1,DBLE(SS1),DIMAG(SS1)
         SS2=SS2+DBLE(SS1)
      ENDDO
      WRITE(6,*)
      WRITE(6,*)'NORM.'
      WRITE(6,*)
      WRITE(6,*)SS2
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)
      DO KK=1,NRTP
         DO L=1,2*JTOT+1
            PSINIT(KK,L)=PSINIT(KK,L)/SQRT(SS2)
         ENDDO
      ENDDO
      WRITE(6,*)
      WRITE(6,*)'WEIGHT FACTORS AFTER NORMALIZATION'
      WRITE(6,*)
      SS2=0.0D0
      DO L=1,2*JTOT+1
         SS1=DCMPLX(0.0D0,0.0D0)
         DO KK=1,NRTP
            SS1=SS1+CONJG(PSINIT(KK,L))*PSINIT(KK,L)*DR*DT*DP
         ENDDO
         WRITE(6,44)L-JTOT-1,DBLE(SS1),DIMAG(SS1)
         SS2=SS2+DBLE(SS1)
      ENDDO
      WRITE(6,*)
      WRITE(6,*)'NORM.'
      WRITE(6,*)
      WRITE(6,*)SS2
  44  FORMAT(1X,I4,1X,4F12.6)
  55  FORMAT(1X,4F12.6)
      STOP
      END

C----------------------------------------------------------------------
C MORSE WAVEFUNCTION
C----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION AFI(N,S,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION RL(50)
      RN=DFLOAT(N)
      IF (Y.EQ.0.0D0) GO TO 5
      FN=0.0D0
      IF (N.LT.2) GO TO 1
C     DO 2 I=2,N
C   2 FN=FN/I
C LOG N! IS CALCULATED
      IF (N.GT.20) GO TO 13
      DO I=2,N
         FN=FN-DLOG(DFLOAT(I))
      ENDDO
      FN=0.5D0*FN
      GO TO 1
   13 XX=DFLOAT(N+1)
      IFAIL=0
      DLNF=GAMMLN(XX)
      FN=-0.5D0*DLNF
    1 XX=S-RN
      IFAIL=0
      DLNG=GAMMLN(XX)
      XX=XX-RN
      IFAIL=0
      DLNG1=GAMMLN(XX)
      XX=XX-1.0D0
      IFAIL=0
      DLNG2=GAMMLN(XX)
      ARG=0.5D0*(DLNG-DLNG1-DLNG2)
      ARG=ARG-0.5D0*S*Y+(0.5D0*S-1.0D0*RN-0.5D0)*DLOG(S*Y)
      ARG=ARG+FN
      IF (N.GT.-1) GO TO 10
      FA=DSQRT(FN)*DEXP(ARG)
      N1=N+1
      SIG=1.0D0
      SU=1.0D0
      IF (N.EQ.0) GO TO 4
      DO I=2,N1
         I1=I-1
         SIG=-SIG
         F=1.0D0
         DO K=1,I1
            F=F*(N1-K)/K/(1.-(2.0D0*RN+1.0D0-DFLOAT(K))/S)
         ENDDO
         SU=SU+F*SIG*Y**(1.0D0*I1)
      ENDDO
    4 AFI=FA*SU
      GO TO 7
    5 AFI=0.0D0
    7 RETURN
   10 ALF=S-2.0D0*RN-1.0D0
      RL(1)=1.0D0
      Z=S*Y
      RL(2)=ALF+1.0D0-Z
      NS=N+1
      DO I=3,NS
         RN1=FLOAT(I-1)
         RL(I)=((2.*RN1-1.+ALF-Z)*RL(I-1)-(RN1-1.+ALF)*RL(I-2))/RN1
      ENDDO
      XX=1.0D0*RN+ALF+1.0D0
      IFAIL=0
      DLNG=GAMMLN(XX)
      ALF1=ALF+1.0D0
      IFAIL=0
      DLNG1=GAMMLN(ALF1)
      RN1=DFLOAT(N+1)
      IFAIL=0
      DLNG2=GAMMLN(RN1)
      FARG=DLNG2+DLNG1-DLNG
      IF (RL(NS).LT.2.0D-9) GO TO 12
      DLNG=DLOG(RL(NS))
      AFI=DEXP(ARG+FARG+DLNG)
      RETURN
   12 AFI=DEXP(ARG+FARG)*RL(NS)
      RETURN
   14 AN=1.0D0
      Z=S*Y
      ALF=S-2.0D0*RN-1.0D0
      DO M=N,1,-1
         RM=DFLOAT(M)
         AN=1.D0-Z*AN*DFLOAT(N-M+1)/RM/(ALF+RM)
      ENDDO
      IF (ARG.GT.10.0D0) WRITE(6,100) Y,ARG,AN
  100 FORMAT(1X,5E10.3)
      AFI=AN*DEXP(ARG)
      RETURN
      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION GAMMLN(XX)
      DOUBLE PRECISION XX,ZZ,TERM,RZ2
C      WRITE(6,*) 'ARG=',XX
      IER=0
      ZZ=XX
      IF(XX-1.0D10) 2,2,1
    1 IF(XX-1.0D30) 8,9,9
C
C SEE IF XX IS NEAR ZERO OR NEGATIVE
C
    2 IF(XX-1.0D-9) 3,3,4
    3 IER=-1
      GAMMLN=-1.0D30
      GO TO 10
C
C XX GREATER THAN ZERO AND LESS THAN OR EQUAL TO 1.D10
C
    4 TERM=1.0D0
    5 IF(ZZ-18.D0) 6,6,7
    6 TERM=TERM*ZZ
      ZZ=ZZ+1.0D0
      GO TO 5
    7 RZ2=1.0D0/ZZ**2
      GAMMLN=(ZZ-0.5D0)*DLOG(ZZ)-ZZ+0.9189385332046727D0-DLOG(TERM)+
     1(1.0D0/ZZ)*(0.83333333333333333D-1 -(RZ2*(0.27777777777777777D-2
     2+(RZ2*(0.793650793650793650D-3-(RZ2*(0.5952380952380952D0)))))))
      GO TO 10
C
C XX GREATER THAN 1.D+10 AND LESS THAN 1.D+70
C
    8 GAMMLN=ZZ*(DLOG(ZZ)-1.0D0)
      GO TO 10
C
C XX GREATER THAN OR EQUAL TO 1.D+70
    9 IER=+1
      WRITE(6,*)'ARGUMENT LNGAM > 1.D30'
      GAMMLN=1.0D+30
   10 RETURN
      END

C----------------------------------------------------------------------  
C----------------------------------------------------------------------  
      SUBROUTINE CLD(J_1,M_1,J_2,M_2,J,M,CLGD)
      IMPLICIT REAL*8(A-H,O-Z)
      IF (M.EQ.M_1+M_2) THEN
         DELTA=1.0D0
      ELSE
         DELTA=0.0D0
      ENDIF
      AJ_1=DFLOAT(J_1)
      AM_1=DFLOAT(M_1)
      AJ_2=DFLOAT(J_2)
      AM_2=DFLOAT(M_2)
      AJ=DFLOAT(J)
      AM=DFLOAT(M)
      JJ1=DINT(AJ+AJ_1-AJ_2)
      JJ2=DINT(AJ-AJ_1+AJ_2)
      JJ3=DINT(-AJ+AJ_1+AJ_2)
      JJ4=DINT(AJ_1+AJ_2+AJ+1)
      FAC1=DSQRT((2.0D0*AJ+1.0D0)*FACT(JJ1)*FACT(JJ2)
     1           *FACT(JJ3)/FACT(JJ4))
      JM1=DINT(AJ+AM)
      JM2=DINT(AJ-AM)
      JM3=DINT(AJ_1-AM_1)
      JM4=DINT(AJ_1+AM_1)
      JM5=DINT(AJ_2-AM_2)
      JM6=DINT(AJ_2+AM_2)
      FAC2=DSQRT(FACT(JM1)*FACT(JM2)*FACT(JM3)
     1           *FACT(JM4)*FACT(JM5)*FACT(JM6))
      FAC3=0.0D0
      DO K=0,20
         JMK1=INT(AJ_1+AJ_2-AJ-DFLOAT(K))
         JMK2=INT(AJ_1-AM_1-DFLOAT(K))
         JMK3=INT(AJ_2+AM_2-DFLOAT(K))
         JMK4=INT(AJ-AJ_2+AM_1+DFLOAT(K))
         JMK5=INT(AJ-AJ_1-AM_2+DFLOAT(K))
         IF(JMK1.LT.0) GO TO 2
         IF(JMK2.LT.0) GO TO 2
         IF(JMK3.LT.0) GO TO 2
         IF(JMK4.LT.0) GO TO 2
         IF(JMK5.LT.0) GO TO 2
         VAL=FACT(K)*FACT(JMK1)*FACT(JMK2)*FACT(JMK3)
     1       *FACT(JMK4)*FACT(JMK5)
         FAC3=FAC3+(-1.0D0)**K/VAL
    2    CONTINUE
      ENDDO
      CLGD=DELTA*FAC1*FAC2*FAC3
      CLGD=CLGD*(-1.0D0)**(AJ_2-AJ_1-AM)/DSQRT(2.0D0*AJ+1.0D0)
      RETURN 
      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------  
      FUNCTION FACT(N) 
      IMPLICIT REAL*8(A-H,O-Z)
      IF (N.EQ.0) THEN
         VAL=1.0D0
      ELSE
         VAL=1.0D0
         DO I=1,N
            VAL=VAL*DFLOAT(I)
         ENDDO
      ENDIF
      FACT=VAL
      RETURN
      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------  
      SUBROUTINE AKMU
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,AA,FA,BB
      PARAMETER (N_T=64,N_P=128)
      PARAMETER (JTOT=4)
      PARAMETER (JROT=1)
      COMMON/AMAT/AA(2*JTOT+1,2*JROT+1)
      DIMENSION PLT(30,30)
      DIMENSION PLP(30,30)
      DIMENSION BB(2*JROT+1,2*JROT+1)
      PI=4.0D0*DATAN(1.0D0)
      DT=PI/DFLOAT(N_T)
      DP=2.0D0*PI/DFLOAT(N_P)
      ZI=DCMPLX(0.0D0,1.0D0)
      K1=1
      DO KK=-JTOT,JTOT,1
         M1=1
         DO MM=-JROT,JROT,1
            AA(K1,M1)=DCMPLX(0.0D0,0.0D0)
            MPJ=ABS(KK)
            NPJ=ABS(MM)
            DO II=1,N_T
               T=(II-1)*DT+0.5D0*DT
               X=DCOS(T)
               CALL BLEG(X,PLT,JTOT+1,1) 
               DO JJ=1,N_P
                  P=(JJ-1)*DP+0.5D0*DP
                  Y=-DSIN(T)*DCOS(P)
                  CALL BLEG(Y,PLP,JTOT+1,1) 
                  DENO=SQRT(1.0D0-DSIN(T)*DSIN(T)*DCOS(P)*DCOS(P))
                  CSXI=DSIN(T)*DSIN(P)/DENO
                  SIXI=DCOS(T)/DENO
                  CSXI2=2.0D0*CSXI*CSXI-1.0D0
                  SIXI2=2.0D0*SIXI*CSXI
CC
C                  XI=DASIN(DCOS(T)/SQRT(1.0D0-Y*Y))
C                  FA=(-1)**KK*EXP(-ZI*KK*P+ZI*MM*XI)
C    1                *PLT(JTOT+1,MPJ+1)*PLP(JTOT+1,NPJ+1)
CC
                  IF (MM.EQ.-2) THEN
                     FA=(-1)**KK*EXP(-ZI*KK*P)*(CSXI2-ZI*SIXI2)
     1                  *PLT(JTOT+1,MPJ+1)*PLP(JTOT+1,NPJ+1)
                  ENDIF
                  IF (MM.EQ.-1) THEN
                     FA=(-1)**KK*EXP(-ZI*KK*P)*(CSXI-ZI*SIXI)
     1                  *PLT(JTOT+1,MPJ+1)*PLP(JTOT+1,NPJ+1)
                  ENDIF
                  IF (MM.EQ.0) THEN
                     FA=(-1)**KK*EXP(-ZI*KK*P)
     1                  *PLT(JTOT+1,MPJ+1)*PLP(JTOT+1,NPJ+1)
                  ENDIF
                  IF (MM.EQ.1) THEN
                     FA=(-1)**KK*EXP(-ZI*KK*P)*(CSXI+ZI*SIXI)
     1                  *PLT(JTOT+1,MPJ+1)*PLP(JTOT+1,NPJ+1)
                  ENDIF
                  IF (MM.EQ.2) THEN
                     FA=(-1)**KK*EXP(-ZI*KK*P)*(CSXI2+ZI*SIXI2)
     1                  *PLT(JTOT+1,MPJ+1)*PLP(JTOT+1,NPJ+1)
                  ENDIF
                  AA(K1,M1)=AA(K1,M1)+FA*DSIN(T)*DT*DP
               ENDDO
            ENDDO
            M1=M1+1 
         ENDDO
         K1=K1+1 
      ENDDO
      DO MM=1,2*JROT+1
         DO NN=1,2*JROT+1
            BB(MM,NN)=DCMPLX(0.0D0,0.0D0)
            DO KK=1,2*JTOT+1
               BB(MM,NN)=BB(MM,NN)+CONJG(AA(KK,MM))*AA(KK,NN)
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*)'RE(A_K,mu)'
      WRITE(6,*)
      WRITE(6,55)((DBLE(AA(II,JJ)),JJ=1,2*JROT+1),II=1,2*JTOT+1)
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'Im(A_K,mu)'
      WRITE(6,*)
      WRITE(6,55)((DIMAG(AA(II,JJ)),JJ=1,2*JROT+1),II=1,2*JTOT+1)
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'RE(A_K,mu*A_k,mu`)'
      WRITE(6,*)
      WRITE(6,55)((DBLE(BB(II,JJ)),JJ=1,2*JROT+1),II=1,2*JROT+1)
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*)'Im(A_K,mu*A_k,mu`)'
      WRITE(6,*)
      WRITE(6,66)((DIMAG(BB(II,JJ)),JJ=1,2*JROT+1),II=1,2*JROT+1)
  44  FORMAT(1X,f12.6) 
  55  FORMAT(1X,3f12.6) 
  66  FORMAT(1X,5f12.6) 
      RETURN
      END

C----------------------------------------------------------------------
C THIS PROGRAM CALCULATES ASSOCIATED LEGENDRE POLYNOMIALS MULTIPLIED 
C BY ((2L+1)(L-M)!/(L+M)!/2)**.5/SQRT(2*PHI), I.E. SPHERICAL HARMONICS
C ARE
C             Y(L,M,TETA,PHI)=P(L+1,M+1)*EXP(I*M*PHI)
C   IOP=0  :  P(L)*FAC
C   IOP=1  :  P(L,M)*FAC
C   IOP=2  :  P(L)
C   IOP=3  :  P(L,M) 
C  
C   WHERE  FAC=((2*L+1)(L-M)!/(2(L+M)!)**0.5 / SQRT(2PI) 
C  
C NOTE THAT X=COS(TE) IN BLEG IS ARGUMENT WHEREAS IN ALEG TE IS
C ARGUMENT
C----------------------------------------------------------------------
      SUBROUTINE BLEG(X,P,N,IOP) 
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(30,30)
      DO I=1,N
         DO J=1,N
            P(I,J)=0.0D0
         ENDDO 
      ENDDO 
      P(1,1)=1.0D0
      SX=DSQRT(1.-X**2)
      P(2,1)=X 
      N1=N+1
      DO I=3,N1
C         P(I,1)=((2*I-3)*P(I-1,1)*X-(I-2)*P(I-2,1))/(I-1)
         P(I,1)=2*X*P(I-1,1)-P(I-2,1)-(X*P(I-1,1)-P(I-2,1))/(I-1) 
      ENDDO
      IF (IOP.EQ.2) RETURN 
      IF (IOP.EQ.0) GO TO 4
      DO M=1,N
         M1=M+1
         DO I=1,N
            I1=I+1
            IF (I1.LT.M1) GO TO 2
            I2=I1-2  
            P(I1,M1)=(2*I1-3)*SX*P(I1-1,M1-1)
            IF (M1.LE.I2) P(I1,M1)=P(I1,M1)+P(I2,M1)
    2       CONTINUE
         ENDDO
      ENDDO 
      IF(IOP.EQ.3) RETURN  
      DO I=1,N1
         DO J=1,I
            L=I-1
            M=J-1
            P(I,J)=0.3989422804D0*P(I,J)*DSQRT((2.0D0*L+1.0D0)*
     1             FAC(L-M)/2.0D0/FAC(L+M))
         ENDDO
      ENDDO 
      RETURN
    4 DO I=1,N1
         P(I,1)=0.2820947918D0*DSQRT(2.0D0*DFLOAT(I)-1.0D0)*P(I,1)
      ENDDO
      RETURN
      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------  
      FUNCTION FAC(N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARR(10)
      DATA ARR/1.0D0,2.0D0,6.0D0,24.0D0,120.0D0,720.0D0,5040.0D0,
     140320.0D0,362880.0D0,3628800.0D0/
      IF (N.LT.1) GO TO 1  
      IF (N.GT.10) GO TO 4 
      FAC=ARR(N)
      GO TO 5  
    4 SU=1.0D0 
      DO I=1,N
         SU=SU*FLOAT(I)
      ENDDO
      GO TO 3  
    1 SU=1.0D0 
    3 FAC=SU
    5 RETURN
      END
