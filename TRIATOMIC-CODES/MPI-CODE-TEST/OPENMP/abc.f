C----------------------------------------------------------------------
CCC
CC  OCTOBER 19, 2010
CC This code is intended for A+BC time-dependent quantum reactive 
CC scattering in hyperspherical coordinates. The current working 
CC version is meant for adibatic processes that occur on a single-
CC sheeted potential energy surface, with the H3 SLTH being given 
CC as an illustration to the study of the D+H2 reaction. The theory
CC comes mostly from the work of Billing, Markovic and Muckerman: 
CC JCP 91, 6830 (1989); CPL 172, 509 (1990);JCP 99, 2674 (1993); 
CC JCP 97, 8201 (1992); CPL 195, 53 (1992); JCP 100, 1085 (1994).
CC The program consists basically of three parts: 
CC
CC I   - Initialization of the wave function by expressing the 
CC       K-component wave packet as a product of Morse functions, 
CC       Arthurs-Dalgarno states, and the translational wave function.
CC  
CC II  - Propagation of the wave packet by employing FFT [D. Kosloff 
CC       and R. Kosloff J. Comput. Phys. 52, 35 (1983)] followed by a 
CC       Lanczos iteration technique (R. T. Pack and J.C. Light, 
CC       JCP 85, 5870 (1986).
CC    
CC III - Projection of propagated wave packet into asymptotic 
CC       channels. This is done by first moving from hyperspherical 
CC       to Jacobi coordinates so as to make the process as efficient
CC       as possible. In order to do this, we require to interpolate 
CC       the time dependent wave function from hyperspherical to 
CC       Jacobi coordiantes such as to get the equilibrium values to 
CC       project on the pure diatomic state.
CCC

CCC
CC Unless not specified, all parameters have their usual meaning. 
CC The others have the following:
CC
CC N_R, N_T and N_P specify the 3D grid size in (rho,theta,phi).
CC
CC JTOT and JROT define to total angular momentum of the triatom 
CC and rotational quantum number of the diatom.
CC
CC NV and NJ are the highest vibrational and rotational quantum 
CC numbers used to perform the analysis.
CC
CC NCH refers the available total number of channels.
CC
CC NSTEP is the total number of time steps used for the propagation.
CC
CC NSTP is the number of time steps used for the analysis.
CC
CC NE is the energy grid used for the analysis.
CC
CC LAD is the max
CC
CC LANMIN is the minimum number of Lanczos iterations for propagation.
CC
CC LANMAX is the maximum number of Lanczos iterations for propagation.
CCC
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 PLAN1,PLAN2,PLAN3,PLAN4,PLAN5,PLAN6,PLAN7,PLAN8
      DOUBLE COMPLEX ARR1,ARR,ART,ARP,ARP1,ARR6,ARR7,ARR8
      COMPLEX*16 PSINIT
      COMPLEX*16 AA,AAP
      COMPLEX*16 HAN
      INTEGER TIME
      PARAMETER (N_R=128,N_T=64,N_P=128)
      PARAMETER (NT2=2*N_T)
      PARAMETER (NRTP=N_R*N_T*N_P)
      PARAMETER (NRTP2=N_R*NT2*N_P)
      PARAMETER (JTOT=6)
      PARAMETER (JROT=0)
      PARAMETER (JTR=JTOT+1)    ! FOR EVEN JROT
C      PARAMETER (JTR=JTOT)      ! FOR ODD JROT
      PARAMETER (NJRTP=NRTP*JTR)
      PARAMETER (NSTATE=3)
      PARAMETER (NSJRTP=NSTATE*NJRTP)
      PARAMETER (NV=10)
      PARAMETER (NJ=12)
      PARAMETER (NCH=3)
      PARAMETER (NSTEP=200)
      PARAMETER (NSTP=16384)
      PARAMETER (NE=1000)
      PARAMETER (LANMIN=14)
      PARAMETER (LANMAX=30)
      PARAMETER (LAD=30)

CCC
CC The following physical quantities are utilized:
CC
CC  DD(3),REQ(3),BETA(3): Mose parameters for the three diatoms
CC
CC  EM(3),EMU,DM(3),EPS(3),RM12(3),RM123(3): atomic masses, 
CC      triatomic reduced mass, channel dependent constant, 
CC      channel dependent angles, diatomic reduced masses,
CC      channel dependent reduced masses [Eq 61 of JCP 97, 8206 (1992)].
CC
CC  EKIN,RMIN,RMAX: kinetic energy, min and max values of rho.
CC
CC  RMST,RINV,SNG: hbar**2/2*emu, 1/hbar,1/sqrt(nstep)
CC
CC  RHO(N_R),THE(N_T),PHI(N_P): grid elements
CC
CC  RHOD(N_R),SINN(N_T),COSN(N_T): terms in 1/rho2, sin2, cos2 in
CC      Eq. 38 of JCP 99, 2676 (1993)
CC
CC  SIND(N_T),COSD(N_T),SIN2D(N_T): terms in 1/sin2, 1/cos2 and 
CC      1/sin(2theta)2 in Eq 38 of previous reference
CC  
CC  PSINIT(NRTP,2*JTOT+1): initial wave function having JTOT or JTOT+1
CC      elements for odd or even JROT. 
CC
CC  PSI(NRTP):total time dependent wf used to calculate norm, draw
CC      features of wf, etc 
CC
CC  AKR(N_R),AKT(N_T),AKP(N_P): Fourier frequencies for rho, 
CC      theta and phi
CC
CC  V0(NRTP): potential energy surface, in this case LSTH.
CC
CC  VIM(N_R): optical potential.
CC
CC  CKI(NE),ET(NE): weight of incoming gaussian wave packet at a given
CC      point of the energy grid (ET).
CC
CC  AT(NSTP): frequency in the transforme energy2time grid.
CC
CC  AKF(NE,0:NV,0:NJ): weight of outgoing wave function for all 
CC      rovibrational states.
CC
CC  NNV(NE),NNJ(NE,0:NV)storing variables for vibrational and 
CC      rotational quantum numbers at a given ET energy grid element.
CC
CC  ISTEP: ith time step.
CC
CC  AK_J,EJK,EJKP,EJKM: elements in Eq. 38 of JCP 99, 2676 (1993).
CC
CC  X1(NRTP),Y1(NRTP): real and imaginary parts of phik-2
CC  X2(NRTP),Y2(NRTP): real and imaginary parts of phik
CC  X3(NRTP),Y3(NRTP): real and imaginary parts of phik+2
CC
CC  VAL1,VAL2,VAL3,VAL4,VAL5,VAL6,VAL7,VAL8,VAL9,VAl10: contributions
CC      of different components of the Hmiltonian to the total energy as 
CC      in Eq 38 of JCP 99, 2676 (1993).
CC
CC  RSTAR: point chosen for analysis in Jacobi coordinates.
CC
CC  PRE(0:NV,0:NJ),PIM(0:NV,0:NJ): real and imaginary parts of 
CC      scattering amplitude in Eq. 26 of JCP 100.1087 (1994).
CC
CC  PSIRE(NJRTP): real part of time dependent wf in Eq 38 (ref above) 
CC  PSIIM(NJRTP): imag part of time dependent wf in Eq 38 (ref above)
CC  PSIRE1(NJRTP): real part of time dependent wf in Eq 38 (ref above)
CC  PSIIM1(NJRTP): imag part of time dependent wf in Eq 38 (ref above)
CC
CC  RNUC(3),DVDR(3): internuclear distances and their derivatives
CC
CC  CKF(NCH,NE,0:NV,0:NJ): channel-dependent weight of outgoing wave 
CC      function for all rovibrational states.
CC
CC  NLAN: number of Lanczos iterations for various k-component 
CC      of time-dependent wave function during propagation.
CC
CC  NR(3),NT(3): The absolute of the wavefunction is written
CC  as function of theta and phi for three fixed rho and
CC  as function of rho and phi for three fixed theta at different time
CC
CC  NN(1): NN(1) is as NSTEP the number of time steps used 
CC      for the time ener FT.
CC
CC  PROJ(0:NV,0:NJ): transition probability for a given vibrational-
CC      rotational state.
CC
CC  QRE(0:NV,0:NJ,NCH):real part of amplitude for a given vibrational-
CC      rotational state during propagation.
CC  QIM(0:NV,0:NJ,NCH):imaginary part of amplitude for a given vibrational-
CC      rotational state during propagation.
CC
CC  URE(NSTEP,0:NV,0:NJ,NCH): real part of amplitude for a given vibrational-
CC      rotational state during analysis.
CC  UIM(NSTEP,0:NV,0:NJ,NCH): imaginary part of amplitude for a given vibrational-
CC      rotational state during analysis.
CC  URT(2*NSTP): real/imaginary parts of time-dependent amplitude for a 
CC      particular (v,j) state.
CC
CC  URRE(NSTP,0:NV,0:NJ,NCH): real part of transition amplitude for a 
CC      particular (v,j) state as a function of energy.
CC  URIM(NSTP,0:NV,0:NJ,NCH): imaginary part of transition amplitude for a 
CC      particular (v,j) state as a function of energy.
CC
CC  URR(NE,0:NV,0:NJ,NCH): transition probability for each channel and all 
CC      rovibrational states; see Eq. 36 of JPC 100, 1088 (1994).
CC
CC  SE(NE),RR(NE,0:NV,NCH): normalizing factor at each energy, and 
CC      transition probability for each channel and all 
CC      vibrational states; see Eq. 36 of JPC 100, 1088 (1994).
CCC
      DIMENSION DD(3),REQ(3),BETA(3),EM(3)
      DIMENSION DM(3),EPS(3),RM12(3),RM123(3)
      DIMENSION AKR(N_R),AKT(NT2),AKP(N_P)
      DIMENSION RHO(N_R),THE(N_T),PHI(N_P)
      DIMENSION RHOD(N_R),SINN(N_T),COSN(N_T)
      DIMENSION SIND(N_T),COSD(N_T),SIN2D(N_T)
      DIMENSION WIN1(N_T),WIN2(N_T)
      DIMENSION PSINIT(NRTP,2*JTOT+1)
C      DIMENSION PSI(NSTATE*NRTP)
      DIMENSION AT(NSTP),ET(NE),CKI(NE)
      DIMENSION AKF1(NE,0:NV,0:NJ)
      DIMENSION AKF2(NE,0:NV,0:NJ)
      DIMENSION AKF3(NE,0:NV,0:NJ)
      DIMENSION NNV1(NE),NNJ1(NE,0:NV)
      DIMENSION NNV2(NE),NNJ2(NE,0:NV)
      DIMENSION NNV3(NE),NNJ3(NE,0:NV)
      DIMENSION V11(NRTP),V12(NRTP),V13(NRTP)
      DIMENSION V21(NRTP),V22(NRTP),V23(NRTP)
      DIMENSION V31(NRTP),V32(NRTP),V33(NRTP)
      DIMENSION VIM(N_R)
      DIMENSION CKF1(NCH,NE,0:NV,0:NJ)
      DIMENSION CKF2(NCH,NE,0:NV,0:NJ)
      DIMENSION CKF3(NCH,NE,0:NV,0:NJ)
      DIMENSION HANRE1(NCH,NE,0:NV,0:NJ,0:NJ+JTOT)
      DIMENSION HANIM1(NCH,NE,0:NV,0:NJ,0:NJ+JTOT)
      DIMENSION HANRE2(NCH,NE,0:NV,0:NJ,0:NJ+JTOT)
      DIMENSION HANIM2(NCH,NE,0:NV,0:NJ,0:NJ+JTOT)
      DIMENSION HANRE3(NCH,NE,0:NV,0:NJ,0:NJ+JTOT)
      DIMENSION HANIM3(NCH,NE,0:NV,0:NJ,0:NJ+JTOT)
      DIMENSION NNV(NE),NNJ(NE,0:NV)
      DIMENSION AAP(0:NJ,JTR,2*NJ+1),VALP(0:NJ,0:NJ+JTOT,2*NJ+1)
      DIMENSION PSIREIN(NJRTP),PSIIMIN(NJRTP)
C      DIMENSION PSIREIN1(NJRTP),PSIIMIN1(NJRTP)
      DIMENSION PSIRE(NSJRTP),PSIIM(NSJRTP)
      DIMENSION AK(JTR),EK(JTR),EKP(JTR),EKM(JTR)
      DIMENSION RNUC(3),DVDR(3)
      DIMENSION EER1(NSTATE*JTR),EER2(NSTATE*JTR),EER3(NSTATE*JTR)
      DIMENSION EER4(NSTATE*JTR),EER5(NSTATE*JTR),EER6(NSTATE*JTR)
      DIMENSION EER7(NSTATE*JTR),EER8(NSTATE*JTR),EER9(NSTATE*JTR)
      DIMENSION EER10(NSTATE*JTR),EER(NSTATE*JTR)
      DIMENSION NRR(3),NTT(3),NPP(3)
      DIMENSION PROJ(0:NV,0:NJ)
      DIMENSION ZTEMP(NSTATE,NSTATE),ZZZ(NSTATE)
      DIMENSION ZP(NSTATE,NSTATE,NRTP)
      DIMENSION PRE1(0:NV,0:NJ),PIM1(0:NV,0:NJ)
      DIMENSION PRE2(0:NV,0:NJ),PIM2(0:NV,0:NJ)
      DIMENSION PRE3(0:NV,0:NJ),PIM3(0:NV,0:NJ)
      DIMENSION QRE1(0:NV,0:NJ,NCH),QIM1(0:NV,0:NJ,NCH)
      DIMENSION QRE2(0:NV,0:NJ,NCH),QIM2(0:NV,0:NJ,NCH)
      DIMENSION QRE3(0:NV,0:NJ,NCH),QIM3(0:NV,0:NJ,NCH)
      DIMENSION URE1(NSTEP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION UIM1(NSTEP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URE2(NSTEP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION UIM2(NSTEP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URE3(NSTEP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION UIM3(NSTEP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URE11(NSTP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION UIM11(NSTP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URE21(NSTP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION UIM21(NSTP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URE31(NSTP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION UIM31(NSTP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URT1(2*NSTP),URT2(2*NSTP),URT3(2*NSTP)
      DIMENSION URRE1(NSTP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URIM1(NSTP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URRE2(NSTP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URIM2(NSTP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URRE3(NSTP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URIM3(NSTP,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION ARE1(NSTP,0:NV,0:NJ,NCH),AIM1(NSTP,0:NV,0:NJ,NCH)
      DIMENSION ARE2(NSTP,0:NV,0:NJ,NCH),AIM2(NSTP,0:NV,0:NJ,NCH)
      DIMENSION ARE3(NSTP,0:NV,0:NJ,NCH),AIM3(NSTP,0:NV,0:NJ,NCH)
      DIMENSION URR1(NE,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URR2(NE,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION URR3(NE,0:NV,0:NJ,NCH,0:NJ+JTOT)
      DIMENSION SUR1(NE,0:NV,0:NJ,NCH)
      DIMENSION SUR2(NE,0:NV,0:NJ,NCH)
      DIMENSION SUR3(NE,0:NV,0:NJ,NCH)
      DIMENSION SE(NE)
      DIMENSION RR1(NE,0:NV,NCH),RR2(NE,0:NV,NCH),RR3(NE,0:NV,NCH)
      DIMENSION POT(3,3)
      DIMENSION SS(NSTATE),RHAV(NSTATE)
      DIMENSION VALRC11(NCH,0:NV,0:NV)
      DIMENSION VALRC22(NCH,0:NV,0:NV)
      DIMENSION VALINF12(NCH,0:NV,0:NV)
      DIMENSION VALINF21(NCH,0:NV,0:NV)
      DIMENSION UU(NSTATE,NSTATE),W(NSTATE),ZS(NSTATE,NSTATE)
      DIMENSION FV1(NSTATE),FV2(NSTATE)
!
      DIMENSION U0A(NSJRTP),U0B(NSJRTP)
      DIMENSION U1A(NSJRTP),U1B(NSJRTP)
      DIMENSION WA(NSJRTP,LAD),WB(NSJRTP,LAD)
      DIMENSION TT(LAD,2),D(LAD),E(LAD),Z(LAD,LAD),RA(LAD),RB(LAD)
      DIMENSION EE(LAD)
      DIMENSION PRE(NRTP,JTR,NSTATE),PIM(NRTP,JTR,NSTATE)
      DIMENSION X1(NRTP),Y1(NRTP)
      DIMENSION X2(NRTP),Y2(NRTP)
      DIMENSION X3(NRTP),Y3(NRTP)
      DIMENSION ARR1(N_R,NT2,N_P)
      DIMENSION ARR(N_R,NT2,N_P)
      DIMENSION ART(N_R,NT2,N_P)
      DIMENSION ARP(N_R,NT2,N_P)
      DIMENSION ARP1(N_R,NT2,N_P)
      DIMENSION ARR6(NSTP)
      DIMENSION ARR7(NSTP)
      DIMENSION ARR8(NSTP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FOR FFTW ROUTINES                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
      OPEN(2,FILE='data.dat',FORM='UNFORMATTED',
     $STATUS='UNKNOWN')
      OPEN (5,FILE='abc.dat')
CCC
CC ALL QUANTITIES AND VARIABLES IN THIS PROGRAM WILL BE DEFINED IN  
CC SI UNITS (UNLESS SPECIFIED), E.G., THE VALUE OF HBAR APPEARS IN 
CC AMU*ANG**2/(10^-14 SEC). However, for easy use, the user may specify
CC them in the most conventional units, with the transformation being
CC done inside the code itself with the following conversion factors:
CCC   
      EVEPS=0.9648533822D0
      AUANG=0.5291770644D0
      AUEPS=26.25499575D0
CCC
CC NATURAL CONSTANTS
CCC
      PI=4.0D0*DATAN(1.0D0)
      HBAR=0.06350781278D0
CCC
CC  ORBITAL ANGULAR MOMENTUM = Total angular momentum
CCC
      LTOT=JTOT
CCC
CC  Read input data
CCC
CC  LANCZOS PARAMETERS: STARTING NO OF ITERATION, ACCURACY AND INCREMENT 
CC  OF LANCZOS ITERATION 
CC  INITIAL KE,MORSE QUANTUM NUMBER AND CHANNEL
CC  IN CASE OF PROPAGATION, INA=0 AND ANALYSIS, INA=1
CCC
      READ(5,*)LAITER
      READ(5,*)EPSL1,EPSL2,NLDEL
      READ(5,*)EKIN,KMORSE,ICH
      READ(5,*)INA,INIT
      READ(5,*)RMIN,RMAX
      READ(5,*)NSTART,IDX
!
!     HYPERSPHERICAL CONSTANTS
! 
      CALL MOSC(REQ,BETA,DD,EM)
      CALL HYPSET(DD,REQ,BETA,EM,EMU,DM,EPS,RM12,RM123)
!
!     CALCULATE INITIAL VIB/ROT ENERGY
!
      XDE=DD(ICH)
      XBETA=BETA(ICH)
      XREQ=REQ(ICH)
      AMR=RM12(ICH)
      EVJ=EVIBRO(XDE,XBETA,XREQ,AMR,KMORSE,JROT)
      EVJ=EVJ/EVEPS
      WRITE(21,*)
      WRITE(21,223)KMORSE,JROT,EVJ,EKIN+EVJ
      WRITE(21,*)
      CALL FLUSH(21)
  223 FORMAT('V=',I2,5X,'J=',I2,5X,'EVJ=',F8.6,5X,'ETOT=',F8.6)
!
!    THE HYPERSPHERICAL GRID (RHO, TH AND PHI).
!    FEW PRE-FACTORS IN THE OPERATORS OF HYPERSPHERICAL HAMILTONIAN.
!   
      DR=(RMAX-RMIN)/DFLOAT(N_R)
      DT=PI/DFLOAT(N_T)/2.0D0
      DP=2.0D0*PI/DFLOAT(N_P)
      VOLEM=DR*DT*DP
      RMST = HBAR*HBAR/EMU
      RINV=1.0D0/HBAR
      SNG=1.0D0/DSQRT(DFLOAT(NRTP2))
      BTA=1.5D0
      ANT1=DFLOAT(N_T/2)-8.0D0
      ANT2=DFLOAT(N_T/2)+9.0D0
!
      DO I=1,N_R
        RHO(I)=RMIN+DFLOAT(I-1)*DR+0.5D0*DR
        RH=RHO(I)
        RHOD(I)=1.0D0/RH/RH
        WRITE(11,'(2F12.6)')RH,RHOD(I)
        CALL FLUSH(11)
      ENDDO
      DO I=1,N_T
        THE(I)=DFLOAT(I-1)*DT+0.5D0*DT
        TH=THE(I)
        SINN(I)=DSIN(TH)
        COSN(I)=DCOS(TH)
        SIND(I)=1.0D0/DSIN(TH)/DSIN(TH)
        COSD(I)=1.0D0/DCOS(TH)/DCOS(TH)
        SIN2D(I)=1.0D0/DSIN(2.0D0*TH)/DSIN(2.0D0*TH)
        ARG1=2.0D0*BTA*(DFLOAT(I)-ANT1)
        ARG2=2.0D0*BTA*(ANT2-DFLOAT(I))
        WIN1(I)=1.0D0/(1.0D0+EXP(-ARG1))
        WIN2(I)=1.0D0/(1.0D0+EXP(-ARG2))
        WRITE(11,'(6F12.6)')TH,SIND(I),COSD(I),WIN1(I),WIN2(I),SIN2D(I)
        CALL FLUSH(11)
c        IF (SIND(I).GT.2.5D0) THEN
c           SIND(I)=2.5D0
c        ENDIF
c        IF (COSD(I).GT.2.5D0) THEN
c           COSD(I)=2.5D0
c        ENDIF
        IF (SIN2D(I).GT.2.5D0) THEN
          SIN2D(I)=2.5D0
        ENDIF
      ENDDO
      DO I=1,N_P
         PHI(I)=DFLOAT(I-1)*DP+0.5D0*DP
      ENDDO
!
!     INITIALIZATION OF THE WAVEFUNCTION
!
      KM=KMORSE
      IC=ICH
      JT=JTOT
      JR=JROT
      CALL PSIHYP(KM,JR,IC,EKIN,RMIN,RMAX,
     1DD,REQ,BETA,DM,EPS,RM12,RM123,PSINIT)
      IF(JTOT.NE.0)THEN
        IF(MOD(JROT,2).EQ.0)THEN
          J=0
          DO K=1,2*JTOT+1,2
            J=J+1
            DO II=1,NRTP
              JJ=II+(J-1)*NRTP
              PSIREIN(JJ)=DBLE(PSINIT(II,K))
              PSIIMIN(JJ)=DIMAG(PSINIT(II,K))
            ENDDO
          ENDDO
          JTT=J
        ENDIF
        IF(MOD(JROT,2).EQ.1)THEN
          J=0
          DO K=2,2*JTOT+1,2
            J=J+1
            DO II=1,NRTP
              JJ=II+(J-1)*NRTP
              PSIREIN(JJ)=DBLE(PSINIT(II,K))
              PSIIMIN(JJ)=DIMAG(PSINIT(II,K))
            ENDDO
          ENDDO
          JTT=J
        ENDIF
      ELSE
        J=0
        DO K=1,2*JTOT+1,2
          J=J+1
          DO II=1,NRTP
            JJ=II+(J-1)*NRTP
            PSIREIN(JJ)=DBLE(PSINIT(II,K))
            PSIIMIN(JJ)=DIMAG(PSINIT(II,K))
          ENDDO
        ENDDO
        JTT=J
      ENDIF
C      CALL PSIHYP(KM+1,JR,IC,EKIN,RMIN,RMAX,
C     1DD,REQ,BETA,DM,EPS,RM12,RM123,PSINIT)
C      IF(JTOT.NE.0)THEN
C        IF(MOD(JROT,2).EQ.0)THEN
C          J=0
C          DO K=1,2*JTOT+1,2
C            J=J+1
C            DO II=1,NRTP
C              JJ=II+(J-1)*NRTP
C              PSIREIN1(JJ)=DBLE(PSINIT(II,K))
C              PSIIMIN1(JJ)=DIMAG(PSINIT(II,K))
C            ENDDO
C          ENDDO
C        ENDIF
C        IF(MOD(JROT,2).EQ.1)THEN
C          J=0
C          DO K=2,2*JTOT+1,2
C            J=J+1
C            DO II=1,NRTP
C              JJ=II+(J-1)*NRTP
C              PSIREIN1(JJ)=DBLE(PSINIT(II,K))
C              PSIIMIN1(JJ)=DIMAG(PSINIT(II,K))
C            ENDDO
C          ENDDO
C        ENDIF
C      ELSE
C        J=0
C        DO K=1,2*JTOT+1,2
C          J=J+1
C          DO II=1,NRTP
C            JJ=II+(J-1)*NRTP
C            PSIREIN1(JJ)=DBLE(PSINIT(II,K))
C            PSIIMIN1(JJ)=DIMAG(PSINIT(II,K))
C          ENDDO
C        ENDDO
C        JTT=J
C      ENDIF
C      WRITE(166,*)
C      WRITE(166,*)'ORTHONORMALITY BETWEEN TWO WAVEFUNCTION'
C      WRITE(166,*)
C      SS2=0.0D0
C      DO L=1,JTT
C        SS1=0.0D0
C        DO II=1,NRTP
C          JJ=II+(L-1)*NRTP
C          SS1=SS1+(PSIREIN(JJ)*PSIREIN1(JJ)
C     1       +PSIIMIN(JJ)*PSIIMIN1(JJ))*DR*DT*DP
C        ENDDO
C        WRITE(166,'(I4,2X,F12.6)')L,SS1
C        SS2=SS2+SS1
C      ENDDO
C      WRITE(166,*)SS2
!
      IF(INA.EQ.1) GO TO 10
        WRITE(166,*)
        WRITE(166,*)'THE WORKABILITY OF PROJECTION TECHNIQUE'
        WRITE(166,*)
! 
!       A_k,mu and G_j,l,mu ARE SET TO ZERO
!       NOTE THE DIMENSION OF mu IS DEFINED BY THE HIGHEST JR(=NJ)
!       TO AVOID DYNAMICAL DIMENSION OF mu FOR DIFFERENT JR
! 
        IF(JT.NE.0) THEN
          DO JR=0,NJ
            DO L=1,JTT
              DO M=1,2*NJ+1
                AAP(JR,L,M)=DCMPLX(0.0D0,0.0D0)
              ENDDO
            ENDDO
            DO LP=0,NJ+JTOT
              DO M=1,2*NJ+1
                VALP(JR,LP,M)=0.0D0
              ENDDO
            ENDDO
          ENDDO
          DO JR=0,NJ
            DO L=1,JTT
              IF(MOD(JROT,2).EQ.0)THEN
                IL=-JT+2*L-2
              ENDIF
              IF(MOD(JROT,2).EQ.1)THEN
                IL=-JT+2*L-1
              ENDIF
              DO M=1,2*JR+1
                MP=M-JR-1
                CALL AKU(IL,MP,AA)
                AAP(JR,L,M)=AA
              ENDDO
            ENDDO
            DO LP=0,NJ+JTOT
              DO M=1,2*JR+1
                NP=ABS(M-JR-1)
                CALL CLD(JR,NP,LP,0,JT,NP,VAL)
                VALP(JR,LP,M)=VAL
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
      DDR=7.4D0/DFLOAT(N_R)
      DO KM=0,NV
        DO JR=0,NJ
          PROJ(KM,JR)=0.0D0
        ENDDO
      ENDDO
      DO I=1,N_R
        RSTAR=2.5D0+(I-1)*DDR+0.5D0*DDR
        CALL PSIJAC(JTT,LTOT,IC,DD,REQ,BETA,DM,EPS,RM12,
     1  RHO,THE,PHI,RSTAR,AAP,VALP,PSIREIN,PSIIMIN,PRE1,PIM1)
        DO KM=0,NV
          DO JR=0,NJ
            PROJ(KM,JR)=PROJ(KM,JR)+(PRE1(KM,JR)**2+PIM1(KM,JR)**2)*DDR
          ENDDO
        ENDDO
      ENDDO
      PR=0.0D0
      DO KM=0,NV
        DO JR=0,NJ
          WRITE(166,'(1X,2I4,2F12.6)')KM,JR,PROJ(KM,JR)
          PR=PR+PROJ(KM,JR)
        ENDDO
        WRITE(166,*)
      ENDDO
      WRITE(166,*)
      WRITE(166,'(2F12.6)')RSTAR,PR
      WRITE(166,*)
  10  CONTINUE
!************************************************************************************
!                                                                                   *
!     DEFINING THE THREE ARMS OF THE TRIATOM EMPLOYING HYPERSPHERICAL COORDINATES.  *
!     THE DOUBLET GROUND ADIABATIC PES.                                             *
!                                                                                   *
!************************************************************************************
      FACTR=1.D0/DSQRT(2.0D0)
      DO I=1,N_R
        RH=RHO(I)
        DO J=1,N_T
          TH=THE(J)
          JJ=J+(I-1)*N_T
          DO K=1,N_P
            PH=PHI(K)
            KK=K+(JJ-1)*N_P
            RNUC(1)=FACTR*RH*DM(3)*DSQRT(1.D0+DSIN(TH)*
     1              DCOS(PH-EPS(3)))
            RNUC(2)=FACTR*RH*DM(1)*DSQRT(1.D0+DSIN(TH)*DCOS(PH))
            RNUC(3)=FACTR*RH*DM(2)*DSQRT(1.D0+DSIN(TH)*
     1              DCOS(PH-EPS(2)))
            R1=RNUC(1)/AUANG
            R2=RNUC(2)/AUANG
            R3=RNUC(3)/AUANG
            CALL SINGLET(R1,R2,R3,POT)
            V11(KK)=POT(1,1)*AUEPS+30.833395D0
            V12(KK)=POT(1,2)*AUEPS+30.833395D0
            V13(KK)=POT(1,3)*AUEPS+30.833395D0
            V21(KK)=POT(2,1)*AUEPS+30.833395D0
            V22(KK)=POT(2,2)*AUEPS+30.833395D0
            V23(KK)=POT(2,3)*AUEPS+30.833395D0
            V31(KK)=POT(3,1)*AUEPS+30.833395D0
            V32(KK)=POT(3,2)*AUEPS+30.833395D0
            V33(KK)=POT(3,3)*AUEPS+30.833395D0
            IF(V11(KK).GT.4.0D0)THEN
              V11(KK)=0.0D0
            ENDIF
            IF(V12(KK).GT.4.0D0)THEN
              V12(KK)=0.0D0
            ENDIF
            IF(V13(KK).GT.4.0D0)THEN
              V13(KK)=0.0D0
            ENDIF
            IF(V21(KK).GT.4.0D0)THEN
              V21(KK)=0.0D0
            ENDIF
            IF(V22(KK).GT.4.0D0)THEN
              V22(KK)=0.0D0
            ENDIF
            IF(V23(KK).GT.4.0D0)THEN
              V23(KK)=0.0D0
            ENDIF
            IF(V31(KK).GT.4.0D0)THEN
              V31(KK)=0.0D0
            ENDIF
            IF(V32(KK).GT.4.0D0)THEN
              V32(KK)=0.0D0
            ENDIF
            IF(V33(KK).GT.4.0D0)THEN
              V33(KK)=0.0D0
            ENDIF
          ENDDO
        ENDDO
      ENDDO
CCCCCC
      NLAN=LAITER
      T=0.0D0
      ISTEP=0
      DDT=0.005D0
      KM=KMORSE
      JR=JROT
      IC=ICH
      CALL DEFAT(DDT,AT)
      CALL INITIAL_WEIGHT(EKIN,IC,DD,REQ,BETA,RM12,
     1KM,JR,RM123,AT,IDX,CKI,ET)
      WRITE(166,*)
      WRITE(166,*)'WEIGHT FACTORS FOR THE OUT GOING WAVE FUNCTIONS'
      WRITE(166,*)
      RSTAR=4.75D0
      DO IFCH=1,NCH
        RMCM=RM123(IFCH)
        CALL FINAL_WEIGHT(IFCH,DD,REQ,BETA,RM12,
     1  RM123,IDX,AT,AKF1,AKF2,AKF3,
     2  NNV1,NNV2,NNV3,NNJ1,NNJ2,NNJ3)
        DO IE=1,NE
          DO KM=0,NNV1(IE)
            DO JR=0,NNJ1(IE,KM)
              CKF1(IFCH,IE,KM,JR)=AKF1(IE,KM,JR)
              DO LP=0,NJ+JTOT
                XH=RSTAR*CKF1(IFCH,IE,KM,JR)*RMCM
                CALL HANKEL(LP,XH,HAN)
                HAN=1.0D0/HAN
                HANRE1(IFCH,IE,KM,JR,LP)=DBLE(HAN)
                HANIM1(IFCH,IE,KM,JR,LP)=DIMAG(HAN)
C               WRITE(101,'(5I4,2X,2F12.6)')IFH,IEE,KMM,JRR,LP,
C    1          HANRE(IFH,IEE,KMM,JRR,LP),
C    2          HANIM(IFH,IEE,KMM,JRR,LP)
              ENDDO
            ENDDO
          ENDDO
          DO KM=0,NNV2(IE)
            DO JR=0,NNJ2(IE,KM)
              CKF2(IFCH,IE,KM,JR)=AKF2(IE,KM,JR)
              DO LP=0,NJ+JTOT
                XH=RSTAR*CKF2(IFCH,IE,KM,JR)*RMCM
                CALL HANKEL(LP,XH,HAN)
                HAN=1.0D0/HAN
                HANRE2(IFCH,IE,KM,JR,LP)=DBLE(HAN)
                HANIM2(IFCH,IE,KM,JR,LP)=DIMAG(HAN)
              ENDDO
            ENDDO
          ENDDO
          DO KM=0,NNV3(IE)
            DO JR=0,NNJ3(IE,KM)
              CKF3(IFCH,IE,KM,JR)=AKF3(IE,KM,JR)*RMCM
              DO LP=0,NJ+JTOT
                XH=RSTAR*CKF3(IFCH,IE,KM,JR)
                CALL HANKEL(LP,XH,HAN)
                HAN=1.0D0/HAN
                HANRE3(IFCH,IE,KM,JR,LP)=DBLE(HAN)
                HANIM3(IFCH,IE,KM,JR,LP)=DIMAG(HAN)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CALL FLUSH(166)
!
      RHO1=1.5D0
      RHO2=2.5D0
      RHO3=3.5D0
      NRR(1)=INT((RHO1-RMIN)/DR)
      NRR(2)=INT((RHO2-RMIN)/DR)
      NRR(3)=INT((RHO3-RMIN)/DR)
      THE1=PI/6.0D0
      THE2=2.0D0*PI/6.0D0
      THE3=3.0D0*PI/6.0D0
      NTT(1)=INT(THE1/DT)
      NTT(2)=INT(THE2/DT)-1
      NTT(3)=INT(THE3/DT)-2
      PHI1=2.0D0*PI/3.0D0
      PHI2=4.0D0*PI/3.0D0
      PHI3=6.0D0*PI/3.0D0
      NPP(1)=INT(PHI1/DP)
      NPP(2)=INT(PHI2/DP)-1
      NPP(3)=INT(PHI3/DP)-2
!
!     THE FOURIER FREQUENCIES USED WHILE EVALUATING:  H |PSI>  (FFT)
!
      CALL DEFAK(DR,DT,DP,AKR,AKT,AKP)
      WRITE(21,*)'NORMALIZATION OF THE TOTAL WAVEFUNCTION AT EACH TIME'
      WRITE(21,*)
      WRITE(21,*)'AND TOTAL ENERGY OF THE SYSTEM AT EACH TIME'
      WRITE(21,*)
      CALL FLUSH(21)
      WRITE(13,*)
      WRITE(13,*)'WEIGHT FACTORS OF THREE PARTIAL WAVES AT EACH TIME'
      WRITE(13,*)
      CALL FLUSH(13)
CCCCCC
      DO I=1,N_R
        DO J=1,N_T
          II=J+(I-1)*N_T
          DO K=1,N_P
            JJ=K+(II-1)*N_P
            UU(1,1)=V11(JJ)
            UU(1,2)=V12(JJ)
            UU(1,3)=V13(JJ)
            UU(2,1)=V21(JJ)
            UU(2,2)=V22(JJ)
            UU(2,3)=V23(JJ)
            UU(3,1)=V31(JJ)
            UU(3,2)=V32(JJ)
            UU(3,3)=V33(JJ)
            MATZ=1
            CALL RS(NSTATE,NSTATE,UU,W,MATZ,ZS,FV1,FV2,IERR)
            IF(JJ.EQ.1)THEN
              DO I1=1,NSTATE
                DO J1=1,NSTATE
                  ZTEMP(I1,J1)=ZS(I1,J1)
                ENDDO
              ENDDO
            ENDIF
            IF(JJ.GT.1)THEN
              DO L1=1,NSTATE
                ZZZ(L1)=0.0D0 
                DO K1=1,NSTATE
                  ZZZ(L1)=ZZZ(L1)+ZTEMP(K1,L1)*ZS(K1,L1)
                ENDDO
                IF(ZZZ(L1).LT.0.0D0)THEN
                  DO K1=1,NSTATE
                    ZS(K1,L1)=-ZS(K1,L1) 
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
            DO I1=1,NSTATE
              DO J1=1,NSTATE
                ZTEMP(I1,J1)=ZS(I1,J1)
                ZP(I1,J1,JJ)=ZS(I1,J1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
      DO KK=1,JTT
        DO I=1,NRTP
          IS1=1
          II1=KK+(IS1-1)*JTT
          JJ1=I+(II1-1)*NRTP
          IS2=2
          II2=KK+(IS2-1)*JTT
          JJ2=I+(II2-1)*NRTP
          IS3=3
          II3=KK+(IS3-1)*JTT
          JJ3=I+(II3-1)*NRTP
          KK1=I+(KK-1)*NRTP
          PSIRE(JJ1)=ZP(1,1,I)*PSIREIN(KK1)
          PSIIM(JJ1)=ZP(1,1,I)*PSIIMIN(KK1)
          PSIRE(JJ2)=ZP(2,1,I)*PSIREIN(KK1)
          PSIIM(JJ2)=ZP(2,1,I)*PSIIMIN(KK1)
          PSIRE(JJ3)=ZP(3,1,I)*PSIREIN(KK1)
          PSIIM(JJ3)=ZP(3,1,I)*PSIIMIN(KK1)
        ENDDO
      ENDDO
CCCCCC
      SS1=0.0D0
      DO I=1,NSJRTP
        SS1=SS1+(PSIRE(I)*PSIRE(I)
     1     +PSIIM(I)*PSIIM(I))*VOLEM
      ENDDO
      WRITE(166,*)'here1',SS1
CCCCCC
      AJ_T=DFLOAT(JTOT)
      DO KK=1,JTT
        IF(JTOT.EQ.0)THEN
          AK(KK)=0.0D0
        ELSE
          IF(MOD(JROT,2).EQ.1)THEN
            AK(KK)=-DFLOAT(JTOT)+2.0D0*DFLOAT(KK)-1.0D0
          ENDIF
          IF(MOD(JROT,2).EQ.0)THEN
            AK(KK)=-DFLOAT(JTOT)+2.0D0*(DFLOAT(KK)-1.0D0)
          ENDIF
        ENDIF
        AKJ=AK(KK)
        AK_JM=ABS(AKJ-2)
        AK_JP=ABS(AKJ+2)
        IF(AK_JM.GT.AJ_T)THEN
          EKM(KK)=0.0D0
        ELSE
          EKM(KK)=DSQRT((AJ_T+AKJ)*(AJ_T-AKJ+1)
     1           *(AJ_T+AKJ-1)*(AJ_T-AKJ+2))
        ENDIF
        EK(KK)=AJ_T*(AJ_T+1.0D0)-AKJ*AKJ
        IF(AK_JP.GT.AJ_T)THEN
          EKP(KK)=0.0D0
        ELSE
          EKP(KK)=DSQRT((AJ_T-AKJ)*(AJ_T+AKJ+1)
     1           *(AJ_T-AKJ-1)*(AJ_T+AKJ+2))
        ENDIF
        WRITE(12,'(I3,5F12.6)')KK,AK_JM,AK_JP,EK(KK),EKM(KK),EKP(KK)
      CALL FLUSH(12)
      ENDDO
C      WRITE(12,'(I3,5F12.6)')
C      DO KK=1,JTT
C        IF(ABS(AK(KK)).GT.11)THEN
C          AK(KK)=0.0D0
C          EK(KK)=0.0D0
C          EKM(KK)=0.0D0
C          EKP(KK)=0.0D0
C        ENDIF
C        WRITE(12,'(I3,5F12.6)')KK,AK(KK),EK(KK),EKM(KK),EKP(KK)
C      ENDDO
      IST=0
      IST1=0
!
      CALL DFFTW_PLAN_DFT_3D(PLAN1,N_R,NT2,N_P,ARR1,ARR1,
     &                       FFTW_FORWARD,FFTW_MEASURE)
      CALL DFFTW_PLAN_DFT_3D(PLAN2,N_R,NT2,N_P,ARR,ARR,
     &                       FFTW_BACKWARD,FFTW_MEASURE)
      CALL DFFTW_PLAN_DFT_3D(PLAN3,N_R,NT2,N_P,ART,ART,
     &                       FFTW_BACKWARD,FFTW_MEASURE)
      CALL DFFTW_PLAN_DFT_3D(PLAN4,N_R,NT2,N_P,ARP,ARP,
     &                       FFTW_BACKWARD,FFTW_MEASURE)
      CALL DFFTW_PLAN_DFT_3D(PLAN5,N_R,NT2,N_P,ARP1,ARP1,
     &                       FFTW_BACKWARD,FFTW_MEASURE)
      CALL ENERGY(AK,EK,EKP,EKM,PSIRE,PSIIM,
     1AKR,AKT,AKP,WIN1,WIN2,RHOD,SINN,COSN,SIND,COSD,SIN2D,
     2V11,V12,V13,V21,V22,V23,V31,V32,V33,
     2EER1,EER2,EER3,EER4,EER5,EER6,EER7,EER8,EER9,EER10,EER,
     3RMST,SNG,VOLEM,EVEPS,
     4PLAN1,PLAN2,PLAN3,PLAN4,PLAN5)
      ENER1=0.0D0
      ENER2=0.0D0
      ENER3=0.0D0
      ENER4=0.0D0
      ENER5=0.0D0
      ENER6=0.0D0
      ENER7=0.0D0
      ENER8=0.0D0
      ENER9=0.0D0
      ENER10=0.0D0
      ENER=0.0D0
      DO KK=1,NSTATE*JTR
        ENER1=ENER1+EER1(KK)
        ENER2=ENER2+EER2(KK)
        ENER3=ENER3+EER3(KK)
        ENER4=ENER4+EER4(KK)
        ENER5=ENER5+EER5(KK)
        ENER6=ENER6+EER6(KK)
        ENER7=ENER7+EER7(KK)
        ENER8=ENER8+EER8(KK)
        ENER9=ENER9+EER9(KK)
        ENER10=ENER10+EER10(KK)
        ENER=ENER+EER(KK)
      ENDDO
      WRITE(21,'(12F12.8)')T,ENER1,ENER2,ENER3,ENER4,ENER5,ENER6,
     1ENER7,ENER8,ENER9,ENER10,ENER
      CALL FLUSH(21)
      IF(INA.EQ.1) GO TO 20
********************************************************************
      IF(INIT.EQ.1)THEN                                            !
        READ(2) ISTEP,IST,IST1,T                                   !
        READ(2) NLAN
        READ(2) (PSIRE(I),I=1,NSJRTP)                               ! 
        READ(2) (PSIIM(I),I=1,NSJRTP)                               ! 
      ENDIF                                                        !
********************************************************************
 50   ISTEP=ISTEP+1
      T=T+DDT
      CALL VPOTIM(RHO,VIM)
!
      CALL CPU_TIME(START1)
      OTIMES1=OMP_GET_WTIME()
      CALL ENERGY(AK,EK,EKP,EKM,PSIRE,PSIIM,
     1AKR,AKT,AKP,WIN1,WIN2,RHOD,SINN,COSN,SIND,COSD,SIN2D,
     2V11,V12,V13,V21,V22,V23,V31,V32,V33,
     2EER1,EER2,EER3,EER4,EER5,EER6,EER7,EER8,EER9,EER10,EER,
     3RMST,SNG,VOLEM,EVEPS,
     4PLAN1,PLAN2,PLAN3,PLAN4,PLAN5)
      OTIMEF1=OMP_GET_WTIME()
      ENTIME=OTIMEF1-OTIMES1
      CALL CPU_TIME(FINISH1)
      WRITE(111,*)'ENERGY',FINISH1-START1,ENTIME
      CALL FLUSH(111)
      ENER1=0.0D0
      ENER2=0.0D0
      ENER3=0.0D0
      ENER4=0.0D0
      ENER5=0.0D0
      ENER6=0.0D0
      ENER7=0.0D0
      ENER8=0.0D0
      ENER9=0.0D0
      ENER10=0.0D0
      ENER=0.0D0
      DO KK=1,NSTATE*JTR
        ENER1=ENER1+EER1(KK)
        ENER2=ENER2+EER2(KK)
        ENER3=ENER3+EER3(KK)
        ENER4=ENER4+EER4(KK)
        ENER5=ENER5+EER5(KK)
        ENER6=ENER6+EER6(KK)
        ENER7=ENER7+EER7(KK)
        ENER8=ENER8+EER8(KK)
        ENER9=ENER9+EER9(KK)
        ENER10=ENER10+EER10(KK)
        ENER=ENER+EER(KK)
      ENDDO
!
      WRITE(22,'(7F12.6)')T,ENER1,ENER2,ENER3,ENER4,ENER5
      WRITE(23,'(7F12.6)')T,ENER6,ENER7,ENER8,ENER9,ENER10
!
      CALL CPU_TIME(START2)
      OTIMES2=OMP_GET_WTIME()
!
      MMM=NLAN
      SCALEA=DSQRT(VOLEM)
      SCALEI=1.0D0/SCALEA
      DO I=1,NSJRTP
        U0A(I)=PSIRE(I)*SCALEA
        U0B(I)=PSIIM(I)*SCALEA
        U1A(I)=0.0D0
        U1B(I)=0.0D0
      ENDDO
!
      BN=0.0D0
      IT=0
   11 IT=IT+1
      DO I=1,NSJRTP
        WA(I,IT)=U0A(I)
        WB(I,IT)=U0B(I)
        PSIRE(I)=U0A(I)
        PSIIM(I)=U0B(I)
      ENDDO
!
      NC1=JTR
      DO IS=1,NSTATE
        DO JT=1,JTR
          II=JT+(IS-1)*JTR
          DO JJ=1,NRTP
            LL=JJ+(II-1)*NRTP
            PRE(JJ,JT,IS)=PSIRE(LL)
            PIM(JJ,JT,IS)=PSIIM(LL)
          ENDDO
        ENDDO
      ENDDO
!
      ICHUNK=1
!$OMP PARALLEL 
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(ICHUNK,NC1,AK,EK,EKP,EKM,PRE,PIM,PSIRE,PSIIM,
!$OMP$ AKR,AKT,AKP,WIN1,WIN2,RHOD,SINN,COSN,SIND,COSD,SIN2D,
!$OMP$ V11,V12,V13,V21,V22,V23,V31,V32,V33,
!$OMP$ RMST,RINV,SNG,VOLEM,
!$OMP$ PLAN1,PLAN2,PLAN3,PLAN4,PLAN5)
!$OMP DO SCHEDULE(DYNAMIC,ICHUNK)
      DO III=1,NSTATE*JTR
        IIT=III
        IS=(IIT-1)/DFLOAT(NC1)
        IIT=IIT-IS*NC1
        IS=IS+1
        JT=IIT
        II=JT+(IS-1)*JTR
        AK_J=AK(JT)
        EJK=EK(JT)
        EJKP=EKP(JT)
        EJKM=EKM(JT)
        IF(JT-1.LT.1)THEN
          S1=0.0D0
          DO JJ=1,NRTP
            X1(JJ)=0.0D0
            Y1(JJ)=0.0D0
            S1=S1+X1(JJ)**2+Y1(JJ)**2
          ENDDO
          S1=S1*VOLEM
        ELSE
          S1=0.0D0
          DO JJ=1,NRTP
            X1(JJ)=PRE(JJ,JT-1,IS)
            Y1(JJ)=PIM(JJ,JT-1,IS)
            S1=S1+X1(JJ)**2+Y1(JJ)**2
          ENDDO
          S1=S1*VOLEM
        ENDIF
        S2=0.0D0
        DO JJ=1,NRTP
          X2(JJ)=PRE(JJ,JT,IS)
          Y2(JJ)=PIM(JJ,JT,IS)
          S2=S2+X2(JJ)**2+Y2(JJ)**2
        ENDDO
        S2=S2*VOLEM
        IF(JT+1.GT.JTR)THEN
          S3=0.0D0
          DO JJ=1,NRTP
            X3(JJ)=0.0D0
            Y3(JJ)=0.0D0
            S3=S3+X3(JJ)**2+Y3(JJ)**2
          ENDDO
          S3=S3*VOLEM
        ELSE
          S3=0.0D0
          DO JJ=1,NRTP
            X3(JJ)=PRE(JJ,JT+1,IS)
            Y3(JJ)=PIM(JJ,JT+1,IS)
            S3=S3+X3(JJ)**2+Y3(JJ)**2
          ENDDO
          S3=S3*VOLEM
        ENDIF
!
        DO I=1,N_R
          DO J=1,N_T
            JJ=J+(I-1)*N_T
            DO K=1,N_P
              KK=K+(JJ-1)*N_P
              ARR1(I,J,K)=DCMPLX(X2(KK),Y2(KK))
              JP=NT2-J+1
              ARR1(I,JP,K)=-DCMPLX(X2(KK),Y2(KK))
            ENDDO
          ENDDO
        ENDDO
        CALL DFFTW_EXECUTE_DFT(PLAN1,ARR1,ARR1)
        DO I=1,N_R
          DO J=1,NT2
            DO K=1,N_P
              ARR(I,J,K)=ARR1(I,J,K)*AKR(I)*AKR(I)*SNG
              ART(I,J,K)=ARR1(I,J,K)*AKT(J)*AKT(J)*SNG
              ARP(I,J,K)=ARR1(I,J,K)*AKP(K)*AKP(K)*SNG
              ARP1(I,J,K)=-ARR1(I,J,K)*AKP(K)*SNG
            ENDDO
          ENDDO
        ENDDO
        CALL DFFTW_EXECUTE_DFT(PLAN2,ARR,ARR)
        CALL DFFTW_EXECUTE_DFT(PLAN3,ART,ART)
        CALL DFFTW_EXECUTE_DFT(PLAN4,ARP,ARP)
        CALL DFFTW_EXECUTE_DFT(PLAN5,ARP1,ARP1)
        DO I=1,N_R
          DO J=1,N_T
            JJ=J+(I-1)*N_T
            DO K=1,N_P
              KK=K+(JJ-1)*N_P
              LL=KK+(II-1)*NRTP
              HPSIRE1=0.5D0*RMST*DBLE(ARR(I,J,K))*SNG
              HPSIRE2=2.0D0*RMST*RHOD(I)*DBLE(ART(I,J,K))*SNG
              HPSIRE3=2.0D0*RMST*RHOD(I)*SIND(J)*DBLE(ARP(I,J,K))
     1               *WIN1(J)*SNG
              HPSIRE4=2.0D0*RMST*RHOD(I)*COSN(J)*SIND(J)
     1               *AK_J*DBLE(ARP1(I,J,K))*WIN1(J)*SNG
              HPSIRE5=0.5D0*RMST*RHOD(I)*SIND(J)*AK_J*AK_J
     1               *X2(KK)*WIN1(J)
              HPSIRE6=RMST*RHOD(I)*COSD(J)*EJK*X2(KK)*WIN2(J)
              HPSIRE7=-0.5D0*RMST*RHOD(I)
     1               *(0.25D0+4.0D0*SIN2D(J))*X2(KK)
              HPSIRE8=0.5D0*RMST*RHOD(I)*SINN(J)
     1               *COSD(J)*EJKM*X1(KK)*WIN2(J)
              HPSIRE9=0.5D0*RMST*RHOD(I)*SINN(J)
     1               *COSD(J)*EJKP*X3(KK)*WIN2(J)
              IF(IS.EQ.1)THEN
              HPSIRE10=V11(KK)*PRE(KK,JT,1)
     1                +V12(KK)*PRE(KK,JT,2)
     2                +V13(KK)*PRE(KK,JT,3)
              ENDIF
              IF(IS.EQ.2)THEN
              HPSIRE10=V21(KK)*PRE(KK,JT,1)
     1                +V22(KK)*PRE(KK,JT,2)
     2                +V23(KK)*PRE(KK,JT,3)
              ENDIF
              IF(IS.EQ.3)THEN
              HPSIRE10=V31(KK)*PRE(KK,JT,1)
     1                +V32(KK)*PRE(KK,JT,2)
     2                +V33(KK)*PRE(KK,JT,3)
              ENDIF
!
              HPSIIM1=0.5D0*RMST*DIMAG(ARR(I,J,K))*SNG
              HPSIIM2=2.0D0*RMST*RHOD(I)*DIMAG(ART(I,J,K))*SNG
              HPSIIM3=2.0D0*RMST*RHOD(I)*SIND(J)*DIMAG(ARP(I,J,K))
     1               *WIN1(J)*SNG
              HPSIIM4=2.0D0*RMST*RHOD(I)*COSN(J)*SIND(J)
     1               *AK_J*DIMAG(ARP1(I,J,K))*WIN1(J)*SNG
              HPSIIM5=0.5D0*RMST*RHOD(I)*SIND(J)*AK_J*AK_J
     1               *Y2(KK)*WIN1(J)
              HPSIIM6=RMST*RHOD(I)*COSD(J)*EJK*Y2(KK)*WIN2(J)
              HPSIIM7=-0.5D0*RMST*RHOD(I)
     1               *(0.25D0+4.0D0*SIN2D(J))*Y2(KK)
              HPSIIM8=0.5D0*RMST*RHOD(I)*SINN(J)
     1               *COSD(J)*EJKM*Y1(KK)*WIN2(J)
              HPSIIM9=0.5D0*RMST*RHOD(I)*SINN(J)
     1               *COSD(J)*EJKP*Y3(KK)*WIN2(J)
              IF(IS.EQ.1)THEN
              HPSIIM10=V11(KK)*PIM(KK,JT,1)
     1                +V12(KK)*PIM(KK,JT,2)
     2                +V13(KK)*PIM(KK,JT,3)
              ENDIF
              IF(IS.EQ.2)THEN
              HPSIIM10=V21(KK)*PIM(KK,JT,1)
     1                +V22(KK)*PIM(KK,JT,2)
     2                +V23(KK)*PIM(KK,JT,3)
              ENDIF
              IF(IS.EQ.3)THEN
              HPSIIM10=V31(KK)*PIM(KK,JT,1)
     1                +V32(KK)*PIM(KK,JT,2)
     2                +V33(KK)*PIM(KK,JT,3)
              ENDIF
!
              PSIRE(LL)=HPSIRE1+HPSIRE2+HPSIRE3+HPSIRE4+HPSIRE5
     1              +HPSIRE6+HPSIRE7+HPSIRE8+HPSIRE9+HPSIRE10
              PSIRE(LL)=PSIRE(LL)*RINV
              PSIIM(LL)=HPSIIM1+HPSIIM2+HPSIIM3+HPSIIM4+HPSIIM5
     1              +HPSIIM6+HPSIIM7+HPSIIM8+HPSIIM9+HPSIIM10
              PSIIM(LL)=PSIIM(LL)*RINV
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$OMP END DO 
!$OMP END PARALLEL
      IF(MMM.EQ.1) GO TO 1
      AN=0.0D0
      DO 3 I=1,NSJRTP
    3   AN=AN+U0A(I)*PSIRE(I)+U0B(I)*PSIIM(I)
        TT(IT,1)=AN
      DO I=1,NSJRTP
        PSIRE(I)=PSIRE(I)-AN*U0A(I)-BN*U1A(I)
        PSIIM(I)=PSIIM(I)-AN*U0B(I)-BN*U1B(I)
      ENDDO
      BN=0.0D0
      DO I=1,NSJRTP
        BN=BN+PSIRE(I)*PSIRE(I)+PSIIM(I)*PSIIM(I)
      ENDDO
      BN=DSQRT(BN)
      BN1=1.0D0/BN
      DO I=1,NSJRTP
        U1A(I)=U0A(I)
        U1B(I)=U0B(I)
        U0A(I)=PSIRE(I)*BN1
        U0B(I)=PSIIM(I)*BN1
      ENDDO
      TT(IT,2)=BN
      IF(DABS(BN).LT.1.D-6)THEN
        WRITE(166,'(I4,E10.3)')'BN',IT,BN
      ENDIF
      IF(IT.LT.MMM) GO TO 11
    1 CONTINUE
!
      IF(MMM.EQ.1) GO TO 2
      DO I=1,NSJRTP
        PSIRE(I)=0.0D0
        PSIIM(I)=0.0D0
      ENDDO
      DO I=2,MMM
        D(I)=TT(I,1)
        E(I)=TT(I-1,2)
      ENDDO
      D(1)=TT(1,1)
      DO I=1,MMM
        DO J=1,MMM
          Z(I,J)=0.0D0
          IF(I.EQ.J) Z(I,J)=1.0D0
        ENDDO
      ENDDO
      CALL TQLI(D,E,MMM,Z)
      DO I=1,MMM
        EE(I)=D(I)*HBAR
      ENDDO
      CALL ODR(MMM,EE)
      DO I=1,MMM
        WRITE(9,'(2I6,F16.6)')ISTEP,I,EE(I)
        CALL FLUSH(9)
      ENDDO
      DO I=1,MMM
        SUA=0.0D0
        SUB=0.0D0
        DO J=1,MMM
          SUA=SUA+Z(I,J)*Z(1,J)*DCOS(D(J)*DDT)
          SUB=SUB-Z(I,J)*Z(1,J)*DSIN(D(J)*DDT)
        ENDDO
        RA(I)=SUA*SCALEI
        RB(I)=SUB*SCALEI
      ENDDO
***************************************************
*     AVERAGE OF THE FIVE LAST COMPONENTS         *
***************************************************
      RAV=0.0D0
      DO I=MMM-4,MMM
        RAV=RAV+RA(I)*RA(I)+RB(I)*RB(I)
      ENDDO
      RAV=RAV*0.2D0
      DO J=1,MMM
        DO I=1,NSJRTP
          PSIRE(I)=PSIRE(I)+WA(I,J)*RA(J)-WB(I,J)*RB(J)
          PSIIM(I)=PSIIM(I)+WA(I,J)*RB(J)+WB(I,J)*RA(J)
        ENDDO
      ENDDO
    2 CONTINUE
!
      OTIMEF2=OMP_GET_WTIME()
      ALNTIME=OTIMEF2-OTIMES2
      CALL CPU_TIME(FINISH2)
      WRITE(111,*)'LANCZOS',FINISH2-START2,ALNTIME
      CALL FLUSH(111)
!
      RAV=RAV*VOLEM
      IF(EPSL1.GT.0.0D0)THEN
        IF(RAV.LT.EPSL1) NLAN=NLAN-NLDEL
        IF(RAV.GT.EPSL2) NLAN=NLAN+NLDEL
        IF(NLAN.GT.LANMAX) NLAN=LANMAX
        IF(NLAN.LT.LANMIN) NLAN=LANMIN
      ENDIF
      WRITE(14,*)ISTEP,RAV,NLAN
      CALL FLUSH(14)
!
      SS1=0.0D0
      DO IS=1,NSTATE
        SS(IS)=0.0D0
        RHAV(IS)=0.0D0
        DO IJT=1,JTR
          II=IJT+(IS-1)*JTR
          DO IR=1,N_R
            DO IT=1,N_T
              IRT=IT+(IR-1)*N_T
              DO IP=1,N_P
                IRTP=IP+(IRT-1)*N_P
                JJ=IRTP+(II-1)*NRTP
C                PSI(JJ)=PSI(JJ)
C     1                 +PSIRE(JJ)**2+PSIIM(JJ)**2
                SS(IS)=SS(IS)+(PSIRE(JJ)**2+PSIIM(JJ)**2)*VOLEM
                RHAV(IS)=RHAV(IS)
     1                  +RHO(IR)*(PSIRE(JJ)**2+PSIIM(JJ)**2)*VOLEM
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        SS1=SS1+SS(IS)
      ENDDO
      RHAVT=RHAV(1)+RHAV(2)+RHAV(3)
      WRITE(21,'(8F16.8)')T,SS(1),SS(2),SS(3),SS1,ENER,RHAVT,RHAVT/SS1
!
      CALL FLUSH(21)
      CALL FLUSH(22)
      CALL FLUSH(23)
!
C      DO IS=1,NSTATE
C        DO IRTP=1,NRTP
C          LL=IRTP+(IS-1)*NRTP
C          PSI(LL)=0.0D0
C          DO IJT=1,JTR
C            II=IJT+(IS-1)*JTR
C            JJ=IRTP+(II-1)*NRTP
C            PSI(LL)=PSI(LL)+PSIRE(JJ)**2+PSIIM(JJ)**2
C          ENDDO
C        ENDDO
C      ENDDO
C!
C      IF(MOD(ISTEP-1,100).EQ.0) THEN
C        IST=IST+1
C        DO II=1,3
C          ISTT=II+(IST-1)*3
C          WRITE(199+ISTT,*)NRR(II)*DR+RMIN,T 
C          CALL FLUSH(199+ISTT)
C          I=NRR(II)
C          DO J=1,N_T
C            DO K=1,N_P
C              JJ=J+(I-1)*N_T
C              KK=K+(JJ-1)*N_P
C              WRITE(199+ISTT,'(7F12.6)')THE(J),PHI(K),
C     1        PSI(KK),PSI(KK+NRTP),PSI(KK+NRTP+NRTP)
C              CALL FLUSH(199+ISTT)
C            ENDDO
C            WRITE(199+ISTT,*) 
C            CALL FLUSH(199+ISTT)
C          ENDDO
C          WRITE(299+ISTT,*)NTT(II)*DT,T
C          CALL FLUSH(299+ISTT)
C          J=NTT(II)
C          DO I=1,N_R
C            DO K=1,N_P
C              JJ=J+(I-1)*N_T
C              KK=K+(JJ-1)*N_P
C              WRITE(299+ISTT,'(7F12.6)')RHO(I),PHI(K),
C     1        PSI(KK),PSI(KK+NRTP),PSI(KK+NRTP+NRTP)
C              CALL FLUSH(299+ISTT)
C            ENDDO
C            WRITE(299+ISTT,*)
C            CALL FLUSH(299+ISTT)
C          ENDDO
C          WRITE(399+ISTT,*)NPP(II)*DP,T
C          CALL FLUSH(399+ISTT)
C          K=NPP(II)
C          DO I=1,N_R
C            DO J=1,N_T
C              JJ=J+(I-1)*N_T
C              KK=K+(JJ-1)*N_P
C              WRITE(399+ISTT,'(7F12.6)')RHO(I),THE(J),
C     1        PSI(KK),PSI(KK+NRTP),PSI(KK+NRTP+NRTP)
C              CALL FLUSH(399+ISTT)
C            ENDDO
C            WRITE(399+ISTT,*)
C            CALL FLUSH(399+ISTT)
C          ENDDO
C        ENDDO
C      ENDIF
!
C      IF(MOD(ISTEP-1,100).EQ.0)THEN
C        IST1=IST1+1
C        CALL PLOT(IST1,RHO,THE,PHI,PSI)
C      ENDIF
!
      CALL CPU_TIME(START3)
      OTIMES3=OMP_GET_WTIME()
      ICHUNK=1
!$OMP PARALLEL
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(ICHUNK,JTT,DD,REQ,BETA,DM,EPS,RM12,RHO,THE,
!$OMP$ PHI,RSTAR,AAP,VALP,ZP,PSIRE,PSIIM,ISTEP)
!$OMP DO SCHEDULE(DYNAMIC,ICHUNK)
      DO LP=0,NJ+JTOT
        DO IFCH=1,NCH
          CALL PSIJAC1(JTT,LP,IFCH,DD,REQ,BETA,DM,EPS,RM12,
     1    RHO,THE,PHI,RSTAR,AAP,VALP,ZP,PSIRE,PSIIM,
     2    PRE1,PIM1,PRE2,PIM2,PRE3,PIM3)
          DO KM=0,NV
            DO JR=0,NJ
              QRE1(KM,JR,IFCH)=PRE1(KM,JR)
              QIM1(KM,JR,IFCH)=PIM1(KM,JR)
              QRE2(KM,JR,IFCH)=PRE2(KM,JR)
              QIM2(KM,JR,IFCH)=PIM2(KM,JR)
              QRE3(KM,JR,IFCH)=PRE3(KM,JR)
              QIM3(KM,JR,IFCH)=PIM3(KM,JR)
            ENDDO
          ENDDO
        ENDDO
        DO KM=0,NV
          DO JR=0,NJ
            WRITE(24+LP,77)ISTEP,KM,JR,(QRE1(KM,JR,IFCH),IFCH=1,NCH),
     1                                 (QRE2(KM,JR,IFCH),IFCH=1,NCH),
     2                                 (QRE3(KM,JR,IFCH),IFCH=1,NCH)
            WRITE(24+LP,77)ISTEP,KM,JR,(QIM1(KM,JR,IFCH),IFCH=1,NCH),
     1                                 (QIM2(KM,JR,IFCH),IFCH=1,NCH),
     2                                 (QIM3(KM,JR,IFCH),IFCH=1,NCH)
            CALL FLUSH(24+LP)
          ENDDO
          WRITE(24+LP,*)
          WRITE(24+LP,*)
          CALL FLUSH(24+LP)
        ENDDO
        WRITE(24+LP,*)
        WRITE(24+LP,*)
        WRITE(24+LP,*)
        WRITE(24+LP,*)
        CALL FLUSH(24+LP)
      ENDDO
!$OMP END DO 
!$OMP END PARALLEL
      OTIMEF3=OMP_GET_WTIME()
      TRTIME=OTIMEF3-OTIMES3
      CALL CPU_TIME(FINISH3)
      WRITE(111,*)'TRANSITION AMPLITUDE',FINISH3-START3,TRTIME
      CALL FLUSH(111)
!
      DO IS=1,NSTATE
        DO IJT=1,JTR
          II=IJT+(IS-1)*JTR
          DO IR=1,N_R
            DO IT=1,N_T
              IRT=IT+(IR-1)*N_T
              DO IP=1,N_P
                IRTP=IP+(IRT-1)*N_P
                LL=IRTP+(II-1)*NRTP
                PSIRE(LL)=PSIRE(LL)*VIM(IR)
                PSIIM(LL)=PSIIM(LL)*VIM(IR)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
****************************************************************
      IF(MOD(ISTEP,100).EQ.0)THEN                              !
        WRITE(2)ISTEP,IST,IST1,T                               !
        WRITE(2)NLAN                                           !
        WRITE(2)(PSIRE(I),I=1,NSJRTP)                          !
        WRITE(2)(PSIIM(I),I=1,NSJRTP)                          !
        CALL FLUSH(2)                                          !
        REWIND(2)                                              !
      ENDIF                                                    !
****************************************************************
      IF (ISTEP.GE.NSTEP) GO TO 40
      GO TO 50
  40  CONTINUE
  20  CONTINUE
      CALL DFFTW_DESTROY_PLAN(PLAN1)
      CALL DFFTW_DESTROY_PLAN(PLAN2)
      CALL DFFTW_DESTROY_PLAN(PLAN3)
      CALL DFFTW_DESTROY_PLAN(PLAN4)
      CALL DFFTW_DESTROY_PLAN(PLAN5)
!
      IF(INA.EQ.0) GO TO 60
!*****
      CALL DFFTW_PLAN_DFT_1D(PLAN6,NSTP,ARR6,ARR6,
     &                       FFTW_FORWARD,FFTW_MEASURE)
      CALL DFFTW_PLAN_DFT_1D(PLAN7,NSTP,ARR7,ARR7,
     &                       FFTW_FORWARD,FFTW_MEASURE)
      CALL DFFTW_PLAN_DFT_1D(PLAN8,NSTP,ARR8,ARR8,
     &                       FFTW_FORWARD,FFTW_MEASURE)
      SN_T=1.0D0/DSQRT(DFLOAT(NSTP))
!*****
      DO LP=0,NJ+JTOT
        DO IS=1,NSTEP
          DO KM=0,NV
            DO JR=0,NJ
              READ(24+LP,77)ISTP,K,J,
     1        AAR1,AAR2,AAR3,BBR1,BBR2,BBR3,CCR1,CCR2,CCR3
              READ(24+LP,77)ISTP,K,J,
     1        AAI1,AAI2,AAI3,BBI1,BBI2,BBI3,CCI1,CCI2,CCI3
              URE1(IS,KM,JR,1,LP)=AAR1
              URE1(IS,KM,JR,2,LP)=AAR2
              URE1(IS,KM,JR,3,LP)=AAR3
              UIM1(IS,KM,JR,1,LP)=AAI1
              UIM1(IS,KM,JR,2,LP)=AAI2
              UIM1(IS,KM,JR,3,LP)=AAI3
!
              URE2(IS,KM,JR,1,LP)=BBR1
              URE2(IS,KM,JR,2,LP)=BBR2
              URE2(IS,KM,JR,3,LP)=BBR3
              UIM2(IS,KM,JR,1,LP)=BBI1
              UIM2(IS,KM,JR,2,LP)=BBI2
              UIM2(IS,KM,JR,3,LP)=BBI3
!
              URE3(IS,KM,JR,1,LP)=CCR1
              URE3(IS,KM,JR,2,LP)=CCR2
              URE3(IS,KM,JR,3,LP)=CCR3
              UIM3(IS,KM,JR,1,LP)=CCI1
              UIM3(IS,KM,JR,2,LP)=CCI2
              UIM3(IS,KM,JR,3,LP)=CCI3
            ENDDO
            READ(24+LP,*)
            READ(24+LP,*)
          ENDDO
          READ(24+LP,*)
          READ(24+LP,*)
          READ(24+LP,*)
          READ(24+LP,*)
        ENDDO
      ENDDO
!*****
      EFFPOT=DFLOAT(JTOT)*(DFLOAT(JTOT)+1.0D0)*HBAR*HBAR/EMU/RSTAR/RSTAR
      DO LP=0,NJ+JTOT
        DO IS=1,NSTEP
          DO KM=0,NV
            DO JR=0,NJ
              DO IFCH=1,NCH
!
                URE11(IS,KM,JR,IFCH,LP)=
     1          URE1(IS,KM,JR,IFCH,LP)*DCOS(EFFPOT*DDT*IS/HBAR)
     2         +UIM1(IS,KM,JR,IFCH,LP)*DSIN(EFFPOT*DDT*IS/HBAR)
                UIM11(IS,KM,JR,IFCH,LP)=
     1          UIM1(IS,KM,JR,IFCH,LP)*DCOS(EFFPOT*DDT*IS/HBAR)
     2         -URE1(IS,KM,JR,IFCH,LP)*DSIN(EFFPOT*DDT*IS/HBAR)
!
                URE21(IS,KM,JR,IFCH,LP)=
     1          URE2(IS,KM,JR,IFCH,LP)*DCOS(EFFPOT*DDT*IS/HBAR)
     2         +UIM2(IS,KM,JR,IFCH,LP)*DSIN(EFFPOT*DDT*IS/HBAR)
                UIM21(IS,KM,JR,IFCH,LP)=
     1          UIM2(IS,KM,JR,IFCH,LP)*DCOS(EFFPOT*DDT*IS/HBAR)
     2         -URE2(IS,KM,JR,IFCH,LP)*DSIN(EFFPOT*DDT*IS/HBAR)
!
                URE31(IS,KM,JR,IFCH,LP)=
     1          URE3(IS,KM,JR,IFCH,LP)*DCOS(EFFPOT*DDT*IS/HBAR)
     2         +UIM3(IS,KM,JR,IFCH,LP)*DSIN(EFFPOT*DDT*IS/HBAR)
                UIM31(IS,KM,JR,IFCH,LP)=
     1          UIM3(IS,KM,JR,IFCH,LP)*DCOS(EFFPOT*DDT*IS/HBAR)
     2         -URE3(IS,KM,JR,IFCH,LP)*DSIN(EFFPOT*DDT*IS/HBAR)
!
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!*****
      DO LP=0,NJ+JTOT
        DO IFCH=1,NCH
          DO KM=0,NV
            DO JR=0,NJ
              ISS=0
              DO IS=NSTART+1,NSTP+NSTART
                ISS=ISS+1
                IF(IS.LE.NSTEP)THEN
                  ARR6(ISS)=DCMPLX(URE11(IS,KM,JR,IFCH,LP)
     1                            ,UIM11(IS,KM,JR,IFCH,LP))
                  ARR7(ISS)=DCMPLX(URE21(IS,KM,JR,IFCH,LP)
     1                            ,UIM21(IS,KM,JR,IFCH,LP))
                  ARR8(ISS)=DCMPLX(URE31(IS,KM,JR,IFCH,LP)
     1                            ,UIM31(IS,KM,JR,IFCH,LP))
                ELSE
                  ARR6(ISS)=DCMPLX(0.0D0,0.0D0)
                  ARR7(ISS)=DCMPLX(0.0D0,0.0D0)
                  ARR8(ISS)=DCMPLX(0.0D0,0.0D0)
                ENDIF
              ENDDO
              WRITE(166,*)LP,IFCH,KM,JR,ISS
              CALL DFFTW_EXECUTE_DFT(PLAN6,ARR6,ARR6)
              CALL DFFTW_EXECUTE_DFT(PLAN7,ARR7,ARR7)
              CALL DFFTW_EXECUTE_DFT(PLAN8,ARR8,ARR8)
              DO IS=1,NSTP
                URRE1(IS,KM,JR,IFCH,LP)=DBLE(ARR6(IS))*SN_T
                URIM1(IS,KM,JR,IFCH,LP)=DIMAG(ARR6(IS))*SN_T
                URRE2(IS,KM,JR,IFCH,LP)=DBLE(ARR7(IS))*SN_T
                URIM2(IS,KM,JR,IFCH,LP)=DIMAG(ARR7(IS))*SN_T
                URRE3(IS,KM,JR,IFCH,LP)=DBLE(ARR8(IS))*SN_T
                URIM3(IS,KM,JR,IFCH,LP)=DIMAG(ARR8(IS))*SN_T
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!*****
      CALL VIBMATRC(RM12,VALRC11,VALRC22)
      CALL VIBMATINF(RM12,VALINF12,VALINF21)
!*****
      DO IFCH=1,NCH
        IS=IDX-1
        DO IE=1,NE
          IS=IS+1
          DO KM=0,NV
            DO JR=0,NJ
              ARE1(IS,KM,JR,IFCH)=0.0D0
              ARE2(IS,KM,JR,IFCH)=0.0D0
              ARE3(IS,KM,JR,IFCH)=0.0D0
              AIM1(IS,KM,JR,IFCH)=0.0D0
              AIM2(IS,KM,JR,IFCH)=0.0D0
              AIM3(IS,KM,JR,IFCH)=0.0D0
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO LP=0,NJ+JTOT
        DO IFCH=1,NCH
          IS=IDX-1
          DO IE=1,NE
            IS=IS+1
            DO KM=0,NNV1(IE)
              DO JR=0,NNJ1(IE,KM)
                ARE1(IS,KM,JR,IFCH)=URRE1(IS,KM,JR,IFCH,LP)
     1                             *HANRE1(IFCH,IE,KM,JR,LP)
     1                             -URIM1(IS,KM,JR,IFCH,LP)
     1                             *HANIM1(IFCH,IE,KM,JR,LP)
                AIM1(IS,KM,JR,IFCH)=URRE1(IS,KM,JR,IFCH,LP)
     1                             *HANIM1(IFCH,IE,KM,JR,LP)
     1                             +URIM1(IS,KM,JR,IFCH,LP)
     1                             *HANRE1(IFCH,IE,KM,JR,LP)
              ENDDO
            ENDDO
            DO KM=0,NNV2(IE)
              DO JR=0,NNJ2(IE,KM)
                ARE2(IS,KM,JR,IFCH)=URRE2(IS,KM,JR,IFCH,LP)
     1                             *HANRE2(IFCH,IE,KM,JR,LP)
     1                             -URIM2(IS,KM,JR,IFCH,LP)
     1                             *HANIM2(IFCH,IE,KM,JR,LP)
                AIM2(IS,KM,JR,IFCH)=URRE2(IS,KM,JR,IFCH,LP)
     1                             *HANIM2(IFCH,IE,KM,JR,LP)
     1                             +URIM2(IS,KM,JR,IFCH,LP)
     1                             *HANRE2(IFCH,IE,KM,JR,LP)
              ENDDO
            ENDDO
            DO KM=0,NNV3(IE)
              DO JR=0,NNJ3(IE,KM)
                ARE3(IS,KM,JR,IFCH)=URRE3(IS,KM,JR,IFCH,LP)
     1                             *HANRE3(IFCH,IE,KM,JR,LP)
     1                             -URIM3(IS,KM,JR,IFCH,LP)
     1                             *HANIM3(IFCH,IE,KM,JR,LP)
                AIM3(IS,KM,JR,IFCH)=URRE3(IS,KM,JR,IFCH,LP)
     1                             *HANIM3(IFCH,IE,KM,JR,LP)
     1                             +URIM3(IS,KM,JR,IFCH,LP)
     1                             *HANRE3(IFCH,IE,KM,JR,LP)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO IFCH=1,NCH
          IS=IDX-1
          DO IE=1,NE
            IS=IS+1
            DO KM=0,NV
              DO JR=0,NJ
                DSRE1=0.0D0
                DSIM1=0.0D0
                DSRE2=0.0D0
                DSIM2=0.0D0
                DO KM1=0,NV
                  DSRE1=DSRE1+ARE1(IS,KM1,JR,IFCH)*VALRC11(IFCH,KM1,KM)
     1                       +ARE2(IS,KM1,JR,IFCH)*VALINF21(IFCH,KM1,KM)
                  DSIM1=DSIM1+AIM1(IS,KM1,JR,IFCH)*VALRC11(IFCH,KM1,KM)
     1                       +AIM2(IS,KM1,JR,IFCH)*VALINF21(IFCH,KM1,KM)
                  DSRE2=DSRE2+ARE2(IS,KM1,JR,IFCH)*VALRC22(IFCH,KM1,KM)
     1                       +ARE1(IS,KM1,JR,IFCH)*VALINF12(IFCH,KM1,KM)
                  DSIM2=DSIM2+AIM2(IS,KM1,JR,IFCH)*VALRC22(IFCH,KM1,KM)
     1                       +AIM1(IS,KM1,JR,IFCH)*VALINF12(IFCH,KM1,KM)
                ENDDO
                URR1(IE,KM,JR,IFCH,LP)=
     1             +(CKF1(IFCH,IE,KM,JR)/CKI(IE))*(DSRE1**2+DSIM1**2)
                URR2(IE,KM,JR,IFCH,LP)=
     1             +(CKF2(IFCH,IE,KM,JR)/CKI(IE))*(DSRE2**2+DSIM2**2)
                URR3(IE,KM,JR,IFCH,LP)=
     1             +(CKF3(IFCH,IE,KM,JR)/CKI(IE))
     2             *(ARE3(IS,KM,JR,IFCH)**2+AIM3(IS,KM,JR,IFCH)**2)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!*****
      DO IE=1,NE
        SE(IE)=0.0D0
        DO LP=0,NJ+JTOT
          DO IFCH=1,NCH
            DO KM=0,NV
              DO JR=0,NJ
                SE(IE)=SE(IE)+URR1(IE,KM,JR,IFCH,LP)
     1                       +URR2(IE,KM,JR,IFCH,LP)
     2                       +URR3(IE,KM,JR,IFCH,LP)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO LP=0,NJ+JTOT
          DO IFCH=1,NCH
            DO KM=0,NV
              DO JR=0,NJ
                URR1(IE,KM,JR,IFCH,LP)=URR1(IE,KM,JR,IFCH,LP)/SE(IE)
                URR2(IE,KM,JR,IFCH,LP)=URR2(IE,KM,JR,IFCH,LP)/SE(IE)
                URR3(IE,KM,JR,IFCH,LP)=URR3(IE,KM,JR,IFCH,LP)/SE(IE)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO IFCH=1,NCH
          DO KM=0,NV
            DO JR=0,NJ
              SUR1(IE,KM,JR,IFCH)=0.0D0
              SUR2(IE,KM,JR,IFCH)=0.0D0
              SUR3(IE,KM,JR,IFCH)=0.0D0
            ENDDO
          ENDDO
        ENDDO
        DO IFCH=1,NCH
          DO KM=0,NV
            DO JR=0,NJ
              DO LP=0,NJ+JTOT
                SUR1(IE,KM,JR,IFCH)=SUR1(IE,KM,JR,IFCH)
     1                             +URR1(IE,KM,JR,IFCH,LP)
                SUR2(IE,KM,JR,IFCH)=SUR2(IE,KM,JR,IFCH)
     1                             +URR2(IE,KM,JR,IFCH,LP)
                SUR3(IE,KM,JR,IFCH)=SUR3(IE,KM,JR,IFCH)
     1                             +URR3(IE,KM,JR,IFCH,LP)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!*****
      WRITE(16,*)
      WRITE(16,*)'TRANSITION PROBABILITY AS FUNCTIONS OF INITIAL'
      WRITE(16,*)'KE OF THE ATOM AND FINAL VIB. AND ROT. STATES'
      WRITE(16,*)'OF THE DIATOM'
      WRITE(16,*)
      AMINUS=ENER-(EKIN+EVJ)
      DO IE=1,NE
        ET(IE)=ET(IE)/EVEPS-AMINUS-EVJ
        DO KM=0,NV
          DO JR=0,NJ
            WRITE(16,88)ET(IE),KM,JR,
     1                 (SUR1(IE,KM,JR,IFCH),IFCH=1,NCH),
     2                  SUR1(IE,KM,JR,ICH+1)+SUR1(IE,KM,JR,ICH+2),
     3                 (SUR2(IE,KM,JR,IFCH),IFCH=1,NCH),
     4                  SUR2(IE,KM,JR,ICH+1)+SUR2(IE,KM,JR,ICH+2),
     5                 (SUR3(IE,KM,JR,IFCH),IFCH=1,NCH),
     6                  SUR3(IE,KM,JR,ICH+1)+SUR3(IE,KM,JR,ICH+2)
          ENDDO
          WRITE(16,*)
        ENDDO
        WRITE(16,*)
        WRITE(16,*)
      ENDDO
!
      DO IE=1,NE
        DO IFCH=1,NCH
          DO KM=0,NV
            RR1(IE,KM,IFCH)=0.0D0
            RR2(IE,KM,IFCH)=0.0D0
            RR3(IE,KM,IFCH)=0.0D0
            DO JR=0,NJ
              RR1(IE,KM,IFCH)=RR1(IE,KM,IFCH)+SUR1(IE,KM,JR,IFCH)
              RR2(IE,KM,IFCH)=RR2(IE,KM,IFCH)+SUR2(IE,KM,JR,IFCH)
              RR3(IE,KM,IFCH)=RR3(IE,KM,IFCH)+SUR3(IE,KM,JR,IFCH)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!
      WRITE(17,*)
      WRITE(17,*)'TRANSITION PROBABILITY AS FUNCTIONS OF INITIAL KE'
      WRITE(17,*)'OF THE ATOM AND FINAL VIB. STATES OF THE DIATOM'
      WRITE(17,*)
!
      WRITE(17,*)
      DO IE=1,NE
        DO KM=0,NV
          WRITE(17,99)ET(IE),KM,
     1                (RR1(IE,KM,IFCH),IFCH=1,NCH),
     2                RR1(IE,KM,ICH+1)+RR1(IE,KM,ICH+2),
     3                (RR2(IE,KM,IFCH),IFCH=1,NCH),
     4                RR2(IE,KM,ICH+1)+RR2(IE,KM,ICH+2),
     5                (RR3(IE,KM,IFCH),IFCH=1,NCH),
     6                RR3(IE,KM,ICH+1)+RR3(IE,KM,ICH+2)
        ENDDO
        WRITE(17,*)
      ENDDO
      WRITE(18,*)'TRANSITION PROBABILITY AS FUNCTIONS OF INITIAL KE'
      WRITE(19,*)'TRANSITION PROBABILITY AS FUNCTIONS OF INITIAL KE'
      WRITE(20,*)'TRANSITION PROBABILITY AS FUNCTIONS OF INITIAL KE'
      WRITE(18,*)'FINAL VIBRATIONAL STATE,v = 0'
      WRITE(19,*)'FINAL VIBRATIONAL STATE,v = 1'
      WRITE(20,*)'FINAL VIBRATIONAL STATE,v = 2'
      WRITE(18,*)
      WRITE(19,*)
      WRITE(20,*)
      DO KM=0,NV
        DO IE=1,NE
         WRITE(118+KM,99)ET(IE),KM,
     1   (RR1(IE,KM,IFCH),IFCH=1,NCH),RR1(IE,KM,ICH+1)+RR1(IE,KM,ICH+2),
     2   (RR2(IE,KM,IFCH),IFCH=1,NCH),RR2(IE,KM,ICH+1)+RR2(IE,KM,ICH+2),
     3   (RR3(IE,KM,IFCH),IFCH=1,NCH),RR3(IE,KM,ICH+1)+RR3(IE,KM,ICH+2)
        ENDDO
      ENDDO
      DO IE=1,NE
        SUM1=0.0D0
        SUM2=0.0D0
        SUM3=0.0D0
        SUM4=0.0D0
        SUM5=0.0D0
        SUM6=0.0D0
        DO KM=0,NV
          SUM1=SUM1+RR1(IE,KM,ICH)
          SUM2=SUM2+RR1(IE,KM,ICH+1)+RR1(IE,KM,ICH+2)
          SUM3=SUM3+RR2(IE,KM,ICH)
          SUM4=SUM4+RR2(IE,KM,ICH+1)+RR2(IE,KM,ICH+2)
          SUM5=SUM5+RR3(IE,KM,ICH)
          SUM6=SUM6+RR3(IE,KM,ICH+1)+RR3(IE,KM,ICH+2)
        ENDDO
        WRITE(117,'(1X,8F12.6)')ET(IE),
     1                          SUM1,SUM2,SUM3+SUM5,SUM4+SUM6
      ENDDO
      CALL DFFTW_DESTROY_PLAN(PLAN6)
      CALL DFFTW_DESTROY_PLAN(PLAN7)
      CALL DFFTW_DESTROY_PLAN(PLAN8)
  60  CONTINUE
  44  FORMAT(1X,12F12.6)
  55  FORMAT(1X,2I4,1X,4F12.6)
  66  FORMAT(1X,I4,2X,F12.6)
  77  FORMAT(1X,3I4,9F12.8)
  88  FORMAT(1X,F12.6,2I4,2X,12F12.6)
  99  FORMAT(1X,F8.4,I4,12F10.6)
      STOP
      END
!********************************************************
!                                                       *
!     INITIALIZATION OF THE WAVEFUNCTION                *
!                                                       *
!********************************************************
      SUBROUTINE PSIHYP(KMORSE,JROT,ICH,EKIN,RMIN,RMAX,
     1DD,REQ,BETA,DM,EPS,RM12,RM123,PSINIT)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,A1,AA,S3
      COMPLEX*16 XI,FA,SS
      COMPLEX*16 PSINIT,SS1
      PARAMETER (N_R=128,N_T=64,N_P=128)
      PARAMETER (NRTP=N_R*N_T*N_P)
      PARAMETER (JTOT=6)
      DIMENSION DD(3),REQ(3),BETA(3)
      DIMENSION DM(3),EPS(3),RM12(3),RM123(3)
      DIMENSION AA(2*JTOT+1,2*JROT+1)
      DIMENSION RHO(N_R),THE(N_T),PHI(N_P)
      DIMENSION CLGD(-JROT:JROT),CLGD1(2*JROT+1)
      DIMENSION PLNITA(30,30)
      DIMENSION PSINIT(NRTP,2*JTOT+1)
CCC
      ZI=DCMPLX(0.0D0,1.0D0)
      PI=4.0D0*DATAN(1.0D0)
      HBAR=0.06350781278D0
      EVEPS=0.9648533822D0
      DR=(RMAX-RMIN)/DFLOAT(N_R) 
      DT=PI/DFLOAT(N_T)/2.0D0
      DP=2.0D0*PI/DFLOAT(N_P)
CCC   
CC  ORBITAL ANGULAR MOMENTUM IS EQUAL TO TOTAL ANGULAR MOMENTUM
CCC
      LTOT=JTOT
CCC
CC  THE MORSE PARAMETER FOR THE REACTANT CHANNEL: ICH
CCC
      DE=DD(ICH)
      RE=REQ(ICH)
      BE=BETA(ICH)
      DICH=DM(ICH)
      EPSICH=EPS(ICH)
      AMR=RM12(ICH)
      RMCM=RM123(ICH)
CCC
CCC   PARAMETERS DEFINING GWP FOR THE TRANSLATIONAL MOTION
CCC 
      EK=EKIN*EVEPS
      AK_0=SQRT(2.0D0*RMCM*EK)/HBAR
      S=SQRT(2.0D0*AMR*DE)*2.0D0/(HBAR*BE)
      R2_0=5.5D0
      R2_F=1.7D0
      SIGMA=0.20D0
      QQ=2.0D0*(R2_0-R2_F)/AK_0
      A1=1.0D0/(4.0D0*SIGMA*SIGMA-ZI*QQ)
      A2=8.0D0*SIGMA*SIGMA/(16.0D0*SIGMA**4+QQ*QQ)
      AN=(A2/PI)**0.25D0
CCC
CC  EVALUATING THE A_k,mu MATRIX
CCC
      JR=JROT
      JT=JTOT
      ISIGN=1
      CALL AKMU(JR,JT,AA,ISIGN)
CCC
Cc  CLEBSCH – GORDAN COEFFICIENTS
CCC
      WRITE(166,*)
      WRITE(166,*)' ','JROT',' ','M',' ','LTOT',' ','0'
     1,' ','JTOT',' ','M',' ','VAL',' ','CG COFF.'
      WRITE(166,*)
      DO M=0,JROT,1
         CALL CLD(JROT,M,LTOT,0,JTOT,M,VAL)
         WRITE(166,*)JROT,M,LTOT,0,JTOT,M,VAL
         WRITE(166,*)
         IF (M.EQ.0.0D0) THEN
            CLGD(M)=VAL
         ELSE
            CLGD(M)=VAL
            CLGD(-M)=VAL*(-1.0D0)**(-JTOT+JROT+LTOT)
         ENDIF
      ENDDO
      WRITE(166,*)
      WRITE(166,*) 'CLEBSCH – GORDAN COEFFICIENTS'
      WRITE(166,*)
      DO I=-JROT,JROT,1
      WRITE(166,*)CLGD(I)
      ENDDO
CCC
CC NORMALIZATION OF THE CLEBSCH – GORDAN COEFFICIENTS
CCC
      S1=0.0D0
      II=1
      DO M=-JROT,JROT,1
         CLGD1(II)=CLGD(M)
         S1=S1+(2.0D0*LTOT+1.0D0)*CLGD1(II)*CLGD1(II)
         II=II+1
      ENDDO
      WRITE(166,*) 
      WRITE(166,*)'NORM. OF THE CLEBSCH – GORDAN COEFFICIENTS'
      WRITE(166,*) 
      WRITE(166,*)'S1=',S1
      WRITE(166,*) 
      S2=0.0D0
      DO L=1,2*JTOT+1
         DO M=1,2*JROT+1
           S2=S2+CONJG(AA(L,M))*AA(L,M)*(2*LTOT+1)*CLGD1(M)*CLGD1(M)
         ENDDO
      ENDDO
      WRITE(166,*) 
      WRITE(166,*)'NORM. OF THE PRODUCT BETWEEN' 
      WRITE(166,*)'CLEBSCH – GORDAN AND A_k,mu'
      WRITE(166,*) 
      WRITE(166,*)'S2=',S2
      WRITE(166,*) 
CC
      WRITE(166,*)
      WRITE(166,*)'ELEMENT OF THE PRODUCT BETWEEN' 
      WRITE(166,*)'CLEBSCH – GORDAN AND A_k,mu'
      WRITE(166,*)
      DO L=1,2*JTOT+1
         S3=DCMPLX(0.0D0,0.0D0)
         DO M=1,2*JROT+1
            S3=S3+AA(L,M)*CLGD1(M)
         ENDDO
         WRITE(166,*)L,DBLE(S3),DIMAG(S3)
         WRITE(166,*)
      ENDDO
      WRITE(166,*) 
      WRITE(166,*) 

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
               R1=DSQRT(0.5D0*RH*RH*(1.0D0+DSIN(TH)*
     1            DCOS(PH-EPSICH)))*DICH
               R2=DSQRT(0.5D0*RH*RH*(1.0D0-DSIN(TH)*
     1            DCOS(PH-EPSICH)))/DICH
               DN=1.0D0-DSIN(TH)*DSIN(TH)*
     1            DCOS(PH-EPSICH)*DCOS(PH-EPSICH)
               DNO=DSQRT(DN) 
               SNT=DCOS(TH)/DNO
               CNT=-DSIN(TH)*DSIN(PH-EPSICH)/DNO
               AJOC=RH*DSIN(TH)/(2.0D0*DNO)
               ARG=DEXP(-BE*(R1-RE))
               CALL SAFI(KMORSE,S,ARG,AFI)
               PSIR=DSQRT(BE)*AFI
CC   
               XI=AN*CDEXP(-ZI*AK_0*R2-A1*(R2-R2_0)**2)
CC
               IF(JTOT.NE.0) THEN
               FA=DSQRT(2.0D0*PI*SNT*AJOC*(2.0D0*LTOT+1))
     1            *PSIR*XI*(-1)**(JROT-LTOT)
               ELSE
               FA=DSQRT(2.0D0*PI*SNT*AJOC)*PSIR*XI
               ENDIF
      
               IF(JTOT.EQ.0) THEN
               CALL BLEG(CNT,PLNITA,JROT+1,0) 
               ELSE
               CALL BLEG(CNT,PLNITA,JROT+1,1) 
               ENDIF
               
CC
               DO L=1,2*JTOT+1
                  IF(JTOT.NE.0) THEN
                  SS=DCMPLX(0.0D0,0.0D0)
                  DO M=1,2*JROT+1
                     NP=ABS(M-JROT-1)
                     SS=SS+CLGD1(M)*AA(L,M)*PLNITA(JROT+1,NP+1) 
                  ENDDO
                  PSINIT(KK,L)=FA*SS
                  ELSE
                  PSINIT(KK,L)=FA*PLNITA(JROT+1,1)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      SS2=0.0D0
      WRITE(166,*)
      WRITE(166,*)'WEIGHT FACTORS BEFORE NORMALIZATION'
      WRITE(166,*)
      DO L=1,2*JTOT+1
         SS1=DCMPLX(0.0D0,0.0D0)
         DO KK=1,NRTP
            SS1=SS1+CONJG(PSINIT(KK,L))*PSINIT(KK,L)*DR*DT*DP
         ENDDO
         WRITE(166,44)L-JTOT-1,DBLE(SS1),DIMAG(SS1)
         SS2=SS2+DBLE(SS1)
      ENDDO
      WRITE(166,*)
      WRITE(166,*)'NORM.'
      WRITE(166,*)
      WRITE(166,*)SS2
      WRITE(166,*)
      WRITE(166,*)
      WRITE(166,*)
      DO KK=1,NRTP
         DO L=1,2*JTOT+1
            PSINIT(KK,L)=PSINIT(KK,L)/SQRT(SS2)
         ENDDO
      ENDDO
      WRITE(166,*)
      WRITE(166,*)'WEIGHT FACTORS AFTER NORMALIZATION'
      WRITE(166,*)

      SS2=0.0D0
      DO L=1,2*JTOT+1
         SS1=DCMPLX(0.0D0,0.0D0)
         DO KK=1,NRTP
            SS1=SS1+CONJG(PSINIT(KK,L))*PSINIT(KK,L)*DR*DT*DP
         ENDDO
         WRITE(166,44)L-JTOT-1,DBLE(SS1),DIMAG(SS1)
         SS2=SS2+DBLE(SS1)
      ENDDO
      WRITE(166,*)
      WRITE(166,*)'NORM.'
      WRITE(166,*)
      WRITE(166,*)SS2
  44  FORMAT(1X,I4,1X,4F12.6)
      RETURN
      END
CCCCCC
      SUBROUTINE PSIJAC(JTT,LP,ICH,DD,REQ,BETA,DM,EPS,RM12,
     1RHO,THE,PHI,RSTAR,AAP,VALP,PSIREIN,PSIIMIN,PRE1,PIM1)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 AAP
      COMPLEX*16 SP,PRIJ
      PARAMETER (N_R=128,N_T=64,N_P=128)
      PARAMETER (NRTP=N_R*N_T*N_P)
      PARAMETER (JTOT=6)
      PARAMETER (JTR=JTOT+1)    ! FOR EVEN JROT
C     PARAMETER (JTR=JTOT)      ! FOR ODD JROT
      PARAMETER (NV=10)
      PARAMETER (NJ=12)
      DIMENSION DD(3),REQ(3),BETA(3)
      DIMENSION DM(3),EPS(3),RM12(3)
      DIMENSION RHO(N_R),THE(N_T),PHI(N_P)
      DIMENSION AAP(0:NJ,JTR,2*NJ+1),VALP(0:NJ,0:NJ+JTOT,2*NJ+1)
      DIMENSION PSIREIN(NRTP,JTOT+1)
      DIMENSION PSIIMIN(NRTP,JTOT+1)
      DIMENSION PRHORE(N_R),PRHOIM(N_R)
      DIMENSION RO(N_R),PRIJ(JTT)
      DIMENSION PLNITA(30,30)
      DIMENSION PRE1(0:NV,0:NJ),PIM1(0:NV,0:NJ)
CC
      PI=4.0D0*DATAN(1.0D0)
      HBAR=0.06350781278D0
CC
      DT=PI/DFLOAT(N_T)/2.0D0
      DP=2.0D0*PI/DFLOAT(N_P)
CC
C  MORSE PARAMETERS AND CONSTANTS
CC
      DE=DD(ICH)
      RE=REQ(ICH)
      BE=BETA(ICH)
      DICH=DM(ICH)
      EPSICH=EPS(ICH)
      AMR=RM12(ICH)
      S=SQRT(2.0D0*AMR*DE)*2.0D0/(HBAR*BE)
CC
      JT=JTOT
      DO KM=0,NV
      DO JR=0,NJ
      PRE1(KM,JR)=0.0D0
      PIM1(KM,JR)=0.0D0
      ENDDO
      ENDDO

      DO J=1,N_T
         TH=THE(J)
         DO K=1,N_P
            PH=PHI(K)
            RH=DSQRT(2.0D0)*DICH*RSTAR/
     1         DSQRT(1.0D0-DSIN(TH)*DCOS(PH-EPSICH))
            IF(RH.GT.9.90D0) GO TO 20
CC
            R1=DICH*RH*DSQRT(1.0D0+DSIN(TH)*DCOS(PH-EPSICH))/
     1         SQRT(2.0D0)
CC  
            DNO=1.0D0-DSIN(TH)*DSIN(TH)*DCOS(PH-EPSICH)*
     1          DCOS(PH-EPSICH)
            CNT=-DSIN(TH)*DSIN(PH-EPSICH)/SQRT(DNO)
            AJOC=0.5D0*DICH*DICH*RSTAR*
     1           SQRT(DSIN(2.0D0*TH))
     1           /(1.0D0-DSIN(TH)*DCOS(PH-EPSICH))/DNO
            DO L=1,JTT
            DO IR=1,N_R
            LL=K+(J-1)*N_P+(IR-1)*N_T*N_P
            PRHORE(IR)=PSIREIN(LL,L)
            PRHOIM(IR)=PSIIMIN(LL,L)
            RO(IR)=RHO(IR)
            ENDDO
CC
CC THE INTERPOLATION HAS BEEN DONE HERE.
CC
            DRE1=ABS((PRHORE(2)-PRHORE(1))/(RO(2)-RO(1)))
            DREN=ABS((PRHORE(N_R)-PRHORE(N_R-1))/(RO(2)-RO(1)))
            CALL SPLIN(RO,PRHORE,N_R,DRE1,DREN,RH,PREVAL)
            DIM1=ABS((PRHOIM(2)-PRHOIM(1))/(RO(2)-RO(1)))
            DIMN=ABS((PRHOIM(N_R)-PRHOIM(N_R-1))/(RO(2)-RO(1)))
            CALL SPLIN(RO,PRHOIM,N_R,DIM1,DIMN,RH,PIMVAL)
            PRIJ(L)=DCMPLX(PREVAL,PIMVAL)
            ENDDO
            DO KM=0,NV
            ARG=DEXP(-BE*(R1-RE))
            CALL SAFI(KM,S,ARG,AFI)
            PSIR=DSQRT(BE)*AFI
            DO JR=0,NJ
CCC
            IF(JT.EQ.0) THEN
            CALL BLEG(CNT,PLNITA,JR+1,0)
            ELSE
            CALL BLEG(CNT,PLNITA,JR+1,1)
            ENDIF
CC
          SP=DCMPLX(0.0D0,0.0D0)
          DO L=1,JTT
          IF(JT.NE.0) THEN
          DO M=1,2*JR+1
            NP=ABS(M-JR-1)
            MP=M-JR-1
            IF(MP.GE.0) THEN
            CLGD1=DSQRT(2.0D0*PI*(2.0D0*DFLOAT(LP)+1.0D0))*
     1            (-1.0D0)**(JR-LP)*VALP(JR,LP,M)
            ELSE
            CLGD1=DSQRT(2.0D0*PI*(2.0D0*DFLOAT(LP)+1.0D0))*
     1            (-1.0D0)**(JR-LP)*VALP(JR,LP,M)*(-1.0D0)**(-JT+JR+LP)
            ENDIF
            SP=SP+CLGD1*CONJG(AAP(JR,L,M))*PLNITA(JR+1,NP+1)
     1           *PRIJ(L)
          ENDDO
          ELSE
            SP=SP+SQRT(2.0D0*PI)*PLNITA(JR+1,1)*PRIJ(L)
          ENDIF
          ENDDO
          PRE1(KM,JR)=PRE1(KM,JR)+
     1       4.0D0*DBLE(SP)*RSTAR*R1*PSIR*AJOC*DT*DP/DSQRT(RH**5)
          PIM1(KM,JR)=PIM1(KM,JR)+
     1       4.0D0*DIMAG(SP)*RSTAR*R1*PSIR*AJOC*DT*DP/DSQRT(RH**5)
          ENDDO
          ENDDO
  20    CONTINUE
      ENDDO
      ENDDO
      RETURN
      END
CCCCCC
      SUBROUTINE PSIJAC1(JTT,LP,ICH,DD,REQ,BETA,DM,EPS,RM12,
     1RHO,THE,PHI,RSTAR,AAP,VALP,ZP,PSIRE,PSIIM,
     2PRE1,PIM1,PRE2,PIM2,PRE3,PIM3)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 AAP
      COMPLEX*16 SP1,SP2,SP3,PRIJ1,PRIJ2,PRIJ3
      PARAMETER (N_R=128,N_T=64,N_P=128)
      PARAMETER (JTOT=6)
      PARAMETER (JTR=JTOT+1)    ! FOR EVEN JROT
C     PARAMETER (JTR=JTOT)      ! FOR ODD JROT
      PARAMETER (NRTP=N_R*N_T*N_P)
      PARAMETER (NJRTP=NRTP*JTR)
      PARAMETER (NSTATE=3)
      PARAMETER (NSJRTP=NSTATE*NJRTP)
      PARAMETER (NV=10)
      PARAMETER (NJ=12)
      DIMENSION DD(3),REQ(3),BETA(3)
      DIMENSION DM(3),EPS(3),RM12(3)
      DIMENSION RHO(N_R),THE(N_T),PHI(N_P)
      DIMENSION AAP(0:NJ,JTR,2*NJ+1),VALP(0:NJ,0:NJ+JTOT,2*NJ+1)
      DIMENSION PSIRE(NSJRTP),PSIIM(NSJRTP)
      DIMENSION ZP(NSTATE,NSTATE,NRTP)
      DIMENSION PSRE1(NRTP,JTR),PSIM1(NRTP,JTR)
      DIMENSION PSRE2(NRTP,JTR),PSIM2(NRTP,JTR)
      DIMENSION PSRE3(NRTP,JTR),PSIM3(NRTP,JTR)
      DIMENSION PRHORE1(N_R),PRHOIM1(N_R)
      DIMENSION PRHORE2(N_R),PRHOIM2(N_R)
      DIMENSION PRHORE3(N_R),PRHOIM3(N_R)
      DIMENSION RO(N_R)
      DIMENSION PRIJ1(JTT),PRIJ2(JTT),PRIJ3(JTT)
      DIMENSION PLNITA(30,30)
      DIMENSION PRE1(0:NV,0:NJ),PIM1(0:NV,0:NJ)
      DIMENSION PRE2(0:NV,0:NJ),PIM2(0:NV,0:NJ)
      DIMENSION PRE3(0:NV,0:NJ),PIM3(0:NV,0:NJ)
!$OMP PARALLEL
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(JTT,LP,ICH,DD,REQ,BETA,DM,EPS,RM12,
!$OMP$ RHO,THE,PHI,RSTAR,AAP,VALP,ZP,PSIRE,PSIIM,
!$OMP$ PRE1,PIM1,PRE2,PIM2,PRE3,PIM3)
      EVEPS=0.9648533822D0
      PI=4.0D0*DATAN(1.0D0)
      HBAR=0.06350781278D0
      DT=PI/DFLOAT(N_T)/2.0D0
      DP=2.0D0*PI/DFLOAT(N_P)
!
      DO KK=1,JTT
        DO I=1,NRTP
          IS1=1
          II1=KK+(IS1-1)*JTT
          JJ1=I+(II1-1)*NRTP
          IS2=2
          II2=KK+(IS2-1)*JTT
          JJ2=I+(II2-1)*NRTP
          IS3=3
          II3=KK+(IS3-1)*JTT
          JJ3=I+(II3-1)*NRTP
          PSRE1(I,KK)=ZP(1,1,I)*PSIRE(JJ1)+ZP(2,1,I)*PSIRE(JJ2)
     1               +ZP(3,1,I)*PSIRE(JJ3)
          PSIM1(I,KK)=ZP(1,1,I)*PSIIM(JJ1)+ZP(2,1,I)*PSIIM(JJ2)
     1               +ZP(3,1,I)*PSIIM(JJ3)
          PSRE2(I,KK)=ZP(1,2,I)*PSIRE(JJ1)+ZP(2,2,I)*PSIRE(JJ2)
     1               +ZP(3,2,I)*PSIRE(JJ3)
          PSIM2(I,KK)=ZP(1,2,I)*PSIIM(JJ1)+ZP(2,2,I)*PSIIM(JJ2)
     1               +ZP(3,2,I)*PSIIM(JJ3)
          PSRE3(I,KK)=ZP(1,3,I)*PSIRE(JJ1)+ZP(2,3,I)*PSIRE(JJ2)
     1               +ZP(3,3,I)*PSIRE(JJ3)
          PSIM3(I,KK)=ZP(1,3,I)*PSIIM(JJ1)+ZP(2,3,I)*PSIIM(JJ2)
     1               +ZP(3,3,I)*PSIIM(JJ3)
        ENDDO
      ENDDO
CC
C  MORSE PARAMETERS AND CONSTANTS
CC
      DE1=DD(ICH)
      DE2=2.7942D0*EVEPS
      DE3=17.1020D0*EVEPS
      RE1=REQ(ICH)
      RE2=1.0600D0
      RE3=0.1821D0
      BE1=BETA(ICH)
      BE2=1.3594D0
      BE3=1.6456D0
      DICH=DM(ICH)
      EPSICH=EPS(ICH)
      AMR=RM12(ICH)
      S1=SQRT(2.0D0*AMR*DE1)*2.0D0/(HBAR*BE1)
      S2=SQRT(2.0D0*AMR*DE2)*2.0D0/(HBAR*BE2)
      S3=SQRT(2.0D0*AMR*DE3)*2.0D0/(HBAR*BE3)
CC
      JT=JTOT
      DO KM=0,NV
        DO JR=0,NJ
          PRE1(KM,JR)=0.0D0
          PIM1(KM,JR)=0.0D0
          PRE2(KM,JR)=0.0D0
          PIM2(KM,JR)=0.0D0
          PRE3(KM,JR)=0.0D0
          PIM3(KM,JR)=0.0D0
        ENDDO
      ENDDO
CC
      DO J=1,N_T
        TH=THE(J)
        DO K=1,N_P
          PH=PHI(K)
          RH=DSQRT(2.0D0)*DICH*RSTAR
     1      /DSQRT(1.0D0-DSIN(TH)*DCOS(PH-EPSICH))
          IF(RH.GT.9.90D0) GO TO 20
            R1=DICH*RH*DSQRT(1.0D0+DSIN(TH)*DCOS(PH-EPSICH))
     1        /SQRT(2.0D0)
            DNO=1.0D0-DSIN(TH)*DSIN(TH)*DCOS(PH-EPSICH)
     1         *DCOS(PH-EPSICH)
            CNT=-DSIN(TH)*DSIN(PH-EPSICH)/SQRT(DNO)
            AJOC=0.5D0*DICH*DICH*RSTAR
     1          *SQRT(DSIN(2.0D0*TH))
     1          /(1.0D0-DSIN(TH)*DCOS(PH-EPSICH))/DNO
            DO L=1,JTT
              DO IR=1,N_R
                LL=K+(J-1)*N_P+(IR-1)*N_T*N_P
                PRHORE1(IR)=PSRE1(LL,L)
                PRHOIM1(IR)=PSIM1(LL,L)
                PRHORE2(IR)=PSRE2(LL,L)
                PRHOIM2(IR)=PSIM2(LL,L)
                PRHORE3(IR)=PSRE3(LL,L)
                PRHOIM3(IR)=PSIM3(LL,L)
                RO(IR)=RHO(IR)
              ENDDO
CC
CC THE INTERPOLATION HAS BEEN DONE HERE.
CC
              DRE11=ABS((PRHORE1(2)-PRHORE1(1))/(RO(2)-RO(1)))
              DREN1=ABS((PRHORE1(N_R)-PRHORE1(N_R-1))/(RO(2)-RO(1)))
              CALL SPLIN1(RO,PRHORE1,N_R,DRE11,DREN1,RH,PREVAL1)
              DIM11=ABS((PRHOIM1(2)-PRHOIM1(1))/(RO(2)-RO(1)))
              DIMN1=ABS((PRHOIM1(N_R)-PRHOIM1(N_R-1))/(RO(2)-RO(1)))
              CALL SPLIN1(RO,PRHOIM1,N_R,DIM11,DIMN1,RH,PIMVAL1)
              PRIJ1(L)=DCMPLX(PREVAL1,PIMVAL1)
CC
              DRE12=ABS((PRHORE2(2)-PRHORE2(1))/(RO(2)-RO(1)))
              DREN2=ABS((PRHORE2(N_R)-PRHORE2(N_R-1))/(RO(2)-RO(1)))
              CALL SPLIN1(RO,PRHORE2,N_R,DRE12,DREN2,RH,PREVAL2)
              DIM12=ABS((PRHOIM2(2)-PRHOIM2(1))/(RO(2)-RO(1)))
              DIMN2=ABS((PRHOIM2(N_R)-PRHOIM2(N_R-1))/(RO(2)-RO(1)))
              CALL SPLIN1(RO,PRHOIM2,N_R,DIM12,DIMN2,RH,PIMVAL2)
              PRIJ2(L)=DCMPLX(PREVAL2,PIMVAL2)
CC
              DRE13=ABS((PRHORE3(2)-PRHORE3(1))/(RO(2)-RO(1)))
              DREN3=ABS((PRHORE3(N_R)-PRHORE3(N_R-1))/(RO(2)-RO(1)))
              CALL SPLIN1(RO,PRHORE3,N_R,DRE13,DREN3,RH,PREVAL3)
              DIM13=ABS((PRHOIM3(2)-PRHOIM3(1))/(RO(2)-RO(1)))
              DIMN3=ABS((PRHOIM3(N_R)-PRHOIM3(N_R-1))/(RO(2)-RO(1)))
              CALL SPLIN1(RO,PRHOIM3,N_R,DIM13,DIMN3,RH,PIMVAL3)
              PRIJ3(L)=DCMPLX(PREVAL3,PIMVAL3)
            ENDDO
CC
            DO KM=0,NV
              ARG1=DEXP(-BE1*(R1-RE1))
              ARG2=DEXP(-BE2*(R1-RE2))
              ARG3=DEXP(-BE3*(R1-RE3))
              CALL SAFI1(KM,S1,ARG1,AFI1)
              CALL SAFI1(KM,S2,ARG2,AFI2)
              CALL SAFI1(KM,S3,ARG3,AFI3)
              PSIR1=DSQRT(BE1)*AFI1
              PSIR2=DSQRT(BE2)*AFI2
              PSIR3=DSQRT(BE3)*AFI3
            DO JR=0,NJ
CCC
            IF(JT.EQ.0) THEN
            CALL BLEG1(CNT,PLNITA,JR+1,0)
            ELSE
            CALL BLEG1(CNT,PLNITA,JR+1,1)
            ENDIF
CC
          SP1=DCMPLX(0.0D0,0.0D0)
          SP2=DCMPLX(0.0D0,0.0D0)
          SP3=DCMPLX(0.0D0,0.0D0)
          DO L=1,JTT
          IF(JT.NE.0) THEN
          DO M=1,2*JR+1
            NP=ABS(M-JR-1)
            MP=M-JR-1
            IF(MP.GE.0) THEN
            CLGD1=DSQRT(2.0D0*PI*(2.0D0*DFLOAT(LP)+1.0D0))*
     1            (-1.0D0)**(JR-LP)*VALP(JR,LP,M)
            ELSE
            CLGD1=DSQRT(2.0D0*PI*(2.0D0*DFLOAT(LP)+1.0D0))*
     1            (-1.0D0)**(JR-LP)*VALP(JR,LP,M)*(-1.0D0)**(-JT+JR+LP)
            ENDIF
            SP1=SP1+CLGD1*CONJG(AAP(JR,L,M))*PLNITA(JR+1,NP+1)
     1           *PRIJ1(L)
            SP2=SP2+CLGD1*CONJG(AAP(JR,L,M))*PLNITA(JR+1,NP+1)
     1           *PRIJ2(L)
            SP3=SP3+CLGD1*CONJG(AAP(JR,L,M))*PLNITA(JR+1,NP+1)
     1           *PRIJ3(L)
          ENDDO
          ELSE
            SP1=SP1+SQRT(2.0D0*PI)*PLNITA(JR+1,1)*PRIJ1(L)
            SP2=SP2+SQRT(2.0D0*PI)*PLNITA(JR+1,1)*PRIJ2(L)
            SP3=SP3+SQRT(2.0D0*PI)*PLNITA(JR+1,1)*PRIJ3(L)
          ENDIF
          ENDDO
          PRE1(KM,JR)=PRE1(KM,JR)+
     1       4.0D0*DBLE(SP1)*RSTAR*R1*PSIR1*AJOC*DT*DP/DSQRT(RH**5)
          PIM1(KM,JR)=PIM1(KM,JR)+
     1       4.0D0*DIMAG(SP1)*RSTAR*R1*PSIR1*AJOC*DT*DP/DSQRT(RH**5)
CC
          PRE2(KM,JR)=PRE2(KM,JR)+
     1       4.0D0*DBLE(SP2)*RSTAR*R1*PSIR2*AJOC*DT*DP/DSQRT(RH**5)
          PIM2(KM,JR)=PIM2(KM,JR)+
     1       4.0D0*DIMAG(SP2)*RSTAR*R1*PSIR2*AJOC*DT*DP/DSQRT(RH**5)
CC
          PRE3(KM,JR)=PRE3(KM,JR)+
     1       4.0D0*DBLE(SP3)*RSTAR*R1*PSIR3*AJOC*DT*DP/DSQRT(RH**5)
          PIM3(KM,JR)=PIM3(KM,JR)+
     1       4.0D0*DIMAG(SP3)*RSTAR*R1*PSIR3*AJOC*DT*DP/DSQRT(RH**5)
          ENDDO
          ENDDO
  20    CONTINUE
      ENDDO
      ENDDO
!$OMP END PARALLEL
      RETURN
      END
C
C HANKEL FUNCTION
C
      SUBROUTINE HANKEL(NL,X,VAL)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,ANUM,VAL
      ZI=DCMPLX(0.0D0,1.0D0)
      VAL=DCMPLX(0.0D0,0.0D0)
      DO L=0,NL
      NLP=NL+L
      NLN=NL-L
      ANUM=FACT(NLP)*ZI**L
      DENO=FACT(L)*FACT(NLN)*(2.0D0*X)**L
      VAL=VAL+ANUM/DENO
      ENDDO
      VAL=VAL*EXP(ZI*X)*(-ZI)**NL
      RETURN
      END
CCCCCC
      SUBROUTINE VIBMATRC(RM12,VALRC11,VALRC22)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NRR=100)
      PARAMETER (NV=10)
      PARAMETER (NCH=3)
      DIMENSION RM12(NCH)
      DIMENSION VALRC11(NCH,0:NV,0:NV)
      DIMENSION VALRC22(NCH,0:NV,0:NV)
      HBAR=0.06350781278D0
      EVEPS=0.9648533822D0
      RMIN=0.02D0
      RMAX=1.34D0
      DR=(RMAX-RMIN)/DFLOAT(NRR)
      DE1=4.580841D0
      DE2=2.7942D0*EVEPS
      RE1=0.741020D0
      RE2=1.0600D0
      BE1=1.950477D0
      BE2=1.3594D0
      DO IFCH=1,NCH
        AMR=RM12(IFCH)
        S1=SQRT(2.0D0*AMR*DE1)*2.0D0/(HBAR*BE1)
        S2=SQRT(2.0D0*AMR*DE2)*2.0D0/(HBAR*BE2)
        DO KM=0,NV
          DO KM1=0,NV
            VALRC11(IFCH,KM,KM1)=0.0D0
            VALRC22(IFCH,KM,KM1)=0.0D0
            DO I=1,NRR
              R1=RMIN+(I-1)*DR
              ARG1=DEXP(-BE1*(R1-RE1))
              ARG2=DEXP(-BE2*(R1-RE2))
              CALL SAFI(KM,S1,ARG1,AFI1KM)
              CALL SAFI(KM,S2,ARG2,AFI2KM)
              CALL SAFI(KM1,S1,ARG1,AFI1KM1)
              CALL SAFI(KM1,S2,ARG2,AFI2KM1)
              PSIR1KM=DSQRT(BE1)*AFI1KM
              PSIR2KM=DSQRT(BE2)*AFI2KM
              PSIR1KM1=DSQRT(BE1)*AFI1KM1
              PSIR2KM1=DSQRT(BE2)*AFI2KM1
              VALRC11(IFCH,KM,KM1)=VALRC11(IFCH,KM,KM1)
     1                            +PSIR1KM*PSIR1KM1*DR
              VALRC22(IFCH,KM,KM1)=VALRC22(IFCH,KM,KM1)
     1                            +PSIR2KM*PSIR2KM1*DR
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
CCCCCC
      SUBROUTINE VIBMATINF(RM12,VALINF12,VALINF21)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NRR=100)
      PARAMETER (NV=10)
      PARAMETER (NCH=3)
      DIMENSION RM12(NCH)
      DIMENSION VALINF12(NCH,0:NV,0:NV)
      DIMENSION VALINF21(NCH,0:NV,0:NV)
      HBAR=0.06350781278D0
      EVEPS=0.9648533822D0
      RMIN=1.34D0
      RMAX=8.0D0
      DR=(RMAX-RMIN)/DFLOAT(NRR)
      DE1=4.580841D0
      DE2=2.7942D0*EVEPS
      RE1=0.741020D0
      RE2=1.0600D0
      BE1=1.950477D0
      BE2=1.3594D0
      DO IFCH=1,NCH
        AMR=RM12(IFCH)
        S1=SQRT(2.0D0*AMR*DE1)*2.0D0/(HBAR*BE1)
        S2=SQRT(2.0D0*AMR*DE2)*2.0D0/(HBAR*BE2)
        DO KM=0,NV
          DO KM1=0,NV
            VALINF12(IFCH,KM,KM1)=0.0D0
            VALINF21(IFCH,KM,KM1)=0.0D0
            DO I=1,NRR
              R1=RMIN+(I-1)*DR
              ARG1=DEXP(-BE1*(R1-RE1))
              ARG2=DEXP(-BE2*(R1-RE2))
              CALL SAFI(KM,S1,ARG1,AFI1KM)
              CALL SAFI(KM,S2,ARG2,AFI2KM)
              CALL SAFI(KM1,S1,ARG1,AFI1KM1)
              CALL SAFI(KM1,S2,ARG2,AFI2KM1)
              PSIR1KM=DSQRT(BE1)*AFI1KM
              PSIR2KM=DSQRT(BE2)*AFI2KM
              PSIR1KM1=DSQRT(BE1)*AFI1KM1
              PSIR2KM1=DSQRT(BE2)*AFI2KM1
              VALINF12(IFCH,KM,KM1)=VALINF12(IFCH,KM,KM1)
     1                             +PSIR1KM*PSIR2KM1*DR
              VALINF21(IFCH,KM,KM1)=VALINF21(IFCH,KM,KM1)
     1                             +PSIR2KM*PSIR1KM1*DR
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!****************************
!                           *
!     MORSE WAVEFUNCTION    *
!                           * 
!****************************
      SUBROUTINE SAFI(N,S,Y,AFI)
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
      CALL SGAMMLN(XX,GAMMLN)
      DLNF=GAMMLN
      FN=-0.5D0*DLNF
    1 XX=S-RN
      CALL SGAMMLN(XX,GAMMLN)
      DLNG=GAMMLN
      XX=XX-RN
      CALL SGAMMLN(XX,GAMMLN)
      DLNG1=GAMMLN
      XX=XX-1.0D0
      CALL SGAMMLN(XX,GAMMLN)
      DLNG2=GAMMLN
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
            F=F*(N1-K)/K/(1.0D0-(2.0D0*RN+1.0D0-DFLOAT(K))/S)
         ENDDO
         SU=SU+F*SIG*Y**(1.0D0*I1)
      ENDDO
    4 AFI=FA*SU
      GO TO 7
    5 AFI=0.0D0
    7 GOTO 21
   10 ALF=S-2.0D0*RN-1.0D0
      RL(1)=1.0D0
      Z=S*Y
      RL(2)=ALF+1.0D0-Z
      NS=N+1
      DO I=3,NS
         RN1=DFLOAT(I-1)
         RL(I)=((2.D0*RN1-1.D0+ALF-Z)*RL(I-1)
     1        -(RN1-1.D0+ALF)*RL(I-2))/RN1
      ENDDO
      XX=1.0D0*RN+ALF+1.0D0
      CALL SGAMMLN(XX,GAMMLN)
      DLNG=GAMMLN
      ALF1=ALF+1.0D0
      CALL SGAMMLN(ALF1,GAMMLN)
      DLNG1=GAMMLN
      RN1=DFLOAT(N+1)
      CALL SGAMMLN(RN1,GAMMLN)
      DLNG2=GAMMLN
      FARG=DLNG2+DLNG1-DLNG
      IF (RL(NS).LT.2.0D-9) GO TO 12
      DLNG=DLOG(RL(NS))
      AFI=DEXP(ARG+FARG+DLNG)
      GOTO 21
   12 AFI=DEXP(ARG+FARG)*RL(NS)
      GOTO 21
   14 AN=1.0D0
      Z=S*Y
      ALF=S-2.0D0*RN-1.0D0
      DO M=N,1,-1
         RM=DFLOAT(M)
         AN=1.0D0-Z*AN*DFLOAT(N-M+1)/RM/(ALF+RM)
      ENDDO
      IF (ARG.GT.10.0D0) WRITE(166,100) Y,ARG,AN
  100 FORMAT(1X,5E10.3)
      AFI=AN*DEXP(ARG)
  21  CONTINUE
      RETURN
      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE SGAMMLN(XX,GAMMLN)
      DOUBLE PRECISION XX,ZZ,TERM,RZ2,GAMMLN
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
      WRITE(166,*)'ARGUMENT LNGAM > 1.D30'
      GAMMLN=1.0D+30
  10  CONTINUE
      RETURN
      END

C----------------------------------------------------------------------  
      SUBROUTINE SAFI1(N,S,Y,AFI)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION RL(50)
!$OMP PARALLEL
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(N,S,Y,AFI)
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
      CALL SGAMMLN1(XX,GAMMLN)
      DLNF=GAMMLN
      FN=-0.5D0*DLNF
    1 XX=S-RN
      CALL SGAMMLN1(XX,GAMMLN)
      DLNG=GAMMLN
      XX=XX-RN
      CALL SGAMMLN1(XX,GAMMLN)
      DLNG1=GAMMLN
      XX=XX-1.0D0
      CALL SGAMMLN1(XX,GAMMLN)
      DLNG2=GAMMLN
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
            F=F*(N1-K)/K/(1.0D0-(2.0D0*RN+1.0D0-DFLOAT(K))/S)
         ENDDO
         SU=SU+F*SIG*Y**(1.0D0*I1)
      ENDDO
    4 AFI=FA*SU
      GO TO 7
    5 AFI=0.0D0
    7 GOTO 21
   10 ALF=S-2.0D0*RN-1.0D0
      RL(1)=1.0D0
      Z=S*Y
      RL(2)=ALF+1.0D0-Z
      NS=N+1
      DO I=3,NS
         RN1=DFLOAT(I-1)
         RL(I)=((2.D0*RN1-1.D0+ALF-Z)*RL(I-1)
     1        -(RN1-1.D0+ALF)*RL(I-2))/RN1
      ENDDO
      XX=1.0D0*RN+ALF+1.0D0
      CALL SGAMMLN1(XX,GAMMLN)
      DLNG=GAMMLN
      ALF1=ALF+1.0D0
      CALL SGAMMLN1(ALF1,GAMMLN)
      DLNG1=GAMMLN
      RN1=DFLOAT(N+1)
      CALL SGAMMLN1(RN1,GAMMLN)
      DLNG2=GAMMLN
      FARG=DLNG2+DLNG1-DLNG
      IF (RL(NS).LT.2.0D-9) GO TO 12
      DLNG=DLOG(RL(NS))
      AFI=DEXP(ARG+FARG+DLNG)
      GOTO 21
   12 AFI=DEXP(ARG+FARG)*RL(NS)
      GOTO 21
   14 AN=1.0D0
      Z=S*Y
      ALF=S-2.0D0*RN-1.0D0
      DO M=N,1,-1
         RM=DFLOAT(M)
         AN=1.0D0-Z*AN*DFLOAT(N-M+1)/RM/(ALF+RM)
      ENDDO
      IF (ARG.GT.10.0D0) WRITE(166,100) Y,ARG,AN
  100 FORMAT(1X,5E10.3)
      AFI=AN*DEXP(ARG)
  21  CONTINUE
!$OMP END PARALLEL
      RETURN
      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE SGAMMLN1(XX,GAMMLN)
      DOUBLE PRECISION XX,ZZ,TERM,RZ2,GAMMLN
!$OMP PARALLEL
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(XX,GAMMLN)
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
      WRITE(166,*)'ARGUMENT LNGAM > 1.D30'
      GAMMLN=1.0D+30
  10  CONTINUE
!$OMP END PARALLEL
      RETURN
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
      JJ1=INT(AJ+AJ_1-AJ_2)
      JJ2=INT(AJ-AJ_1+AJ_2)
      JJ3=INT(-AJ+AJ_1+AJ_2)
      JJ4=INT(AJ_1+AJ_2+AJ+1)
      FAC1=DSQRT((2.0D0*AJ+1.0D0)*FACT(JJ1)*FACT(JJ2)
     1           *FACT(JJ3)/FACT(JJ4))
      JM1=INT(AJ+AM)
      JM2=INT(AJ-AM)
      JM3=INT(AJ_1-AM_1)
      JM4=INT(AJ_1+AM_1)
      JM5=INT(AJ_2-AM_2)
      JM6=INT(AJ_2+AM_2)
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
      CLGD=CLGD*(-1.0D0)**(J_2-J_1-M)/DSQRT(2.0D0*AJ+1.0D0)
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
      SUBROUTINE AKMU(JROT,JTOT,AA,ISIGN)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,AA,VAL,FA,BB
      PARAMETER (NT=64,NP=128)
      DIMENSION AA(2*JTOT+1,2*JROT+1)
      DIMENSION PLT(30,30)
      DIMENSION PLP(30,30)
      DIMENSION BB(2*JROT+1,2*JROT+1)
      PI=4.0D0*DATAN(1.0D0)
      DT=PI/DFLOAT(NT)
      DP=2.0D0*PI/DFLOAT(NP)
      ZI=DCMPLX(0.0D0,1.0D0)
      K1=1
      DO KK=-JTOT,JTOT,1
         M1=1
         DO MM=-JROT,JROT,1
            AA(K1,M1)=DCMPLX(0.0D0,0.0D0)
            MPJ=ABS(KK)
            NPJ=ABS(MM)
            DO II=1,NT
               T=(II-1)*DT+0.5D0*DT
               X=DCOS(T)
               CALL BLEG(X,PLT,JTOT+1,1) 
               DO JJ=1,NP
                  P=(JJ-1)*DP+0.5D0*DP
                  Y=-DSIN(T)*DCOS(P)
                  CALL BLEG(Y,PLP,JTOT+1,1) 
                  DENO=SQRT(1.0D0-DSIN(T)*DSIN(T)*DCOS(P)*DCOS(P))
                  CSXI=DSIN(T)*DSIN(P)/DENO
                  SIXI=DCOS(T)/DENO
                     CALL DEMV(CSXI,SIXI,MM,VAL)
                     FA=(-1)**KK*EXP(-ZI*KK*P)*VAL
     1                  *PLT(JTOT+1,MPJ+1)*PLP(JTOT+1,NPJ+1)
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
      IF(ISIGN.EQ.1) THEN
      WRITE(166,*)'RE(A_K,mu)'
      WRITE(166,*)
      WRITE(166,44)((DBLE(AA(II,JJ)),JJ=1,2*JROT+1),II=1,2*JTOT+1)
      WRITE(166,*)
      WRITE(166,*)
      WRITE(166,*)'Im(A_K,mu)'
      WRITE(166,*)
      WRITE(166,44)((DIMAG(AA(II,JJ)),JJ=1,2*JROT+1),II=1,2*JTOT+1)
      WRITE(166,*)
      WRITE(166,*)
      WRITE(166,*)
      WRITE(166,*)'RE(A_K,mu*A_k,mu)'
      WRITE(166,*)
      WRITE(166,44)((DBLE(BB(II,JJ)),JJ=1,2*JROT+1),II=1,2*JROT+1)
      WRITE(166,*)
      WRITE(166,*)
      WRITE(166,*)'Im(A_K,mu*A_k,mu)'
      WRITE(166,*)
      WRITE(166,44)((DIMAG(BB(II,JJ)),JJ=1,2*JROT+1),II=1,2*JROT+1)
      ENDIF
  44  FORMAT(1X,F12.6) 
      RETURN
      END
C----------------------------------------------------------------------  
C----------------------------------------------------------------------  
      SUBROUTINE AKU(KK,MM,AA)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,AA,VAL,FA
      PARAMETER (NT=64,NP=128)
      PARAMETER (JTOT=6)
      DIMENSION PLT(30,30)
      DIMENSION PLP(30,30)
      PI=4.0D0*DATAN(1.0D0)
      DT=PI/DFLOAT(NT)
      DP=2.0D0*PI/DFLOAT(NP)
      ZI=DCMPLX(0.0D0,1.0D0)
            AA=DCMPLX(0.0D0,0.0D0)
            MPJ=ABS(KK)
            NPJ=ABS(MM)
            DO II=1,NT
               T=(II-1)*DT+0.5D0*DT
               X=DCOS(T)
               CALL BLEG(X,PLT,JTOT+1,1) 
               DO JJ=1,NP
                  P=(JJ-1)*DP+0.5D0*DP
                  Y=-DSIN(T)*DCOS(P)
                  CALL BLEG(Y,PLP,JTOT+1,1) 
                  DENO=SQRT(1.0D0-DSIN(T)*DSIN(T)*DCOS(P)*DCOS(P))
                  CSXI=DSIN(T)*DSIN(P)/DENO
                  SIXI=DCOS(T)/DENO
                     CALL DEMV(CSXI,SIXI,MM,VAL)
                     FA=(-1)**KK*EXP(-ZI*KK*P)*VAL
     1                  *PLT(JTOT+1,MPJ+1)*PLP(JTOT+1,NPJ+1)
                  AA=AA+FA*DSIN(T)*DT*DP
               ENDDO
            ENDDO
        RETURN
        END 
CCC
CCC   DE MOIVRE'S THEOREM 
CCC
      SUBROUTINE DEMV(CS,SN,N,VAL)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ZI,VAL
      ZI=DCMPLX(0.0D0,1.0D0)
CCC
      M=ABS(N)
      IF(N.EQ.0) THEN
      VAL=DCMPLX(1.0D0,0.0D0)
      ENDIF
CCC
      IF(N.LT.0) THEN
      VAL=DCMPLX(1.0D0,0.0D0)
      DO I=1,M 
      VAL=VAL*(CS-ZI*SN)
      ENDDO
      ENDIF
CCC
      IF(N.GT.0) THEN
      VAL=DCMPLX(1.0D0,0.0D0)
      DO I=1,M 
      VAL=VAL*(CS+ZI*SN)
      ENDDO
      ENDIF
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
      SX=DSQRT(1.0D0-X**2)
      P(2,1)=X 
      N1=N+1
      DO I=3,N1
C         P(I,1)=((2*I-3)*P(I-1,1)*X-(I-2)*P(I-2,1))/(I-1)
         P(I,1)=2*X*P(I-1,1)-P(I-2,1)-(X*P(I-1,1)-P(I-2,1))/(I-1) 
      ENDDO
      IF (IOP.EQ.2) GOTO 21
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
      IF(IOP.EQ.3) GOTO 21 
      DO I=1,N1
         DO J=1,I
            L=I-1
            M=J-1
            P(I,J)=0.3989422804D0*P(I,J)*DSQRT((2.0D0*L+1.0D0)*
     1             FAC(L-M)/2.0D0/FAC(L+M))
         ENDDO
      ENDDO 
      GOTO 21
    4 DO I=1,N1
         P(I,1)=0.2820947918D0*DSQRT(2.0D0*DFLOAT(I)-1.0D0)*P(I,1)
      ENDDO
 21   CONTINUE
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE BLEG1(X,P,N,IOP) 
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(30,30),ARR(10)
!$OMP PARALLEL
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(X,P,N,IOP)
      ARR(1)=1.0D0
      ARR(2)=2.0D0
      ARR(3)=6.0D0
      ARR(4)=24.0D0
      ARR(5)=120.0D0
      ARR(6)=720.0D0
      ARR(7)=5040.0D0
      ARR(8)=40320.0D0
      ARR(9)=362880.0D0
      ARR(10)=3628800.0D0
      DO I=1,N
         DO J=1,N
            P(I,J)=0.0D0
         ENDDO 
      ENDDO 
      P(1,1)=1.0D0
      SX=DSQRT(1.0D0-X**2)
      P(2,1)=X 
      N1=N+1
      DO I=3,N1
C         P(I,1)=((2*I-3)*P(I-1,1)*X-(I-2)*P(I-2,1))/(I-1)
         P(I,1)=2*X*P(I-1,1)-P(I-2,1)-(X*P(I-1,1)-P(I-2,1))/(I-1) 
      ENDDO
      IF (IOP.EQ.2) GOTO 21
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
      IF(IOP.EQ.3) GOTO 21 
      DO I=1,N1
        DO J=1,I
          L=I-1
          M=J-1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(L-M.LT.1) GO TO 11
      IF(L-M.GT.10) GO TO 41
      VALLMM=ARR(L-M)
      GO TO 51
   41 SU=1.0D0
      DO I1=1,L-M
        SU=SU*DFLOAT(I1)
      ENDDO
      GO TO 31
   11 SU=1.0D0
   31 VALLMM=SU
   51 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(L+M.LT.1) GO TO 12
      IF(L+M.GT.10) GO TO 42
      VALLMP=ARR(L+M)
      GO TO 52
   42 SU=1.0D0
      DO I1=1,L+M
        SU=SU*DFLOAT(I1)
      ENDDO
      GO TO 32
   12 SU=1.0D0
   32 VALLMP=SU
   52 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          P(I,J)=0.3989422804D0*P(I,J)*DSQRT((2.0D0*L+1.0D0)*
     1             VALLMM/2.0D0/VALLMP)
        ENDDO
      ENDDO 
      GOTO 21
    4 DO I=1,N1
         P(I,1)=0.2820947918D0*DSQRT(2.0D0*DFLOAT(I)-1.0D0)*P(I,1)
      ENDDO
 21   CONTINUE
!$OMP END PARALLEL
      RETURN
      END
CCCCCC
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
         SU=SU*DFLOAT(I)
      ENDDO
      GO TO 3
    1 SU=1.0D0
    3 FAC=SU
    5 RETURN
      END
*********************************************************************
*                                                                   *
*     MORSE PARAMETERS (REQ, BETA, DD) AND ATOMIC MASSES (EM)       *
*                                                                   *
*********************************************************************
      SUBROUTINE MOSC(REQ,BETA,DD,EM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DD(3),REQ(3),BETA(3),EM(3)
      REQ(1)=0.741020D0
      REQ(2)=0.741020D0
      REQ(3)=0.741020D0
      BETA(1)=1.950477D0
      BETA(2)=1.950477D0
      BETA(3)=1.950477D0
      DD(1)=4.580841D0
      DD(2)=4.580841D0
      DD(3)=4.580841D0
      EM(1)=2.016490D0
      EM(2)=1.007825D0
      EM(3)=1.007825D0
      RETURN
      END
****************************************************************************
*                                                                          *
*    TOTAL TRIATOMIC MASS (XM) AND REDUCE MASS OF THE TRIATOM (EMU).       *
*    DIATOMIC MASSES (RM12) AND CENTER OF MASS OF THE TRIATOM (RM123)      *
*    FOR THREE CHANNEL.                                                    *
*    CHANNEL DEPENDENT CONSTANTS DM AND EPS:                               *
*    DMs ARE USED AS CONSTANTS TO DEFINE r AND R IN TERMS OF               *
*    HYPERSPHERICAL COORDINATE (RHO, THETA AND PHI).                       *
*    EPSs ARE INVOLVED TO EXPRESS THE LENGTH OF THE ARMS OF THE TRIANGLE   * 
*    FORMED BY THE THREE ATOMS.                                            *
*                                                                          * 
****************************************************************************
      SUBROUTINE HYPSET(DD,REQ,BETA,EM,EMU,DM,EPS,RM12,RM123)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DD(3),REQ(3),BETA(3),EM(3)
      DIMENSION DM(3),EPS(3),RM12(3),RM123(3)
      XM=EM(1)+EM(2)+EM(3)
      EMU=EM(1)*EM(2)*EM(3)/XM
      EMU=DSQRT(EMU)
      WRITE(166,*)
      WRITE(166,*)'CHANNEL DEPENDENT CONSTANTS DM AND EPS'
      WRITE(166,*)'TO DEFINE r AND R IN HYPERSPHERICAL COORDINATES'
      WRITE(166,*)
      DO I=1,3
        DM(I)=DSQRT(EM(I)*(1.D0-EM(I)/XM)/EMU)
      ENDDO
      EPS(1)=0.0D0
      EPS(2)=2.0D0*DATAN(EM(3)/EMU)
      EPS(3)=-2.0D0*DATAN(EM(2)/EMU)
      DO I=1,3
        WRITE(166,*)DM(I),EPS(I)
      ENDDO
      WRITE(166,*)
      RM12(1)=EM(2)*EM(3)/(EM(2)+EM(3))
      RM12(2)=EM(1)*EM(3)/(EM(1)+EM(3))
      RM12(3)=EM(1)*EM(2)/(EM(1)+EM(2))
      RM123(1)=EM(1)*(EM(2)+EM(3))/XM
      RM123(2)=EM(2)*(EM(1)+EM(3))/XM
      RM123(3)=EM(3)*(EM(1)+EM(2))/XM
      WRITE(166,*)
      WRITE(166,*)'CHANNEL DEPENDENT DIATOMIC REDUCE MASS'
      WRITE(166,*)'AND TRIATOMIC CENTER OF MASS'
      WRITE(166,*)
      DO I=1,3
        WRITE(166,*)RM12(I),RM123(I)
      ENDDO
      WRITE(166,*)
      WRITE(166,*)'TRIATOMIC REDUCE MASS'
      WRITE(166,*)
      WRITE(166,*)EMU
      WRITE(166,*)
      RETURN
      END
CCCCCC
      SUBROUTINE TQLI(D,E,N,Z)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (LAD=30)
      DIMENSION D(LAD),E(LAD),Z(LAD,LAD)
      IF(N.GT.1) THEN
        DO I=2,N
          E(I-1)=E(I)
        ENDDO
        E(N)=0.D0
        DO L=1,N
          ITER=0
    1     DO M=L,N-1
            DD=DABS(D(M))+DABS(D(M+1))
            IF (DABS(E(M))+DD.EQ.DD) GO TO 2
          ENDDO
          M=N
    2     IF(M.NE.L)THEN
            IF(ITER.EQ.30)THEN
              WRITE(166,*)'TOO MANY ITERATIONS'
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
              IF(DABS(F).GE.DABS(G))THEN
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
CCCCCC
      SUBROUTINE VPOTIM(RHO,VIM)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N_R=128,N_T=64,N_P=128)
      DIMENSION RHO(N_R),VIM(N_R)
      EVEPS=0.9648533822D0
      VMAX=0.163D0*EVEPS
      RHO1=8.0D0
      DO I=1,N_R
        IF(RHO(I).LT.RHO1)THEN
          VIM(I)=1.0D0
        ENDIF
        IF(RHO(I).GE.RHO1)THEN
          ARG=(RHO(I)-RHO1)*VMAX
          VIM(I)=EXP(-ARG)
        ENDIF
      ENDDO
      RETURN
      END
CCCCCC
      SUBROUTINE DEFAK(DR,DT,DP,AKR,AKT,AKP)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N_R=128,N_T=64,N_P=128)
      PARAMETER (NT2=2*N_T)
      DIMENSION AKR(N_R),AKT(NT2),AKP(N_P)
      PI=2.0D0*ASIN(1.0D0)
      PI2=PI*2.0D0
      ACONSR=PI2/(N_R*DR)
      ACONST=PI2/(NT2*DT)
      ACONSP=PI2/(N_P*DP)
      NHALFR=N_R/2+1
      NHALFT=NT2/2+1
      NHALFP=N_P/2+1
      DO I=1,NHALFR
        AKR(I)=ACONSR*(I-1)
      ENDDO
      DO I=NHALFR+1,N_R
        AKR(I)=ACONSR*(I-1-N_R)
      ENDDO
      DO I=1,NHALFT
        AKT(I)=ACONST*(I-1)
      ENDDO
      DO I=NHALFT+1,NT2
        AKT(I)=ACONST*(I-1-NT2)
      ENDDO
      DO I=1,NHALFP
        AKP(I)=ACONSP*(I-1)
      ENDDO
      DO I=NHALFP+1,N_P
        AKP(I)=ACONSP*(I-1-N_P)
      ENDDO
      RETURN
      END
****************************************************
*                                                  *
*  HAMILTONIAN OPERATES ON THE WAVE FUNCTION       *
*  BY THE FFT-METHOD TO CALCULATE ENERGY           *
*                                                  *
****************************************************
      SUBROUTINE ENERGY(AK,EK,EKP,EKM,PSIRE,PSIIM,
     1AKR,AKT,AKP,WIN1,WIN2,RHOD,SINN,COSN,SIND,COSD,SIN2D,
     2V11,V12,V13,V21,V22,V23,V31,V32,V33,
     3EER1,EER2,EER3,EER4,EER5,EER6,EER7,EER8,EER9,EER10,EER,
     4RMST,SNG,VOLEM,EVEPS,
     5PLAN1,PLAN2,PLAN3,PLAN4,PLAN5)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 ARR1,ARR,ART,ARP,ARP1
      PARAMETER (N_R=128,N_T=64,N_P=128)
      PARAMETER (NRTP=N_R*N_T*N_P)
      PARAMETER (NT2=2*N_T)
      PARAMETER (JTOT=6)
      PARAMETER (JTR=JTOT+1)
C      PARAMETER (JTR=JTOT)
      PARAMETER (NJRTP=NRTP*JTR)
      PARAMETER (NSTATE=3)
      PARAMETER (NSJRTP=NSTATE*NJRTP)
      DIMENSION AK(JTR),EK(JTR),EKP(JTR),EKM(JTR)
      DIMENSION PSIRE(NSJRTP),PSIIM(NSJRTP)
      DIMENSION AKR(N_R),AKT(NT2),AKP(N_P)
      DIMENSION WIN1(N_T),WIN2(N_T)
      DIMENSION RHOD(N_R),SINN(N_T),COSN(N_T)
      DIMENSION SIND(N_T),COSD(N_T),SIN2D(N_T)
      DIMENSION V11(NRTP),V12(NRTP),V13(NRTP)
      DIMENSION V21(NRTP),V22(NRTP),V23(NRTP)
      DIMENSION V31(NRTP),V32(NRTP),V33(NRTP)
      DIMENSION PRE(NRTP,JTR,NSTATE),PIM(NRTP,JTR,NSTATE)
      DIMENSION X1(NRTP),Y1(NRTP)
      DIMENSION X2(NRTP),Y2(NRTP)
      DIMENSION X3(NRTP),Y3(NRTP)
      DIMENSION EER1(NSTATE*JTR),EER2(NSTATE*JTR),EER3(NSTATE*JTR)
      DIMENSION EER4(NSTATE*JTR),EER5(NSTATE*JTR),EER6(NSTATE*JTR)
      DIMENSION EER7(NSTATE*JTR),EER8(NSTATE*JTR),EER9(NSTATE*JTR)
      DIMENSION EER10(NSTATE*JTR),EER(NSTATE*JTR)
      DIMENSION ARR1(N_R,NT2,N_P)
      DIMENSION ARR(N_R,NT2,N_P)
      DIMENSION ART(N_R,NT2,N_P)
      DIMENSION ARP(N_R,NT2,N_P)
      DIMENSION ARP1(N_R,NT2,N_P)
      INTEGER*8 PLAN1,PLAN2,PLAN3,PLAN4,PLAN5
****************************************************
*                                                  *
*   This program is distributed with permission    *
*                                                  *
****************************************************
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
!
      NC1=JTR
      DO IS=1,NSTATE
        DO JT=1,JTR
          II=JT+(IS-1)*JTR
          DO JJ=1,NRTP
            LL=JJ+(II-1)*NRTP
            PRE(JJ,JT,IS)=PSIRE(LL)
            PIM(JJ,JT,IS)=PSIIM(LL)
          ENDDO
        ENDDO
      ENDDO
!
      ICHUNK=1
!$OMP PARALLEL 
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(ICHUNK,NC1,AK,EK,EKP,EKM,PRE,PIM,PSIRE,PSIIM,
!$OMP$ AKR,AKT,AKP,WIN1,WIN2,RHOD,SINN,COSN,SIND,COSD,SIN2D,
!$OMP$ V11,V12,V13,V21,V22,V23,V31,V32,V33,
!$OMP$ EER1,EER2,EER3,EER4,EER5,EER6,EER7,EER8,EER9,EER10,EER,
!$OMP$ RMST,SNG,VOLEM,EVEPS,
!$OMP$ PLAN1,PLAN2,PLAN3,PLAN4,PLAN5)
!$OMP DO SCHEDULE(DYNAMIC,ICHUNK)
      DO III=1,NSTATE*JTR
        IIT=III
        IS=(IIT-1)/DFLOAT(NC1)
        IIT=IIT-IS*NC1
        IS=IS+1
        JT=IIT
        II=JT+(IS-1)*JTR
        IIM=(JT-1)+(IS-1)*JTR
        IIP=(JT+1)+(IS-1)*JTR
        AK_J=AK(JT)
        EJK=EK(JT)
        EJKP=EKP(JT)
        EJKM=EKM(JT)
        IF(JT-1.LT.1)THEN
          S1=0.0D0
          DO JJ=1,NRTP
            X1(JJ)=0.0D0
            Y1(JJ)=0.0D0
            S1=S1+X1(JJ)**2+Y1(JJ)**2
          ENDDO
          S1=S1*VOLEM
        ELSE
          S1=0.0D0
          DO JJ=1,NRTP
            LL=JJ+(IIM-1)*NRTP
            X1(JJ)=PSIRE(LL)
            Y1(JJ)=PSIIM(LL)
            S1=S1+X1(JJ)**2+Y1(JJ)**2
          ENDDO
          S1=S1*VOLEM
        ENDIF
        S2=0.0D0
        DO JJ=1,NRTP
          LL=JJ+(II-1)*NRTP
          X2(JJ)=PSIRE(LL)
          Y2(JJ)=PSIIM(LL)
          S2=S2+X2(JJ)**2+Y2(JJ)**2
        ENDDO
        S2=S2*VOLEM
        IF(JT+1.GT.JTR)THEN
          S3=0.0D0
          DO JJ=1,NRTP
            X3(JJ)=0.0D0
            Y3(JJ)=0.0D0
            S3=S3+X3(JJ)**2+Y3(JJ)**2
          ENDDO
          S3=S3*VOLEM
        ELSE
          S3=0.0D0
          DO JJ=1,NRTP
            LL=JJ+(IIP-1)*NRTP
            X3(JJ)=PSIRE(LL)
            Y3(JJ)=PSIIM(LL)
            S3=S3+X3(JJ)**2+Y3(JJ)**2
          ENDDO
          S3=S3*VOLEM
        ENDIF
!
        DO I=1,N_R
          DO J=1,N_T
            JJ=J+(I-1)*N_T
            DO K=1,N_P
              KK=K+(JJ-1)*N_P
              ARR1(I,J,K)=DCMPLX(X2(KK),Y2(KK))
              JP=NT2-J+1
              ARR1(I,JP,K)=-DCMPLX(X2(KK),Y2(KK))
            ENDDO
          ENDDO
        ENDDO
        CALL DFFTW_EXECUTE_DFT(PLAN1,ARR1,ARR1)
        DO I=1,N_R
          DO J=1,NT2
            DO K=1,N_P
              ARR(I,J,K)=ARR1(I,J,K)*AKR(I)*AKR(I)*SNG
              ART(I,J,K)=ARR1(I,J,K)*AKT(J)*AKT(J)*SNG
              ARP(I,J,K)=ARR1(I,J,K)*AKP(K)*AKP(K)*SNG
              ARP1(I,J,K)=-ARR1(I,J,K)*AKP(K)*SNG
            ENDDO
          ENDDO
        ENDDO
        CALL DFFTW_EXECUTE_DFT(PLAN2,ARR,ARR)
        CALL DFFTW_EXECUTE_DFT(PLAN3,ART,ART)
        CALL DFFTW_EXECUTE_DFT(PLAN4,ARP,ARP)
        CALL DFFTW_EXECUTE_DFT(PLAN5,ARP1,ARP1)
        VAL1=0.0D0
        VAL2=0.0D0
        VAL3=0.0D0
        VAL4=0.0D0
        VAL5=0.0D0
        VAL6=0.0D0
        VAL7=0.0D0
        VAL8=0.0D0
        VAL9=0.0D0
        VAL10=0.0D0
        DO I=1,N_R
          DO J=1,N_T
            JJ=J+(I-1)*N_T
            DO K=1,N_P
              KK=K+(JJ-1)*N_P
!
              HPSIRE1=0.5D0*RMST*DBLE(ARR(I,J,K))*SNG
              HPSIRE2=2.0D0*RMST*RHOD(I)*DBLE(ART(I,J,K))*SNG
              HPSIRE3=2.0D0*RMST*RHOD(I)*SIND(J)*DBLE(ARP(I,J,K))
     1               *WIN1(J)*SNG
              HPSIRE4=2.0D0*RMST*RHOD(I)*COSN(J)*SIND(J)
     1               *AK_J*DBLE(ARP1(I,J,K))*WIN1(J)*SNG
              HPSIRE5=0.5D0*RMST*RHOD(I)*SIND(J)*AK_J*AK_J
     1               *X2(KK)*WIN1(J)
              HPSIRE6=RMST*RHOD(I)*COSD(J)*EJK*X2(KK)*WIN2(J)
              HPSIRE7=-0.5D0*RMST*RHOD(I)
     1               *(0.25D0+4.0D0*SIN2D(J))*X2(KK)
              HPSIRE8=0.5D0*RMST*RHOD(I)*SINN(J)
     1               *COSD(J)*EJKM*X1(KK)*WIN2(J)
              HPSIRE9=0.5D0*RMST*RHOD(I)*SINN(J)
     1               *COSD(J)*EJKP*X3(KK)*WIN2(J)
              IF(IS.EQ.1)THEN
              HPSIRE10=V11(KK)*PRE(KK,JT,1)
     1                +V12(KK)*PRE(KK,JT,2)
     2                +V13(KK)*PRE(KK,JT,3)
              ENDIF
              IF(IS.EQ.2)THEN
              HPSIRE10=V21(KK)*PRE(KK,JT,1)
     1                +V22(KK)*PRE(KK,JT,2)
     2                +V23(KK)*PRE(KK,JT,3)
              ENDIF
              IF(IS.EQ.3)THEN
              HPSIRE10=V31(KK)*PRE(KK,JT,1)
     1                +V32(KK)*PRE(KK,JT,2)
     2                +V33(KK)*PRE(KK,JT,3)
              ENDIF
!
              HPSIIM1=0.5D0*RMST*DIMAG(ARR(I,J,K))*SNG
              HPSIIM2=2.0D0*RMST*RHOD(I)*DIMAG(ART(I,J,K))*SNG
              HPSIIM3=2.0D0*RMST*RHOD(I)*SIND(J)*DIMAG(ARP(I,J,K))
     1               *WIN1(J)*SNG
              HPSIIM4=2.0D0*RMST*RHOD(I)*COSN(J)*SIND(J)
     1               *AK_J*DIMAG(ARP1(I,J,K))*WIN1(J)*SNG
              HPSIIM5=0.5D0*RMST*RHOD(I)*SIND(J)*AK_J*AK_J
     1               *Y2(KK)*WIN1(J)
              HPSIIM6=RMST*RHOD(I)*COSD(J)*EJK*Y2(KK)*WIN2(J)
              HPSIIM7=-0.5D0*RMST*RHOD(I)
     1               *(0.25D0+4.0D0*SIN2D(J))*Y2(KK)
              HPSIIM8=0.5D0*RMST*RHOD(I)*SINN(J)
     1               *COSD(J)*EJKM*Y1(KK)*WIN2(J)
              HPSIIM9=0.5D0*RMST*RHOD(I)*SINN(J)
     1               *COSD(J)*EJKP*Y3(KK)*WIN2(J)
              IF(IS.EQ.1)THEN
              HPSIIM10=V11(KK)*PIM(KK,JT,1)
     1                +V12(KK)*PIM(KK,JT,2)
     2                +V13(KK)*PIM(KK,JT,3)
              ENDIF
              IF(IS.EQ.2)THEN
              HPSIIM10=V21(KK)*PIM(KK,JT,1)
     1                +V22(KK)*PIM(KK,JT,2)
     2                +V23(KK)*PIM(KK,JT,3)
              ENDIF
              IF(IS.EQ.3)THEN
              HPSIIM10=V31(KK)*PIM(KK,JT,1)
     1                +V32(KK)*PIM(KK,JT,2)
     2                +V33(KK)*PIM(KK,JT,3)
              ENDIF
!
              VAL1=VAL1+HPSIRE1*X2(KK)+HPSIIM1*Y2(KK)
              VAL2=VAL2+HPSIRE2*X2(KK)+HPSIIM2*Y2(KK)
              VAL3=VAL3+HPSIRE3*X2(KK)+HPSIIM3*Y2(KK)
              VAL4=VAL4+HPSIRE4*X2(KK)+HPSIIM4*Y2(KK)
              VAL5=VAL5+HPSIRE5*X2(KK)+HPSIIM5*Y2(KK)
              VAL6=VAL6+HPSIRE6*X2(KK)+HPSIIM6*Y2(KK)
              VAL7=VAL7+HPSIRE7*X2(KK)+HPSIIM7*Y2(KK)
              VAL8=VAL8+HPSIRE8*X2(KK)+HPSIIM8*Y2(KK)
              VAL9=VAL9+HPSIRE9*X2(KK)+HPSIIM9*Y2(KK)
              VAL10=VAL10+HPSIRE10*X2(KK)+HPSIIM10*Y2(KK)
            ENDDO
          ENDDO
        ENDDO
        VAL=VAL1+VAL2+VAL3+VAL4+VAL5+VAL6+VAL7+VAL8+VAL9+VAL10
        EER1(III)=VAL1*VOLEM/EVEPS
        EER2(III)=VAL2*VOLEM/EVEPS
        EER3(III)=VAL3*VOLEM/EVEPS
        EER4(III)=VAL4*VOLEM/EVEPS
        EER5(III)=VAL5*VOLEM/EVEPS
        EER6(III)=VAL6*VOLEM/EVEPS
        EER7(III)=VAL7*VOLEM/EVEPS
        EER8(III)=VAL8*VOLEM/EVEPS
        EER9(III)=VAL9*VOLEM/EVEPS
        EER10(III)=VAL10*VOLEM/EVEPS
        EER(III)=VAL*VOLEM/EVEPS
      ENDDO
!$OMP END DO 
!$OMP END PARALLEL 
      RETURN
      END
******************************************************************
*                                                                *
*    THIS ROUTINE CALCULATES THE VIB/ROT ENERGY (EPSILON) FOR    *
*    A DIATOMIC MOLECULE. INPUT: D(EPSILON), BETA(A(-1)),        *
*    RE(A), RM(AMU) AND THE VIB. AND ROT. QUANTUM NUMBERS (N,J). *
*                                                                *
******************************************************************
      FUNCTION EVIBRO(D,BETA,RE,RM,N,J)
      IMPLICIT REAL*8 (A-H,O-Z)
      HBAR=0.06350781278D+00
      HBAR2=0.5D0*HBAR*HBAR
********************************
*     NO LANGER CORRECTION     *
********************************
      RN=N*1.D0+0.5D0
      RJ=J*(J*1.D0+1.D+00)
      BE=HBAR2/(RM*RE*RE)
      DE=HBAR2**2/RM**2/BETA**2/RE**6/D
      XE=HBAR*BETA/DSQRT(8*RM*D)
      HME=HBAR*BETA*DSQRT(2*D/RM)
      ALE=1.5D0*HBAR2*HME*(1.D0-1.D0/(BETA*RE))/RM/BETA/RE**3/D
      EVIBRO=RJ*(BE-DE*RJ-ALE*RN)+HME*RN*(1.D0-XE*RN)
      RETURN
      END
CCCCCC
      FUNCTION EVIBRO3(D,BETA,RE,RM,N,J)
      IMPLICIT REAL*8 (A-H,O-Z)
      HBAR=0.06350781278D+00
      HBAR2=0.5D0*HBAR*HBAR
CCCCCCC  NO LANGER CORRECTION  CCCCCCC
      RN=N*1.D0+0.5D0
      RJ=J*(J*1.D0+1.D+00)
      BE=HBAR2/(RM*RE*RE)
      DE=HBAR2**2/RM**2/BETA**2/RE**6/D
      XE=HBAR*BETA/DSQRT(8*RM*D)
      HME=HBAR*BETA*DSQRT(2*D/RM)
      ALE=1.5D0*HBAR2*HME*(1.D0-1.D0/(BETA*RE))/RM/BETA/RE**3/D
      EVIBRO3=RJ*BE+HME*RN*(1.D0-XE*RN)
      RETURN
      END
CCCCCC
      SUBROUTINE INITIAL_WEIGHT(EKIN,ICH,DD,REQ,BETA,RM12,
     1KM,JR,RM123,AT,IDX,CKI,ET)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NSTP=16384)
      PARAMETER (NE=1000)
      DIMENSION DD(3),REQ(3),BETA(3),RM12(3),RM123(3)
      DIMENSION AT(NSTP),ET(NE),AKI(NE),CKI(NE)
      PI=4.0D0*DATAN(1.0D0)
      SIGMA=0.20D0
      HBAR=0.06350781278D0
      EVEPS=0.9648533822D0
      EK=EKIN*EVEPS
      XDE=DD(ICH)
      XBETA=BETA(ICH)
      XREQ=REQ(ICH)
      AMR=RM12(ICH)
      EVJ=EVIBRO(XDE,XBETA,XREQ,AMR,KM,JR)
      RMCM=RM123(ICH)
      AK_0=SQRT(2.0D0*RMCM*EK)/HBAR
      WRITE(166,*)
      WRITE(166,*)'AK_0',AK_0
      WRITE(166,*)
      WRITE(166,*)'WEIGHT FACTOR FOR INCOMING WAVE FUNCTION'
      WRITE(166,*)
      IE=0
      DO I=IDX,NE+IDX-1
        IE=IE+1
        ET(IE)=HBAR*AT(I)
        AKI(IE)=SQRT(2.0D0*RMCM*(ET(IE)-EVJ))/HBAR
        CKI(IE)=SQRT(2.0D0/PI)*SIGMA*(RMCM/AKI(IE)/HBAR/HBAR)
     1         *EXP(-2.0D0*SIGMA*SIGMA*(AKI(IE)-AK_0)**2)
        WRITE(166,'(I4,4X,3F16.6)')IE,ET(IE),AKI(IE),CKI(IE)
      ENDDO
      RETURN
      END
CCCCCC
      SUBROUTINE FINAL_WEIGHT(ICH,DD,REQ,BETA,RM12,
     1RM123,IDX,AT,AKF1,AKF2,AKF3,
     2NNV1,NNV2,NNV3,NNJ1,NNJ2,NNJ3)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NSTP=16384)
      PARAMETER (NE=1000)
      PARAMETER (NV=10)
      PARAMETER (NJ=12)
      DIMENSION AT(NSTP)
      DIMENSION DD(3),REQ(3),BETA(3),RM12(3),RM123(3)
      DIMENSION ET(NE),EVJ1(0:NV,0:NJ),EEE1(NE,0:NV,0:NJ)
      DIMENSION EVJ2(0:NV,0:NJ),EEE2(NE,0:NV,0:NJ)
      DIMENSION EVJ3(0:NV,0:NJ),EEE3(NE,0:NV,0:NJ)
      DIMENSION AKF1(NE,0:NV,0:NJ)
      DIMENSION AKF2(NE,0:NV,0:NJ)
      DIMENSION AKF3(NE,0:NV,0:NJ)
      DIMENSION NNV1(NE),NNJ1(NE,0:NV)
      DIMENSION NNV2(NE),NNJ2(NE,0:NV)
      DIMENSION NNV3(NE),NNJ3(NE,0:NV)
      EVEPS=0.9648533822D0
      HBAR=0.06350781278D0
      XDE1=DD(ICH)
      XBETA1=BETA(ICH)
      XREQ1=REQ(ICH)
      XDE2=2.7942D0*EVEPS
      XBETA2=1.3594D0
      XREQ2=1.0600D0
      XDE3=17.1020D0*EVEPS
      XBETA3=1.6456D0
      XREQ3=0.1821D0
      AMR=RM12(ICH)
      RMCM=RM123(ICH)
      IE=0
      DO I=IDX,NE+IDX-1
        IE=IE+1
        ET(IE)=HBAR*AT(I)
      ENDDO
!*****
      DO I=1,NE
        KK=0
        DO KM=0,NV
          JJ=0
          DO JR=0,NJ
            EVJ1(KM,JR)=EVIBRO(XDE1,XBETA1,XREQ1,AMR,KM,JR)
            EEE1(I,KM,JR)=ET(I)-EVJ1(KM,JR)
            IF(EEE1(I,KM,JR).LE.0.0D0) GO TO 10
              JJ=JJ+1
              AKF1(I,KM,JR)=DSQRT(2.0D0*RMCM*EEE1(I,KM,JR))
     1                    /HBAR/RMCM
          ENDDO
  10      CONTINUE
          NNV1(I)=KK
          NNJ1(I,KM)=JJ-1
          KK=KK+1
        ENDDO
      ENDDO
!*****
      DO I=1,NE
        KK=0
        DO KM=0,NV
          JJ=0
          DO JR=0,NJ
            EVJ2(KM,JR)=EVIBRO(XDE2,XBETA2,XREQ2,AMR,KM,JR)
            EEE2(I,KM,JR)=ET(I)-EVJ2(KM,JR)
            IF(EEE2(I,KM,JR).LE.0.0D0) GO TO 20
              JJ=JJ+1
              AKF2(I,KM,JR)=DSQRT(2.0D0*RMCM*EEE2(I,KM,JR))
     1                    /HBAR/RMCM
          ENDDO
  20      CONTINUE
          NNV2(I)=KK
          NNJ2(I,KM)=JJ-1
          KK=KK+1
        ENDDO
      ENDDO
!*****
      DO I=1,NE
        KK=0
        DO KM=0,NV
          JJ=0
          DO JR=0,NJ
            EVJ3(KM,JR)=EVIBRO3(XDE3,XBETA3,XREQ3,AMR,KM,JR)
            EEE3(I,KM,JR)=ET(I)-EVJ3(KM,JR)
            IF(EEE3(I,KM,JR).LE.0.0D0) GO TO 30
              JJ=JJ+1
              AKF3(I,KM,JR)=DSQRT(2.0D0*RMCM*EEE3(I,KM,JR))
     1                    /HBAR/RMCM
          ENDDO
  30      CONTINUE
          NNV3(I)=KK
          NNJ3(I,KM)=JJ-1
          KK=KK+1
        ENDDO
      ENDDO
!*****
      DO I=1,IE
        DO KM=0,NNV1(I)
          DO JR=0,NNJ1(I,KM)
            WRITE(166,'(I4,F12.6,2I4,2F12.6)')
     1      ICH,ET(I),KM,JR,EVJ1(KM,JR),AKF1(I,KM,JR)
          ENDDO
          WRITE(166,*)
        ENDDO
        WRITE(166,*)
        WRITE(166,*)
      ENDDO
        WRITE(166,*)
        WRITE(166,*) 'second state'
        WRITE(166,*)
      DO I=1,IE
        DO KM=0,NNV2(I)
          DO JR=0,NNJ2(I,KM)
            WRITE(166,'(I4,F12.6,2I4,2F12.6)')
     1      ICH,ET(I),KM,JR,EVJ2(KM,JR),AKF2(I,KM,JR)
          ENDDO
          WRITE(166,*)
        ENDDO
        WRITE(166,*)
        WRITE(166,*)
      ENDDO
        WRITE(166,*)
        WRITE(166,*) 'third state'
        WRITE(166,*)
      DO I=1,IE
        DO KM=0,NNV3(I)
          DO JR=0,NNJ3(I,KM)
            WRITE(166,'(I4,F12.6,2I4,2F12.6)')
     1      ICH,ET(I),KM,JR,EVJ3(KM,JR),AKF3(I,KM,JR)
          ENDDO
          WRITE(166,*)
        ENDDO
        WRITE(166,*)
        WRITE(166,*)
      ENDDO
!
      RETURN
      END
CCCCCC
      SUBROUTINE DEFAT(DT,AT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NSTP=16384)
      DIMENSION AT(NSTP)
      PI=4.0D0*DATAN(1.0D0)
      PI2=PI*2.0D0
      ACONST=PI2/(NSTP*DT)
      NHALFT=NSTP/2+1
      DO I=1,NHALFT
        AT(I)=ACONST*(I-1)
      ENDDO
      DO I=NHALFT+1,NSTP
        AT(I)=ACONST*(I-1-NSTP)
      ENDDO
      RETURN
      END
CCCCCC
C      SUBROUTINE PLOT(IST,RHO,THE,PHI,PSI)
C      IMPLICIT REAL*8(A-H,O-Z)
C      PARAMETER (N_R=128,N_T=64,N_P=128)
C      PARAMETER (NRTP=N_R*N_T*N_P)
C      PARAMETER (NSTATE=3)
C      PARAMETER (NSRTP=NSTATE*NRTP)
C      DIMENSION RHO(N_R),THE(N_T),PHI(N_P)
C      DIMENSION PSI(NSRTP)
C      DIMENSION RRH(N_R),PSIRHO1(N_R),PSIRHO2(N_R),PSIRHO3(N_R)
C       DO J=1,N_T
C        TH=THE(J)
C          DO K=1,N_P
C           PH=PHI(K)
C            DO I=1,N_R
C            RRH(I)=RHO(I)
C            KK=K+(J-1)*N_P+(I-1)*N_T*N_P
C            PSIRHO1(I)=PSI(KK)
C            PSIRHO2(I)=PSI(NRTP+KK)
C            PSIRHO3(I)=PSI(2*NRTP+KK)
C            ENDDO
C            CALL ORDER(N_R,RRH,PSIRHO1)
C            CALL ORDER(N_R,RRH,PSIRHO2)
C            CALL ORDER(N_R,RRH,PSIRHO3)
C            RH=RRH(N_R)
C            PS1=PSIRHO1(N_R)
C            PS2=PSIRHO2(N_R)
C            PS3=PSIRHO3(N_R)
C            BETA=DSIN(TH)*DCOS(PH)
C            GAMA=DSIN(TH)*DSIN(PH)
C           WRITE(1399+IST,'(6F12.6)')TH,PH,RH,PS1,PS2,PS3
C           WRITE(1499+IST,'(6F12.6)')BETA,GAMA,RH,PS1,PS2,PS3
C          CALL FLUSH(1399+IST)
C          CALL FLUSH(1499+IST)
C          ENDDO
C          WRITE(1399+IST,*)
C          WRITE(1499+IST,*)
C          CALL FLUSH(1399+IST)
C          CALL FLUSH(1499+IST)
C       ENDDO
C      RETURN
C      END
CCCCCC
C      SUBROUTINE ORDER(N,X,A)
C      IMPLICIT REAL*8(A-H,O-Z)
C      DIMENSION X(N),A(N)
C      DO  40  L=1,N
C      ATEST=1.0E+12
C      DO 41 J=L,N
C      IF(A(J)-ATEST)42,41,41
C   42 ATEST=A(J)
C      XTEST=X(J)
C      JTEST=J
C   41 CONTINUE
C      A(JTEST)=A(L)
C      X(JTEST)=X(L)
C      A(L)=ATEST
C      X(L)=XTEST
C   40 CONTINUE
C      RETURN
C      END
CCCCCC
      SUBROUTINE splin(x,y,n,yp1,ypn,xa,ya)
      INTEGER n
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (N_R=128)
      INTEGER i,k,khi,klo
      REAL*8 p,qn,sig,un,u(n_r)
      REAL*8 a,b,h,xa,ya
      if (yp1.gt.0.99e30) then
           y2(1)=0.0D0
           u(1)=0.0D0
      else
           y2(1)=-0.5D0
           u(1)=(3.0D0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
           sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
           p=sig*y2(i-1)+2.0D0
           y2(i)=(sig-1.0D0)/p
           u(i)=(6.0D0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt.0.99e30) then
           qn=0.0D0
           un=0.0D0
      else
           qn=0.5D0
           un=(3.0D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))
     1              /(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0D0)
      do k=n-1,1,-1
           y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
CCC
      klo=1
      khi=n
   1  if (khi-klo.gt.1) then
           k=(khi+klo)/2
           if(x(k).gt.xa)then
                khi=k
           else
                klo=k
           endif
      goto 1
      endif
      h=x(khi)-x(klo)
      if (h.eq.0.0D0) stop 'bad x input in splin'
      a=(x(khi)-xa)/h
      b=(xa-x(klo))/h
      ya=a*y(klo)+b*y(khi)+
     1           ((a**3-a)*y2(klo)+(b**3-b)*y2(khi))*(h**2)/6.0D0
      return
      END
CCC
      SUBROUTINE SPLIN1(X,Y,N,YP1,YPN,XA,YA)
      INTEGER n
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (N_R=128)
      INTEGER i,k,khi,klo
      REAL*8 p,qn,sig,un,u(n_r)
      REAL*8 a,b,h,xa,ya
!$OMP PARALLEL
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(x,y,n,yp1,ypn,xa,ya)
      if (yp1.gt.0.99e30) then
           y2(1)=0.0D0
           u(1)=0.0D0
      else
           y2(1)=-0.5D0
           u(1)=(3.0D0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
           sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
           p=sig*y2(i-1)+2.0D0
           y2(i)=(sig-1.0D0)/p
           u(i)=(6.0D0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     1           /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt.0.99e30) then
           qn=0.0D0
           un=0.0D0
      else
           qn=0.5D0
           un=(3.0D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))
     1              /(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0D0)
      do k=n-1,1,-1
           y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      klo=1
      khi=n
   1  if (khi-klo.gt.1) then
           k=(khi+klo)/2
           if(x(k).gt.xa)then
                khi=k
           else
                klo=k
           endif
      goto 1
      endif
      h=x(khi)-x(klo)
      if (h.eq.0.0D0) stop 'bad x input in splin'
      a=(x(khi)-xa)/h
      b=(xa-x(klo))/h
      ya=a*y(klo)+b*y(khi)+
     1           ((a**3-a)*y2(klo)+(b**3-b)*y2(khi))*(h**2)/6.0D0
!$OMP END PARALLEL
      return
      END
 
c #####################################################
c
c Multi-valued DMBE potential energy surface reported in
c      L.P. Viegas, A. Alijah & A.J.C. Varandas
c          J. Chem. Phys. 126, 074309 (2007)

c For method, see:
c  A.J.C. Varandas, in Advanced Series in Physical Chemistry, 
c  Special Issue on Conical Intersections: Electronic
c  Structure, Spectroscopy and Dynamics, edited by W. Domcke,
c  D.R. Yarkony and H.  Koppel, (World Scientific, 2004), 
c  Ch. 5., pp. 205-270 (and references therein)

c r1,r2,r3 --> interparticle distances (a0)
c pot      --> 3 x 3 Diabatic Matrix (Eh)

c Call potential with:
c call singlet(r1,r2,r3,pot)
c ######################################################
      

      SUBROUTINE SINGLET(R1,R2,R3,POT)
      implicit none
      integer np,inst,nst,nparam,npol,order
      double precision x0,gam1,r0,beta,pot(3,3)
      double precision tb(4),r1,r2,r3
      common/polorder/npol(4),order(3,4)
      common/switch/gam1(3,4),x0(3,4),r0(3,4),beta(3,4)
      common/param/np(4)
               
      npol(1)=2
      gam1(1,1)=0.3d0
      gam1(2,1)=0.3d0
      r0(1,1)=1.65d0
      r0(2,1)=1.65d0
      beta(1,1)=1.3d0
      beta(2,1)=1.3d0
      x0(1,1)=15.0d0
      x0(2,1)=12.0d0
      order(1,1)=7
      order(2,1)=11

      npol(2)=1
      gam1(1,2)=0.3d0
      r0(1,2)=2.0d0
      beta(1,2)=1.0d0
      x0(1,2)=7.0d0
      order(1,2)=10
      
      npol(3)=1
      gam1(1,3)=0.3d0
      r0(1,3)=2.5d0
      beta(1,3)=1.0d0
      x0(1,3)=10.0d0
      order(1,3)=5
      
      npol(4)=1
      gam1(1,4)=0.3d0
      r0(1,4)=2.5d0
      beta(1,4)=1.0d0
      x0(1,4)=12.0d0
      order(1,4)=7
      
      
      np(1)=nparam(order(1,1))+nparam(order(2,1))
      np(2)=nparam(order(1,2))
      np(3)=nparam(order(1,3)) 
      np(4)=nparam(order(1,4)) 
                           
      inst=1
      
      do nst=1,4
c calculates three-body terms        
        call thrbody(r1,r2,r3,tb(nst),nst,inst)
        
        inst=inst+np(nst)
        
      enddo                 
      
      call dimpot(r1,r2,r3,tb,pot)
      
      return
      end
      
      subroutine thrbody(d1,d2,d3,v,nst,inst)
      implicit none
      double precision A,Q,q1,q2,q3,V,r,x
      double precision g1p,g2p,g3p,damp,rho
      double precision gam1,x0,r0,d1,d2,d3,beta
      integer i,k,j,l,num,nr,order,npol,nst,inst,np
      dimension R(3),A(3,3),Q(3)
      common/polorder/npol(4),order(3,4)
      common/switch/gam1(3,4),x0(3,4),r0(3,4),beta(3,4)
      common/param/np(4)      
      common/coefs/x(228)

      R(1)=d1
      R(2)=d2
      R(3)=d3
      
      rho=dsqrt((R(1)**2+R(2)**2+R(3)**2)/dsqrt(3.0d0))
      A(1,1)=dsqrt(1.d0/3.d0)
      A(1,2)=A(1,1)
      A(1,3)=A(1,1)
      A(2,1)=0.d0
      A(2,2)=dsqrt(1.d0/2.d0)
      A(2,3)=-A(2,2)
      A(3,1)=dsqrt(2.d0/3.d0)
      A(3,2)=-dsqrt(1.d0/6.d0)
      A(3,3)=A(3,2)
      
      num=inst-1
      
      
      V=0.0d0
      
      do nr=1,npol(nst)
        
        do 100 i=1,3
          Q(i)=0.d0
          do 100 j=1,3
            Q(i)=Q(i)+A(i,j)*
     &        (1.d0-dexp(-beta(nr,nst)*(R(j)/R0(nr,nst)-
     &        1.d0)))/beta(nr,nst)
 100      continue
          
          q1=Q(1)
          q2=Q(2)**2+Q(3)**2
          q3=Q(3)**3-3.d0*Q(3)*Q(2)**2
          
          do 101 l=0,order(nr,nst)
            do 101 i=0,l
              if (i.eq.0) then
                g1p=1.0d0
              else
                g1p=q1**i
              end if
              do 101 j=0,(l-i),2
                if (j.eq.0) then
                  g2p=1.0d0
                else
                  g2p=q2**(j/2)
                end if
                k=l-i-j
                if(mod(k,3).EQ.0) then
                  if (k.eq.0) then
                    g3p=1.0d0
                  else
                    g3p=q3**(k/3)
                  end if
                  num=num+1
                  V=V+x(num)*g1p*g2p*g3p*
     &              damp(gam1(nr,nst),rho,x0(nr,nst))
                endif
                
 101          continue
              
            enddo  ! end of npol loop
            
            V=V*(1.0d0-damp(10.0d0,d1,0.8d0))*
     &        (1.0d0-damp(10.0d0,d2,0.8d0))*
     &        (1.0d0-damp(10.0d0,d3,0.8d0))
            
            return
          end
          


c  Builds dressed DIM matrix 
      subroutine dimpot(r1,r2,r3,tb,pot)
      implicit none
      double precision r1,r2,r3,xh2pd,ah2pd,pothhx,G
      double precision tb(4),pot(3,3)
      
      pot(1,1)=xh2pd(r2)+ah2pd(r2)+xh2pd(r3)+ah2pd(r3)
      pot(1,1)=pot(1,1)+(tb(2)+tb(3))*2.0d0
      pot(1,1)=0.5d0*pot(1,1)+pothhx(r1)+tb(1)-1.0d0
      pot(1,1)=pot(1,1)+G(r1,r2,r3)
      
      pot(2,2)=xh2pd(r1)+ah2pd(r1)+xh2pd(r3)+ah2pd(r3)
      pot(2,2)=pot(2,2)+(tb(2)+tb(3))*2.0d0
      pot(2,2)=0.5d0*pot(2,2)+pothhx(r2)+tb(1)-1.0d0
      pot(2,2)=pot(2,2)+G(r1,r2,r3)
      
      pot(3,3)=xh2pd(r1)+ah2pd(r1)+xh2pd(r2)+ah2pd(r2)
      pot(3,3)=pot(3,3)+(tb(2)+tb(3))*2.0d0
      pot(3,3)=0.5d0*pot(3,3)+pothhx(r3)+tb(1)-1.0d0
      pot(3,3)=pot(3,3)+G(r1,r2,r3)
   
      pot(1,2)=0.5d0*(xh2pd(r3)-ah2pd(r3)-tb(4)**2)
      pot(2,1)=pot(1,2)
      
      pot(1,3)=0.5d0*(xh2pd(r2)-ah2pd(r2)-tb(4)**2)
      pot(3,1)=pot(1,3)

      pot(2,3)=0.5d0*(xh2pd(r1)-ah2pd(r1)-tb(4)**2)
      pot(3,2)=pot(2,3)

      return
      end


      FUNCTION POTHHX(R)
C===============================================================
C  NEWEST EHFACE2U FIT FOR H2 GROUND STATE POTENTIAL CURVE 
c  ##2006-06-22##
C  USING AB INITIO POINTS FROM:
C  L. WOLNIEWICZ JCP 99(3),1851 (1993) ---> FIRST POINTS
C  L. WOLNIEWICZ JCP 103(5),1792 (1995) ---> CORRECTIONS    
C  RANGE of R: 0.6 to 12 a0
C
C  Initial fit  
C  RMS(m= 52 )= 35.2987773353410503  cm-1
C  Final fit
C  RMS(m= 52 )= 0.967169906620670844E-01  cm-1
C
C  For method, see:
C  A.J.C. Varandas, S.P. Rodrigues  V.M.O. Batista
C  Chem. Phys. Lett., 424, 425-431 (2006), and references therein;
C  A.J.C. Varandas, Adv. Chem. Phys., 74, 255 (1988)
C===============================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DD(20)
      COMMON/POTEN/D1(20),D2(20),C(20),GAMMA,AGP,AI(9),R0,RM
      COMMON/LIM/NLOW,NUPP
      COMMON/TH/G0,G1,G2
      COMMON/EXPOE/RE,IEXP
C      COMMON/OUT/VHF,ASEXC,DISP,DD(20)
      COMMON/ASYEXC/CATILD,ATILD(2),ALPHT,GTILD
      NLOW=6
      NUPP=16
      DO II=NLOW,NUPP
        D1(II)=AN(II)
        D2(II)=BN(II)
      ENDDO
      CATILD=-0.8205D0
      ATILD(1)=0.0D0
      ATILD(2)=0.0D0
      ALPHT=2.5D0
      GTILD=2.0D0
      IEXP=1
      RE=0.14010000D+01
      C(6)=0.64990000D+01
      C(7)=0.0D0
      C(8)=0.12440000D+03
      C(9)=0.0D0
      C(10)=0.32858000D+04
      C(11)=-0.34750000D+04
      C(12)=0.12150000D+06
      C(13)=-0.29140000D+06
      C(14)=0.60610000D+07
      C(15)=-0.23050000D+08
      C(16)=0.39380000D+09
      GAMMA=2.5D0
      AGP=0.229794389784158d0
      AI(1)=1.74651398886700093d0
      AI(2)=0.631036031819560028d0
      AI(3)=0.747363488024733624d0
      AI(4)=0.956724297662875783d-01
      AI(5)=0.131320504483065703d0
      AI(6)=-0.812200084994067194d-07
      AI(7)=0.119803887928360935d-01
      AI(8)=-0.212584227748381302d-02
      AI(9)=0.509125901134908042d-03
      G0=1.02072511539524680d0
      G1=1.82599688484061118d0
      G2=0.269916332495592104d0
      R0=0.69282032D+01
      RM=0.11000000D+02
      X=R-RE

       pol1=ai(9)
       do i=8,1,-1              
         pol1=pol1*x + ai(i)
       end do
       pol1=pol1*x + 1.0d0 

      GAM=G0*(1.0D0+G1*TANH(G2*X))
      VHF=-AGP/(R**IEXP)*POL1*EXP(-GAM*X)
    
      ASEXC=1.0D0
      DO I=1,2
        ASEXC=ASEXC+ATILD(I)*R**I
      ENDDO
      ASEXC=ASEXC*CATILD*R**ALPHT*EXP(-GTILD*R)
    
      RHH=0.5D0*(RM+GAMMA*R0)
      X=R/RHH
      DEXC=(1.0D0-EXP(-D1(NLOW)*X-D2(NLOW)*X**2))**NLOW
    
      ASEXC=ASEXC*DEXC
      VHF=VHF+ASEXC
    
      DISP=0.0D0
      DO 1 I=NLOW,NUPP
        DAMPI=(1.0D0-EXP(-D1(I)*X-D2(I)*X**2))**I
        DD(I)=DAMPI
        DISP=DISP-C(I)*DAMPI*R**(-I)
 1    CONTINUE
      POTHHX=VHF+DISP
      RETURN
      END

           
      real*8 function ah2pd(r)
C===============================================================
C POTENTIAL CURVE FOR H2+ ( A 2^sigma^(+)_(u) )
C ## 2006-06-22 ##
C USING AB INITIO POINTS FROM:
C J. M. PEEK, JCP 43(9), 3004 (1965)
C RANGE of R: 3.5 to 15 a0 
C
C  Final fit  
C  RMS(m= 24 )= 0.125196743261995869  cm-1
C  For method, see:
C  A.J.C. Varandas, J. Chem. Phys. 107, 867 (1997)
C===============================================================
      implicit none
      integer i
      real*8 r,coef(0:7),xh2pd,v
      coef( 0 )=  1.11773285795729826d0
      coef( 1 )= -1.27592697554394174d0
      coef( 2 )=  0.235612064424508216d0
      coef( 3 )= -0.500203729467869895d-01
      coef( 4 )=  0.568627052480373801d-02
      coef( 5 )= -0.382978465312642114d-03
      coef( 6 )=  0.149149267032670154d-04
      coef( 7 )= -0.267518847221239873d-06
c
      v=coef(7)
      do i=6,0,-1
        v=v*r+coef(i)
      enddo
      ah2pd=dexp(v)+xh2pd(r)
      return
      end



      function XH2PD(R)
C===============================================================
C EHFACE2U POTENTIAL CURVE FOR H2+ ( X 2^sigma^(+)_(g) )
C ## 2006-06-22 ##
C USING AB INITIO POINTS FROM:
C D. M. BISHOP AND R. W. WETMORE, MOL. PHYS. 26(1),145 (1972)
C RANGE of R: 0.6 to 10 a0
C
C  Initial fit  
C  RMS(m= 95 )= 67.6787429191265630  cm-1
C  Final fit  
C  RMS(m= 95 )= 0.675527097585787786E-02  cm-1

C  For method, see:
C  A.J.C. Varandas, S.P. Rodrigues  V.M.O. Batista
C  Chem. Phys. Lett., 424, 425-431 (2006), and references therein;
C  A.J.C. Varandas, Adv. Chem. Phys., 74, 255 (1988)
C===============================================================     
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DD(20)
      COMMON/POTEN1/D1(20),D2(20),C(20),GAMMA,AGP,AI(11),R0,RM
      COMMON/LIM1/NLOW,NUPP
      COMMON/TH1/G0,G1,G2
      COMMON/EXPOE1/RE,IEXP
      COMMON/ASYEXC1/CATILD,ATILD(2),ALPHT,GTILD
      NLOW=4
      NUPP=11
      DO II=NLOW,NUPP
        D1(II)=AN(II)
        D2(II)=BN(II)
      ENDDO
      IEXP=1
      RE=2.0D0
      C(4)=2.250d0
      C(5)=0.0d0
      C(6)=7.5000d0
      C(7)=53.25d0
      C(8)=65.625d0
      C(9)=886.5d0
      C(10)=1063.125d0
      C(11)=21217.5d0
      GAMMA=2.5D0
      AGP=-0.180395506614312d0
      AI(1)=-0.897676265677028185d0
      AI(2)=-0.771599228853606545d0
      AI(3)=-0.245766669963638718d0
      AI(4)=-0.788889284685244524d-01
      AI(5)=-0.252032558464952844d-01
      AI(6)=0.681894227468654839d-02
      AI(7)=0.655940943163255976d-03
      AI(8)=-0.531288172311992135d-03
      AI(9)=0.890418306898401330d-04
      AI(10)=-0.666314834544138477d-05
      AI(11)=0.194182431833699709d-06
      G0=0.960151039243191562d0
      G1=-0.353946173857859037d0
      G2=-0.496213155382123572d0
      R0=3.4641d0
      RM=0.11000000D+02
      X=R-RE

      pol1=ai(11)
      do i=10,1,-1
         pol1=pol1*x + ai(i)
      end do
      pol1=pol1*x + 1.0d0

      GAM=G0*(1.0D0+G1*TANH(G2*X))
      VHF=-AGP/(R**IEXP)*POL1*EXP(-GAM*X)
C
      RHH=0.5D0*(RM+GAMMA*R0)
      X=R/RHH    
C
      DISP=0.0D0
      DO 1 I=NLOW,NUPP
      DAMPI=(1.0D0-EXP(-D1(I)*X-D2(I)*X**2))**I
      DD(I)=DAMPI
      DISP=DISP-C(I)*DAMPI*R**(-I)
    1 CONTINUE
      XH2PD=VHF+DISP
      RETURN
      END


      FUNCTION AN(N)
      IMPLICIT REAL*8(A-H,O-Z)
      ALPH0=16.36606D0
      ALPH1=0.70172D0
      AN=ALPH0/(FLOAT(N))**ALPH1
      RETURN
      END


      FUNCTION BN(N)
      IMPLICIT REAL*8(A-H,O-Z)
      BET0=17.19338D0
      BET1=0.09574D0
      BN=BET0*EXP(-BET1*FLOAT(N))
      RETURN
      END


      function G(r1,r2,r3)
      implicit none
      double precision r1,r2,r3,G,rh,ph,th
      double precision gauss(4),r(3)
      r(1)=r1
      r(2)=r2
      r(3)=r3
      call inte2hyper(r,rh,ph,th)            
      
      gauss(1)=
     &  dexp(-1d4*(r1-1.401d0)**2)+
     &  dexp(-1d4*(r2-1.401d0)**2)+
     &  dexp(-1d4*(r3-1.401d0)**2)
      
      gauss(2)=
     & -0.000195d0*dexp(-0.03d0*(rh-17.0d0)**2)-
     &  0.00099d0*dexp(-0.15d0*(rh-10.0d0)**2)+
     &  0.000252d0*dexp(-0.25d0*(rh-11.0d0)**2)-
     &  0.000138d0*dexp(-0.01d0*(rh-16.0d0)**2)+
     &  0.000078d0*dexp(-0.03d0*(rh-20.0d0)**2)+
     &  0.000003d0*dexp(-0.01d0*(rh-26.0d0)**2)-
     &  0.000135d0*dexp(-0.1d0*(rh-9.0d0)**2)+
     &  0.000005d0*dexp(-0.1d0*(rh-28.0d0)**2)+
     &  0.000002d0*dexp(-0.1d0*(rh-31.0d0)**2)
      
      gauss(3)=
     &  dexp(-1d4*(r1-2.0d0)**2)+
     &  dexp(-1d4*(r2-2.0d0)**2)+
     &  dexp(-1d4*(r3-2.0d0)**2)
      
      gauss(4)=
     & -0.00002d0*dexp(-0.03d0*(rh-19.0d0)**2)+
     &  0.000008d0*dexp(-0.05d0*(rh-21.0d0)**2)+
     &  0.000002d0*dexp(-0.01d0*(rh-25.0d0)**2)
      
      
      G=gauss(1)*gauss(2)+gauss(3)*gauss(4)
      
      return
      end
      

      BLOCK DATA
c coeficients for singlet H3+ surface
      implicit none
      double precision x
      common/coefs/x(228)
      data x/-0.2845459958197D+00, 0.3850540915997D+01,
     & -0.1362727983983D+02,-0.1763127575795D+01,
     &  -0.3997956751079D+02,-0.9478227108607D+01,
     &  -0.6256864319761D+00, 0.6922830117563D+02,
     &  -0.5544008250600D+01,-0.1223243581975D+02,
     &  -0.5565035493903D+00,-0.7990835144213D+02,
     &  0.4937987578000D+01, 0.1032846680854D+01,
     &  -0.1024944938355D+01, 0.8171934328106D-01,
     &  -0.3400630003098D+03, 0.6364663539559D+02,
     &  -0.1294217695606D+03,-0.7241144393993D+02,
     &  -0.5238580628055D+01, 0.1645346915236D+02,
     &  -0.4559177543910D+00,-0.1114852144762D+03,
     &  -0.4160076873226D+02, 0.5646176507749D+02,
     &  0.1927845463688D+03, 0.1034016582997D+02,
     &  -0.2495720708364D+01,-0.2196094831605D+01,
     &  0.4995040605555D+00, 0.2204571641363D+01,
     &  -0.5439237573590D+00, 0.6846403692672D-02,
     &  -0.2277735103524D+01, 0.2759089283649D+02,
     &  0.7936127182406D+00,-0.5871241896630D+00,
     &  -0.1323872092546D+03,-0.1005323892600D+02,
     &  0.1582657184728D+02,-0.1541686382410D+01,
     &  0.5087330884991D+02,-0.6751737369193D+02,
     &  -0.1262521382543D+02, 0.7798725873929D+01,
     &  0.2037607478530D+01, 0.3962102074881D+03,
     &  -0.5555238520074D+02, 0.9548022399071D+02,
     &  0.5696848300022D+02,-0.1947033826314D+02,
     &  -0.1832745716879D+02, 0.8890681073646D+00,
     &  0.1411483776813D+03, 0.7978735114995D+02,
     &  -0.1117319392262D+03,-0.2388272564512D+03,
     &  0.3374095081649D+02, 0.2094562992129D+02,
     &  -0.4003037330580D+01,-0.2856224060853D+00,
     &  0.5881794112860D+01,-0.8982996737073D+01,
     &  0.8046224130870D+01, 0.5281960633581D-01,
     &  -0.3546907369872D+01,-0.6528334773919D+01,
     &  -0.3830505287670D+01,-0.6386051039291D+00,
     &  0.3162332951830D+01,-0.2223083229855D-01,
     &  -0.4491167884772D+00,-0.1188967934541D+01,
     &  0.3489210439438D+01,-0.4869763268661D+01,
     &  0.2175303912647D+01,-0.5191497797669D+01,
     &  0.7828626803682D+01,-0.1128407061209D+01,
     &  -0.1271808573208D+01, 0.3069274516791D+00,
     &  0.1854243578919D+01,-0.8535061206250D-01,
     &  0.5833990722727D+00,-0.1055455848541D+01,
     &  -0.5919355026412D+00, 0.4226971822154D+01,
     &  0.8892733627208D+01, 0.2065687224306D+01,
     &  0.1917895395954D+02,-0.4133476161782D+00,
     &  0.1614975694526D+02, 0.8669361117049D+01,
     &  0.1408495804579D+01, 0.3686771253329D+00,
     &  -0.7860440980033D-01,-0.7080531777858D-01,
     &  0.2785512538230D+00,-0.1016421721250D+00,
     &  0.4176560040398D+01,-0.1568656118207D+01,
     &  0.8007659100381D+00, 0.8788486668834D+01,
     &  0.1268922736683D+02, 0.2261940416494D+01,
     &  0.1795296373183D+02, 0.3515372129617D+01,
     &  0.9815017191585D+01, 0.7562210497515D+01,
     &  0.2722248355571D+01, 0.4294261733772D+00,
     &  0.7691463873126D+00, 0.2481557620244D-02,
     &  0.6026978066922D+00, 0.3295481061855D+01,
     &  0.7015291733970D+01,-0.2387729507482D+01,
     &  -0.6533271953599D+00,-0.1252544784685D+02,
     &  0.1912931010490D+01, 0.6134022859639D+02,
     &  0.1342268001485D+02,-0.2625380523555D+01,
     &  -0.1899038446199D+01, 0.2656123342914D+01,
     &  -0.4389413075873D+02,-0.2729075588350D+01,
     &  0.5828467674841D+01,-0.1520452027360D+01,
     &  -0.4598017470650D+02,-0.1756604916681D+02,
     &  -0.1596671552906D+02,-0.5703262507979D+02,
     &  -0.1600677373244D+02, 0.3881629920764D+01,
     &  0.1614190302676D+01, 0.4952213787993D+01,
     &  0.5812444415360D+02, 0.3511472345525D+02,
     &  0.5099606433991D+01, 0.3049558978461D+02,
     &  -0.5297861997608D+01,-0.6873510179327D+01,
     &  0.1040913734067D+00, 0.1197981571767D+02,
     &  0.4865789229401D+01,-0.7019610038457D+01,
     &  -0.2624867301367D+01,-0.5637503355720D+02,
     &  -0.4716014331310D-01, 0.4098674152271D+02,
     &  0.2104582076046D+02,-0.7065446618818D-01,
     &  -0.1793715098910D+00, 0.1470149524446D+01,
     &  -0.2723240374922D+01,-0.7702762844985D+01,
     &  -0.1085691653328D+02, 0.7878956333729D+01,
     &  -0.1933229959665D+02, 0.4110444708049D+02,
     &  0.9498387330543D+00,-0.3407816727441D+02,
     &  -0.1085082072732D+02, 0.3459417287786D+01,
     &  -0.6360621370667D-01, 0.9314527643683D+00,
     &  -0.2167692738050D+01,-0.1276412781004D+01,
     &  0.1410484578285D+01,-0.3940181087105D+01,
     &  0.1179565629874D+02,-0.1774547832707D+01,
     &  0.7072060093826D+01,-0.1127319089089D+02,
     &  0.1714444782390D+00, 0.7394073327918D+01,
     &  0.1704664527444D+01,-0.1192526565174D+01,
     &  0.2821887364885D-01,-0.3413951575886D+01,
     &  0.1360890787929D+01, 0.8633198410350D+01,
     &  0.4246887005837D+01, 0.9535669384044D+01,
     &  -0.3578154147373D+01,-0.3106769965434D+01,
     &  -0.1004190932219D+01,-0.1123219690535D+02,
     &  -0.5735177734152D+01, 0.1410002623784D+00,
     &  -0.9955301221383D+00,-0.6025440086335D-01,
     &  0.3423771693923D+01, 0.3040350360659D+01,
     &  0.1955754187339D+00,-0.1156220573976D+00,
     &  0.8224556646293D-01, 0.1399371093644D+00,
     &  -0.1318539326582D+00,-0.3842144298995D+00,
     &  -0.4121345438816D+00,-0.1794233281375D-01,
     &  -0.9300050560412D+00,-0.6364712743012D+00,
     &  -0.2602285731992D+00, 0.2229809064704D+00,
     &  0.3062049976404D+00, 0.6744482883143D+00,
     &  0.4858749722998D-01, 0.7297018494878D+00,
     &  -0.5934091032128D-01, 0.2055510177770D+00,
     &  0.1155389702307D+01, 0.1466606391956D+01,
     &  0.6468339605757D+00, 0.3077172810685D+00,
     &  -0.1173646603994D+00,-0.6609279132694D-01,
     &  0.5619836307431D+00,-0.3858574645853D+00,
     &  0.7136936266049D+00,-0.3012225426016D+00,
     &  -0.7345473025896D+00,-0.3083986346985D+00,
     &  -0.2505835865913D+00, 0.2537569624720D-01/
      end
           

c SUBROUTINE THAT TRANSFORMS FROM INTERNAL COORDINATES
c TO HYPERSPHERICAL COORDINATES, USING MASS-SCALED JACOBI COORDINATES
c r is the vector of internal coordinates
c rho is given in au. phi and theta in degrees      
c 0 <= phi <= 4pi , 0 <= theta <= pi/2
      subroutine inte2hyper(r,rho,phi,theta)
      implicit none
      double precision xmass(3),mtot,mu,d(3),rsmall,phi,sphi,ct,angphi
      double precision r1,rbig,r2,rho2,a,xcos,cthet,sthet2,sthet,cphi2
      double precision r(3),rs(3),w(3),theta,rho,ccita(3),cphi,pre

c nuclear masses
      xmass(1)=1.007276d0
      xmass(2)=1.007276d0
      xmass(3)=1.007276d0

      w(1)=xmass(1)
      w(2)=xmass(2)
      w(3)=xmass(3)

c mass scaling factors for jacobi system k=1
      
      mtot=xmass(1)+xmass(2)+xmass(3)
      mu=sqrt(xmass(1)*xmass(2)*xmass(3)/mtot)      
      d(1)=dsqrt( (1.0d0-xmass(1)/mtot)*xmass(1)/mu )

c Deals with problem at linear configurations
c THETA=90 degrees
      pre=5.0d-9
      IF(abs(r(1)+r(2)-r(3)).le.pre.or.abs(r(1)+r(3)-r(2)).
     &    le.pre.or.abs(r(2)+r(3)-r(1)).le.pre) THEN
    
      rho=dsqrt(r(1)**2+r(2)**2+r(3)**2)*3.0d0**(-0.25d0)
      
      theta=90.0d0
      sphi=2.0d0*r(1)**2/(d(1)**2*rho**2)-1.0d0
      cphi2=(r(2)**2-r(3)**2)/rho**2       
               
c Avoid instability in sphi calculation
      if(sphi.gt.1.0d0) then
        sphi=1.0d0
      endif

      ct=cphi2/(sphi+5.0d-12)
      phi=angphi(sphi,ct)*180.0d0/dacos(-1.0d0)
   
      ELSE
      
      call INT_JAC(R,W,RS,CCITA)

      xcos=ccita(1)
      r1=r(1)
      r2=rs(1)
      
c mass scaling of jacobi coordinates
c r1 is the distance between atoms 2 and 3
c r2 is the distance between atom 1 and the diatom (2,3)
      rsmall=r1/d(1)
      rbig=r2*d(1)
c conversion mass-scaled jacobi to hyperspherical coordinates
      rho2=rsmall**2+rbig**2
      a=0.5d0*rsmall*rbig*dsqrt(1.0d0-xcos**2)      
      cthet=4.0d0*a/rho2
      
c Sets a random value for phi at the pole, as phi is undefined there
c THETA=0 degrees
      
c Avoid instability in cthet**2 with gt. eq is for theta=0
      if(cthet**2.ge.1.0d0) then
        theta=0.0d0
        phi=30.0d0
        rho=dsqrt(rho2)
      else
        
c General case, THETA =\= 0 and 90 degrees
        
      sthet2=1.0d0-cthet**2
      sthet=dsqrt(sthet2+5.0d-12)
      rho=dsqrt(rho2)
      theta=dacos(cthet)*180.0d0/dacos(-1.0d0)       
      cphi=2.0d0*rho**(-2.0d0)*rsmall*rbig*xcos/sthet
      sphi=dsqrt(1.0d0-cphi**2)
      ct=cphi/(sphi+5.0d-12)
      phi=angphi(sphi,ct)*180.0d0/dacos(-1.0d0)
    
      endif

      ENDIF
      
      return
      end


C Esta subrutina toma las coordenadas internas de una molecula
C triatomica y halla las correspondientes coordenadas de Jacobi
C Los atomos son enumerados 1,2,3, de masas w1,w2,w3 y las distancias
C son R(1)=D(1,2) R(2)=D(1,3) R(3)=D(2,3)
C Todos los vectores son bidimensionales Z=0
      SUBROUTINE INT_JAC(R,W,RS,CCITA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(3),W(3),RS(3),CCITA(3),X1(2),X2(2),X3(2),
     &  CM12(2),CM13(2),CM23(2)

      DO I=1,2
        X1(I)=0.0D0
        X2(I)=0.0D0
      ENDDO
      X2(1)=R(1)
      CALPH=(R(1)**2+R(2)**2-R(3)**2)/(2.0D0*R(1)*R(2))
      SALPH=DSQRT(1.0D0-CALPH**2+5.d-12)
      X3(1)=R(2)*CALPH
      X3(2)=R(2)*SALPH
      CALL CMASA(X1,X2,W(1),W(2),CM12)
      CALL CMASA(X1,X3,W(1),W(3),CM13)
      CALL CMASA(X2,X3,W(2),W(3),CM23)
      RS(1)=DIST(X3,CM12)
      RS(2)=DIST(X2,CM13)
      RS(3)=DIST(X1,CM23)
      CALL ANG(X3,X2,CM12,CCITA(1))
      CALL ANG(X2,X1,CM13,CCITA(2))
      CALL ANG(X1,X3,CM23,CCITA(3))
      RETURN
      END

C SUBRUTINA PARA CALCULAR EL CENTRO DE MASA
      SUBROUTINE CMASA(X1,X2,W1,W2,CM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X1(2),X2(2),CM(2)
      DO I=1,2
        CM(I)=(W1*X1(I)+W2*X2(I))/(W1+W2)
      ENDDO

      RETURN
      END

C CALCULAR EL COSENO DEL ANGULO ENTRE LOS VECTORES (A-C) Y (B-C)
      SUBROUTINE ANG(A,B,C,COSENO)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2),B(2),C(2),D(2),E(2)
      DO I=1,2
        D(I)=A(I)-C(I)
        E(I)=B(I)-C(I)
      ENDDO
      DD=P_ESC(D,D)
      EE=P_ESC(E,E)
      COSENO=P_ESC(D,E)/(DSQRT(DD*EE)+5.0D-12)
      RETURN
      END

C FUNCION PRODUCTO ESCALAR ENTRE A LOS VECTORES A Y B
      DOUBLE PRECISION FUNCTION P_ESC(A,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2),B(2)
      SUMA=0.0D0
      DO I=1,2
        SUMA=SUMA+A(I)*B(I)
      ENDDO
      P_ESC=SUMA
      RETURN
      END

C FUNCION DISTANCIA ENTRE DOS PUNTOS
      DOUBLE PRECISION FUNCTION DIST(X1,X2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X1(2),X2(2)
      DIST=DSQRT((X1(1)-X2(1))**2+(X1(2)-X2(2))**2)
      RETURN
      END


c     =====================================
      double precision function angphi (si,ct)
      implicit real*8 (a-h,o-z)
c     =====================================
c
c     input : si     sinus
c             ct     cotangens
c
c     output: winkel im bereich
c             von  0 bis 2*pi
c
      data pi/3.14159265358979323846264d0/,zpi/6.28318530717958647692529
     *d0/
      theta=dasin(si)
      if (si*ct.gt.0.0d0) go to 10
      angphi=pi-theta
      return
10    angphi=theta
      if (ct.ge.0.0d0) return
      angphi=zpi+theta
      end

          
      real*8 function damp(gam,q,qq)
c range determining factor
      real*8 q,gam,qq
      damp=1.d0/(1.d0+dexp(gam*(q-qq)))
      return
      end
        
        
      integer function nparam (order)
      implicit none
c calculates number of parameters up to order "order"
      integer i,j,k,l,order
      nparam=0
      do l=0,order
        do i=0,l
          do j=0,(l-i),2
            k=l-i-j
            if(mod(k,3).eq.0) then
              nparam=nparam+1
            endif
            end do
          end do
        end do
        return
      end
CCCCCC
      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
      implicit real*8(a-h,o-z)
C
      INTEGER N,NM,IERR,MATZ
      dimension A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL SYMMETRIC MATRIX.
C
C     ON INPUT
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT.
C
C        N  IS THE ORDER OF THE MATRIX  A.
C
C        A  CONTAINS THE REAL SYMMETRIC MATRIX.
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED.  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER.
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO.
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN ERROR
C           COMPLETION CODE DESCRIBED IN THE DOCUMENTATION FOR TQLRAT
C           AND TQL2.  THE NORMAL COMPLETION CODE IS ZERO.
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     .......... FIND EIGENVALUES ONLY ..........
c      CALL  TRED1(NM,N,A,W,FV1,FV2)
c      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
      END
C
C********************************************************************
C
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
      implicit real*8(a-h,o-z)
C
      INTEGER I,J,K,L,N,II,NM,JP1
       dimension a(NM,N),D(N),E(N),Z(NM,N)
c:      REAL F,G,H,HH,SCALE
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.
C
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX
C          PRODUCED IN THE REDUCTION.
C
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      DO 100 I = 1, N
C
         DO 80 J = I, N
   80    Z(J,I) = A(J,I)
C
         D(I) = A(N,I)
  100 CONTINUE
C
      IF (N .EQ. 1) GO TO 510
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0d0
         SCALE = 0.0d0
         IF (L .LT. 2) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + dABS(D(K))
C
         IF (SCALE .NE. 0.0d0) GO TO 140
  130    E(I) = D(L)
C
         DO 135 J = 1, L
            D(J) = Z(L,J)
            Z(I,J) = 0.0d0
            Z(J,I) = 0.0d0
  135    CONTINUE
C
         GO TO 290
C
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         F = D(L)
         G = -dSIGN(dSQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
C     .......... FORM A*U ..........
         DO 170 J = 1, L
  170    E(J) = 0.0d0
C
         DO 240 J = 1, L
            F = D(J)
            Z(J,I) = F
            G = E(J) + Z(J,J) * F
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + Z(K,J) * D(K)
               E(K) = E(K) + Z(K,J) * F
  200       CONTINUE
C
  220       E(J) = G
  240    CONTINUE
C     .......... FORM P ..........
         F = 0.0d0
C
         DO 245 J = 1, L
            E(J) = E(J) / H
            F = F + E(J) * D(J)
  245    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM Q ..........
         DO 250 J = 1, L
  250    E(J) = E(J) - HH * D(J)
C     .......... FORM REDUCED A ..........
         DO 280 J = 1, L
            F = D(J)
            G = E(J)
C
            DO 260 K = J, L
  260       Z(K,J) = Z(K,J) - F * E(K) - G * D(K)
C
            D(J) = Z(L,J)
            Z(I,J) = 0.0E0
  280    CONTINUE
C
  290    D(I) = H
  300 CONTINUE
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0d0
         H = D(I)
         IF (H .EQ. 0.0d0) GO TO 380
C
         DO 330 K = 1, L
  330    D(K) = Z(K,I) / H
C
         DO 360 J = 1, L
            G = 0.0d0
C
            DO 340 K = 1, L
  340       G = G + Z(K,I) * Z(K,J)
C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * D(K)
  360    CONTINUE
C
  380    DO 400 K = 1, L
  400    Z(K,I) = 0.0d0
C
  500 CONTINUE
C
  510 DO 520 I = 1, N
         D(I) = Z(N,I)
         Z(N,I) = 0.0d0
  520 CONTINUE
C
      Z(N,N) = 1.0d0
      E(1) = 0.0d0
      RETURN
      END
C
C********************************************************************
C
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C
      implicit real*8(a-h,o-z)
C
       dimension D(N),E(N),Z(NM,N)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1.
C
C        E HAS BEEN DESTROYED.
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F =    0.0D0
      TST1 = 0.0D0
      E(N) = 0.0D0
C
      DO 240 L = 1, N
         J = 0
         H =dABS(D(L)) + dABS(E(L))
         IF (TST1 .LT. H) TST1 = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            TST2 = TST1 +dABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 121
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  121    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0d0 * E(L))
         R = PYTHAG(P,1.0d0)
         D(L) = E(L) / (P + dSIGN(R,P))
         D(L1) = E(L) * (P + dSIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0d0
         C2 = C
         EL1 = E(L1)
         S = 0.0d0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            R = PYTHAG(P,E(I))
            E(I+1) = S * R
            S = E(I) / R
            C = P / R
            P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         TST2 = TST1 + dABS(E(L))
         IF (TST2 .GT. TST1) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
C     SORTING START
C      DO 300 II = 2, N
C         I = II - 1
C         K = I
C         P = D(I)
CC
C         DO 260 J = II, N
C            IF (D(J) .GE. P) GO TO 260
C            K = J
C            P = D(J)
C  260    CONTINUE
CC
C         IF (K .EQ. I) GO TO 300
C         D(K) = D(I)
C         D(I) = P
CC
C         DO 280 J = 1, N
C            P = Z(J,I)
C            Z(J,I) = Z(J,K)
C            Z(J,K) = P
C  280    CONTINUE
CC
C  300 CONTINUE
C     SORTING START
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
C
C********************************************************************
C
      REAL*8 FUNCTION PYTHAG(A,B)
C
      implicit real*8(a-h,o-z)
C
C     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
      P = dMAX1(dABS(A),daBS(B))
      IF (P .EQ. 0.0d0) GO TO 20
      R = (dMIN1(dABS(A),dABS(B))/P)**2
   10 CONTINUE
         T = 4.0d0 + R
         IF (T .EQ. 4.0d0) GO TO 20
         S = R/T
         U = 1.0d0 + 2.0d0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
C
C********************************************************************
