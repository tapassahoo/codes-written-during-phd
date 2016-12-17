      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'mpif.h'
      PARAMETER (NN=128)
      PARAMETER (MM=128)
      PARAMETER (NN1=128)
      PARAMETER (MM1=128)
      PARAMETER (NP=9)
      PARAMETER (KKJMP=8)
      DIMENSION XX(NN),YY(MM)
      DIMENSION YYA(NN,MM,NP),Y2A(NN,MM)
      DIMENSION XA(NN1),YA(MM1)
      DIMENSION YAA(NN,MM),TAUR(NN1,MM1,NP)
****************************************************************
      INTEGER STATUS(MPI_STATUS_SIZE)                          !
      DOUBLE PRECISION TRR(NN1,MM1,KKJMP)                      !
      CALL MPI_INIT(IERR)                                      !
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROC,IERR)          !
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)             !
      MASTER=0                                                 !
****************************************************************
      IF(MYID.EQ.MASTER)THEN
        OPEN(1,FILE='DIAH-0.7_128.DAT',STATUS='OLD')
        PI=4.0D0*DATAN(1.0D0)
        DO I=1,MM
          DO J=1,NN
            READ(1,*)XX(I),YY(J),(YYA(I,J,K),K=1,NP)
          ENDDO
          READ(1,*)
        ENDDO
        XMIN=XX(1)
        XMAX=XX(NN)
        YMIN=YY(1)
        YMAX=YY(MM)
        DDX=(XMAX-XMIN)/DFLOAT(NN1)
        DDY=(YMAX-YMIN)/DFLOAT(MM1)
        DO I=1,NN1
          XA(I)=XMIN+(I-1)*DDX
        ENDDO
        DO I=1,MM1
          YA(I)=YMIN+(I-1)*DDY
        ENDDO
      ENDIF
      IF(MOD(NP,KKJMP).EQ.0)THEN
        NTOT=NP/KKJMP
      ELSE
        NTOT=NP/KKJMP+1
      ENDIF
****************************************************************
      CALL MPI_BCAST(YYA,NN*MM*NP,MPI_DOUBLE_PRECISION,MASTER  !  
     1              ,MPI_COMM_WORLD,IERR)                      !  
      CALL MPI_BCAST(XX,NN,MPI_DOUBLE_PRECISION,MASTER         !  
     1              ,MPI_COMM_WORLD,IERR)                      ! MPI BROADCAST
      CALL MPI_BCAST(YY,MM,MPI_DOUBLE_PRECISION,MASTER         ! FROM MASTER
     1              ,MPI_COMM_WORLD,IERR)                      ! NODE
      CALL MPI_BCAST(XA,NN1,MPI_DOUBLE_PRECISION,MASTER        ! 
     1              ,MPI_COMM_WORLD,IERR)                      ! 
      CALL MPI_BCAST(YA,MM1,MPI_DOUBLE_PRECISION,MASTER        ! 
     1              ,MPI_COMM_WORLD,IERR)                      ! 
****************************************************************
      IF(MYID.EQ.MASTER)THEN
        NSENT=1
        KKIN=1
  230 CONTINUE
        DO NRCVR=1,NUMPROC-1
          ITAG=0
          CALL MPI_SEND(KKIN,1,MPI_INTEGER
     1                 ,NRCVR,ITAG,MPI_COMM_WORLD,IERR)
          KKFNL=KKIN+KKJMP-1
          IF(KKFNL.GT.NP) KKFNL=NP
          ITAG=1
          CALL MPI_SEND(KKFNL,1,MPI_INTEGER
     1                 ,NRCVR,ITAG,MPI_COMM_WORLD,IERR)
          KKIN=KKIN+KKJMP
          NSENT=NSENT+1
        ENDDO
        DO II=1,NUMPROC-1
          CALL MPI_RECV(TRR,NN1*MM1*KKJMP,MPI_DOUBLE_PRECISION
     1        ,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,STATUS,IERR)
          ITAG1=STATUS(MPI_TAG)
          NRCVR1=STATUS(MPI_SOURCE)
          KKIN1=ITAG1
          KKFNL1=KKIN1+KKJMP-1
          IF(KKFNL1.GT.NP) KKFNL1=NP
          DO K=KKIN1,KKFNL1
            DO I=1,NN1
              DO J=1,MM1
                TAUR(I,J,K)=TRR(I,J,K-KKIN1+1)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF(NSENT.LE.NTOT)GOTO 230
        ITAG=0
        DO NRCVR=1,NUMPROC-1
          CALL MPI_SEND(-1,1,MPI_INTEGER
     1                 ,NRCVR,ITAG,MPI_COMM_WORLD,IERR)
        ENDDO
      ELSE
  232   ITAG=0
        CALL MPI_RECV(KKIN,1,MPI_INTEGER
     1               ,MASTER,ITAG,MPI_COMM_WORLD,STATUS,IERR)
        IF(KKIN.LT.0)GOTO 231
        ITAG=1
        CALL MPI_RECV(KKFNL,1,MPI_INTEGER
     1               ,MASTER,ITAG,MPI_COMM_WORLD,STATUS,IERR)
****************************************************************
*                                                              *
*     OPENMP PARALLELIZATION FOR ONLY OUTERMOST K - LOOP       *
*     YYA,XX,YY,XA,YA,TAUR - SHARED VARIABLE                   *
*                                                              *
****************************************************************
        ICHUNK=1
        CALL OMP_SET_NUM_THREADS(KKJMP)
!$OMP PARALLEL 
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(KKIN,KKFNL,YYA,XX,YY,XA,YA,TRR)
!$OMP DO SCHEDULE(DYNAMIC,ICHUNK)
        DO K=KKIN,KKFNL
          DO I=1,NN
            DO J=1,MM
              YAA(I,J)=YYA(I,J,K)
            ENDDO
          ENDDO
          CALL splie2(xx,yy,yaa,nn,mm,y2a)
          DO I=1,NN1
            X1=XA(I)
            DO J=1,MM1
              Y1=YA(J)
              CALL splin2(xx,yy,yaa,y2a,nn,mm,x1,y1,yout)
              TRR(I,J,K-KKIN+1)=YOUT
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO 
!$OMP END PARALLEL
        ITAG=KKIN
        CALL MPI_SEND(TRR,NN1*MM1*KKJMP,MPI_DOUBLE_PRECISION
     1               ,MASTER,ITAG,MPI_COMM_WORLD,IERR)
      GOTO 232
      ENDIF
  231 CONTINUE
      IF(MYID.EQ.MASTER)THEN
        DO I=1,NN1
          DO J=1,MM1
            WRITE(20,55)XA(I),YA(J),(TAUR(I,J,K),K=1,NP)
          ENDDO
          WRITE(20,*)
        ENDDO
      ENDIF
  55  FORMAT(1X,12F14.6)
      CALL MPI_FINALIZE(IERR)
      STOP
      END
CCC
      SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NN=128)
      DIMENSION x1a(m),x2a(n),y2a(m,n),ya(m,n)
      DIMENSION y2tmp(NN),ytmp(NN)
!$OMP PARALLEL 
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(x1a,x2a,ya,m,n,y2a)
      XX1=1.d30
      YY1=1.d30
      do j=1,m
      do k=1,n
      ytmp(k)=ya(j,k)
      enddo
      call spline(x2a,ytmp,n,XX1,YY1,y2tmp)
      do k=1,n
      y2a(j,k)=y2tmp(k)
      enddo
      enddo
!$OMP END PARALLEL
      return
      END
CCC
      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NN=128)
      DIMENSION x1a(m),x2a(n),y2a(m,n),ya(m,n)
      DIMENSION y2tmp(NN),ytmp(NN),yytmp(NN)
!$OMP PARALLEL 
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      XX1=1.d30
      YY1=1.d30
      do j=1,m
      do k=1,n
      ytmp(k)=ya(j,k)
      y2tmp(k)=y2a(j,k)
      enddo
      call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
      enddo
      call spline(x1a,yytmp,m,XX1,YY1,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
!$OMP END PARALLEL
      return
      END
CCC
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=128)
      DIMENSION x(n),y(n),y2(n)
      DIMENSION u(NMAX)
!$OMP PARALLEL 
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(x,y,n,yp1,ypn,y2)
      if (yp1.gt..99e30) then
      y2(1)=0.
      u(1)=0.
      else
      y2(1)=-0.5
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
      u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     */(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99e30) then
      qn=0.
      un=0.
      else
      qn=0.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
!$OMP END PARALLEL
      return
      END
CCC
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION xa(n),y2a(n),ya(n)
!$OMP PARALLEL 
!$OMP$ DEFAULT(PRIVATE)
!$OMP$ SHARED(xa,ya,y2a,n,x,y)
      klo=1
      khi=n
    1 if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
      khi=k
      else
      klo=k
      endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) STOP "bad xa input in splint"
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     * ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
!$OMP END PARALLEL
      return
      END
