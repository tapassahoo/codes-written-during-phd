	implicit none
	include 'mpif.h'
	integer n,m,l
	parameter(n=10,m=10,l=10)
	double precision a(n,m),b(m,l),c(n,l),ck
	integer i,j,k,ierr,numtasks,rank,nsend,nrecvr,ij,nrecv,ia,ja,ija
	integer status(MPI_STATUS_SIZE)
	call MPI_INIT(ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	if(rank.eq.0)then
	  do i=1,n
	    do j=1,m
	      a(i,j)=i+j
	    enddo
	  enddo
	  do i=1,m
	    do j=1,l
	      b(i,j)=i*j
	    enddo
	  enddo
	  write(*,*)'a'
	  do i=1,n
	    write(*,*)(a(i,j),j=1,m)
	  enddo
	  write(*,*)'b'
	  do i=1,m
	    write(*,*)(b(i,j),j=1,l)
	  enddo
	  write(*,*)
	endif
	call MPI_BCAST(a,n*m,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(b,m*l,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	if(rank.eq.0)then
	  nsend=0
	  nrecvr=1
	  nrecv=0
	  do i=1,n
	    do j=1,l
	      ij=(i-1)*l+j
	      call MPI_SEND(ij,1,MPI_INTEGER,nrecvr,0,MPI_COMM_WORLD,
     $                      ierr)
	      nsend=nsend+1
	        if(nrecvr.lt.numtasks-1)nrecvr=nrecvr+1
	    enddo
	  enddo
	  do i=nrecv+1,n*l
	    call MPI_RECV(ck,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,
     $               MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
	    nrecv=nrecv+1
	    ija=status(MPI_TAG)
	    ia=((ija-1)/l)+1
	    ja=ija-(ia-1)*l
	    c(ia,ja)=ck
	  enddo
	else
 1	  call MPI_RECV(ij,1,MPI_INTEGER,0,MPI_ANY_TAG,
     $             MPI_COMM_WORLD,status,ierr)
	  ija=status(MPI_TAG)
	  if(ija.eq.0)then
	    i=((ij-1)/l)+1
	    j=ij-(i-1)*l
	    ck=0.0
	    do k=1,m
	      ck=ck+a(i,k)*b(k,j)
	    enddo
	    call MPI_SEND(ck,1,MPI_DOUBLE_PRECISION,0,ij,MPI_COMM_WORLD,
     $                    ierr)
	    goto 1
	  elseif(ija.eq.n*l+1)then
	    goto 2
	  endif
	endif
	if(rank.eq.0)then
	  do i=1,numtasks-1
	    call MPI_SEND(0,1,MPI_INTEGER,i,n*l+1,MPI_COMM_WORLD,ierr)
	  enddo
	  write(*,*)'c'
	  do i=1,n
	    write(*,*)(c(i,j),j=1,l)
	  enddo
	endif
 2	call MPI_FINALIZE(ierr)
	stop
	end
