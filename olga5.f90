 PROGRAM JACOBS
 Implicit none
 include 'mpif.h'
 Integer status(MPI_STATUS_SIZE);
 integer ierr,myid, size;
 real, allocatable :: U(:,:); 
 Integer N, T, Q, iQ
 PARAMETER(N=50, T=1000)
 Real  h, pi,x,y,tau,c, L
 Integer I,J, it, Ibeg,Jbeg,Iend,Jend
 real sec,secB,secE
 integer B(10),E(10)    !вектора для значений времени

 pi=3.14159265358979
 L=5
 h=L/(Real(N)-1.0)
 tau=0.00001; ! шаги по времени и пространству
 c=tau/h/h; ! отношение шагов
 Q=4; ! высота пирамиды

! начало работы с mpi
call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)


allocate (U(1:N,1:N/2+Q)); 


if ( myid==0 )  then

! задание начального условия
 do i=1,N
   x=(float(i)-1.0)*h;
   do j=1,N/2+Q
     y=(float(j)-1.0)*h;
     U(i,j)=sin(pi*x/L)*sin(pi*y/L); ! начальное условие
   enddo
 enddo
  
  call date_and_time(values = B) !засекаем время
  Do it = 1,T/Q
    call MPI_SEND(U(2:N-1,N/2-Q+1:N/2),Q*N-2*Q,MPI_REAL,1,it,MPI_COMM_WORLD, ierr);
    do iQ=Q-1,0,-1
      U(2:N-1,2:N/2+iQ)=U(2:N-1,2:N/2+iQ)+c*(U(1:N-2,2:N/2+iQ)+U(3:N,2:N/2+iQ)+ &
                        U(2:N-1,1:N/2+iQ-1)+U(2:N-1,3:N/2+iQ+1)-4*U(2:N-1,2:N/2+iQ)); ! слой пирамиды iQ
    enddo  
    call MPI_RECV(U(2:N-1,N/2+1:N/2+Q),Q*N-2*Q,MPI_REAL,1,it,MPI_COMM_WORLD, status,ierr);
  enddo
 
  call date_and_time(values = E) !засекаем время
  secB = B(5)*3600+B(6)*60+B(7)+B(8)*0.001
  secE = E(5)*3600+E(6)*60+E(7)+E(8)*0.001		
  sec = secE-secB
  print *, "Q= ",Q, "N= ",N,"T= ",T, "  Running time (in sec.): ",sec

else

! задание начального условия
 do i=1,N
   x=(float(i)-1.0)*h;
   do j=1,N/2+Q
     y=L/2-h*(Q-0.5)+(float(j)-1.0)*h; 
     U(i,j)=sin(pi*x/L)*sin(pi*y/L); ! начальное условие
   enddo
 enddo

  Do it = 1,T/Q
    call MPI_SEND(U(2:N-1,Q+1:2*Q),Q*N-2*Q,MPI_REAL,0,it,MPI_COMM_WORLD, ierr);
    do iQ=2,Q+1
      U(2:N-1,iQ:N/2+Q-1)=U(2:N-1,iQ:N/2+Q-1)+c*(U(1:N-2,iQ:N/2+Q-1)+U(3:N,iQ:N/2+Q-1)+ &
                          U(2:N-1,iQ-1:N/2+Q-2)+U(2:N-1,iQ+1:N/2+Q)-4*U(2:N-1,iQ:N/2+Q-1)); ! слой пирамиды iQ
    enddo  
    call MPI_RECV(U(2:N-1,1:Q),Q*N-2*Q,MPI_REAL,0,it,MPI_COMM_WORLD, status,ierr);
  enddo

!  print*, myid, U(N/2,Q+1:N/2+Q);
endif

deallocate (U); 

call MPI_FINALIZE(ierr);


end
