! vector_LU_decomposition FORTRAN programm,
! FORTRAN program With implemented eye(), triu() tril(), matrix_print() functions 


 Program vector_LU_decomposition
    implicit none
    
    interface ! обеспечиваем явный интерфейс функции eye
    function eye(sz) result(res)
            integer, intent(in) :: sz
            real(8) :: res(sz,sz)
            integer :: t
        end function
    end interface
    
    
    interface ! обеспечиваем явный интерфейс функции tri
    ! Реализует triu и tril путем изменения 3го аргумента больше или меньше нуля соотв.
    ! 2й аргумент - размер одного измерения матрицы
    function tri(a,sz2,dir,offset) result(l)
        real(8), intent(in) :: a(:,:)
        integer, intent(in) :: sz2,dir,offset
        real(8) :: l(sz2,sz2)
        integer :: r,w
        end function
    end interface
    
    interface !  обеспечиваем явный интерфейс функции print_matrix
    ! предназначена для вывода матриц
    function print_matrix(P,sz3) result(void)
        real(8), intent(in) :: P(:,:)
        integer, intent(in) :: sz3
        integer :: d,s,void
        end function
    end interface
    
    !Объявляем переменные
    integer, parameter :: N=500
    integer :: i , j, k ,v
    real(8) :: A(N,N)
    real(8) :: B(N,N)
    real(8) :: L(N,N)
    real(8) :: U(N,N)
    real(8) :: C(N,N)
    real(8) :: conv
    real(8) :: time_start
    real(8) :: time_end
    CALL CPU_TIME(time_start)
    CALL RANDOM_NUMBER(A) 
    A=A+N*eye(N) ! заполнение матрицы А
    B = A
    
    Print *, "Start calculation"
    ! Реализация алгоритма LU разложения
    do k=1,N-1
        A(k+1:N,k)=A(k+1:N,k)/A(k,k);
        do i=k+1,N
            A(i,k+1:n)=A(i,k+1:n)-A(i,k)*A(k,k+1:n);
        enddo
    enddo
    
    !U=triu(A,0); L=tril(A,-1)+eye(n); ! получение U и L матриц
    !norm(B-L*U) !проверка сходимости
    
    L=tri(A,N,-1,-1)+eye(N)
    U=tri(A,N,1,0)
    C=B-matmul(L,U)
    conv = norm2(C)
    ! конечное время
    CALL CPU_TIME(time_end)
    !  вывод матрицы
    !v=print_matrix(A,N)
    !v=print_matrix(L,N)
    !v=print_matrix(U,N)
    !v=print_matrix(C,N)
    
    Print *, "Result(convergence):"
    write (*,*) conv
    Print '("Time:" f6.3," seconds")', (time_end - time_start)
    End Program vector_LU_decomposition
    
    
    function eye(sz) result(res)
            integer, intent(in) :: sz
            real(8) :: res(sz,sz)
            integer :: t
            res = 0
            do t = 1, sz
                res(t,t) = 1
            end do
        end function
    
    
    function tri(a,sz2,dir,offset) result(l)
        real(8), intent(in) :: a(:,:)
        integer, intent(in) :: sz2,dir,offset
        real(8) :: l(sz2,sz2)
        integer :: r,w
        if(dir .GE. 0) then
            do r = 1,sz2
                    do w = 1,sz2
                            if(w .GE. r+offset) then   
                                l(r,w) = a(r,w)
                            else
                                l(r,w) = 0
                            end if
                        
                end do 
            end do
        else
            do r = 1,sz2
                    do w = 1,sz2
                            if(w .LE. r+offset) then   
                                l(r,w) = a(r,w)
                            else
                                l(r,w) = 0
                            end if
                        
                end do 
            end do
        end if
    end function 
    
function print_matrix(P,sz3) result(void)
    real(8), intent(in) :: P(:,:)
    integer, intent(in) :: sz3
    integer :: d,s,void
    void = 0
    Print *, "Matrix:"
        do, d=1,sz3
            write(*,*) ( P(d,s), s=1,sz3)
        enddo    
    end function
