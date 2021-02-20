! vector_LU_decomposition FORTRAN programm,
! FORTRAN program With implemented eye(), triu() tril(), matrix_print() functions 


 Program vector_LU_decomposition
    implicit none
    
    interface ! обеспечиваем явный интерфейс функции eye
    function eye(sz) result(res)
            integer, intent(in) :: sz
            real :: res(sz,sz)
            integer :: t
        end function
    end interface
    
    
    interface ! обеспечиваем явный интерфейс функции tri
    ! Реализует triu и tril путем изменения 3го аргумента больше или меньше нуля соотв.
    ! 2й аргумент - размер одного измерения матрицы
    function tri(a,sz2,dir,offset) result(l)
        real, intent(in) :: a(:,:)
        integer, intent(in) :: sz2,dir,offset
        real :: l(sz2,sz2)
        integer :: r,w
        end function
    end interface
    
    interface !  обеспечиваем явный интерфейс функции print_matrix
    ! предназначена для вывода матриц
    function print_matrix(P,sz3) result(void)
        real, intent(in) :: P(:,:)
        integer, intent(in) :: sz3
        integer :: d,s,void
        end function
    end interface
    
    !Объявляем переменные
    integer, parameter :: N=4
    integer :: i , j, k ,v
    real :: A(N,N)
    real :: B(N,N)
    real :: L(N,N)
    real :: U(N,N)
    real :: C(N,N)
    real :: conv
    
    CALL RANDOM_NUMBER(A) 
    A=A+N*eye(N) ! заполнение матрицы А
    B = A
    
    Print *, "Start calculation"
    ! Реализация алгоритма LU разложения
    do k=1,N-1
        A(k+1:N,k)=A(k+1:N,k)/A(k,k);
        do j=k+1,N
            A(k+1:N,j)=A(k+1:N,j)-A(k+1:N,k)*A(k,j);
        enddo
    enddo
    
    !U=triu(A,0); L=tril(A,-1)+eye(n); ! получение U и L матриц
    !norm(B-L*U) !проверка сходимости
    
    L=tri(A,N,-1,-1)+eye(N)
    U=tri(A,N,1,0)
    C=B-matmul(L,U)
    conv = norm2(C)
    
    !  вывод матрицы
    v=print_matrix(A,N)
    v=print_matrix(L,N)
    v=print_matrix(U,N)
    v=print_matrix(C,N)
    
    Print *, "Result(convergence):"
    write (*,*) conv
    
    End Program vector_LU_decomposition
    
    
    function eye(sz) result(res)
            integer, intent(in) :: sz
            real :: res(sz,sz)
            integer :: t
            res = 0
            do t = 1, sz
                res(t,t) = 1
            end do
        end function
    
    
    function tri(a,sz2,dir,offset) result(l)
        real, intent(in) :: a(:,:)
        integer, intent(in) :: sz2,dir,offset
        real :: l(sz2,sz2)
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
    real, intent(in) :: P(:,:)
    integer, intent(in) :: sz3
    integer :: d,s,void
    void = 0
    Print *, "Matrix:"
        do, d=1,sz3
            write(*,*) ( P(d,s), s=1,sz3)
        enddo    
    end function
