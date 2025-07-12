!qr factorisarion
!using gram smith orthogonalisation
!verified

module m_mod
    implicit none
    save
    real,parameter::pi=3.141592
    contains

    subroutine write_matrix(a)
            !--------------------------------------------------------------------------
            !this subroutine prints n/n matrix 'a' on screen
            !--------------------------------------------------------------------------
            implicit none
            real,dimension(:,:)::a 
            integer::i,j
            write(*,*)
                do i=lbound(a,1),ubound(a,1)
                    write(*,*) (a(i,j), j=lbound(a,2),ubound(a,2))
                enddo
            return
        end subroutine write_matrix


    subroutine qr_decompose(m,n,a,q,r)
    !**********************************************************************************
    !This subroutine performs qr decomposition using gram schimidt orthogonalisation process
    ![input]~   m,n is the dimension of [m * n] input ' a' matrix
    !           a is the input matrix
    ![output]~  q is the orthogonal matrix
    !           r is the upper triangular matrix
    !**********************************************************************************
        implicit none
        integer::m,n,k,i,i1,j1
        real,dimension(1:m,1:n)::a,q,u,e,b      !a is the main [m by n] matrix;; q is the orthogonal [m by n] matrix
        real,dimension(1:n,1:n)::r              !r is the upper triangular [n by n] matrix
        !-------------------------put zero for all matrix (except 'a')-------------------------
        e=0.0
        u=0.0
        b=0.0
        q=0.0
        r=0.0
        !-------------------------first calculate e(:,1) manually-------------------------
        u(:,1)=a(:,1)                   !copy the first column of 'a' to first column of 'u'
        e(:,1)=u(:,1)/norm2(u(:,1))     !normalise u(:,1)
              
        !-------------------------iterate for the other e's-------------------------
        do k=1,n-1
            b=0.0                                               !put b=0.0 each time
            do i=1,k
                b(:,i)=(dot_product(a(:,k+1),e(:,i)))*e(:,i)    !calculate (a.e)*e and stoe it in 'b' columnwise
            enddo

            u(:,k+1)=a(:,k+1)-sum(b, dim=2)                     !find u
            e(:,k+1)=u(:,k+1)/norm2(u(:,k+1))                   !normalise u
        enddo
        !-------------------------transfer e to q matrix-------------------------
        q=e
        !-------------------------find the r matrix-------------------------

        do j1=1,n 
            do i1=1,j1 
            r(i1,j1)=dot_product(a(:,j1),e(:,i1))
            enddo
        enddo

        return
    end subroutine qr_decompose

end module m_mod


program main_matrix
    use m_mod
    implicit none
    integer,parameter::m=3,n=3      ![m by n] is the dimension of the a matrix
    real,dimension(1:m,1:n)::a,q    !a is the main [m by n] matrix;; q is the orthogonal [m by n] matrix
    real,dimension(1:n,1:n)::r      !r is the upper triangular [n by n] matrix
    real, dimension(1:m*n)::xx      !here write the input matrix(x) 
    real,dimension(1:m)::b          !'b' is the RHS of the set of equations

    !xx=[1.0,1.0,1.0,1.0,2.0,4.0,1.0,3.0,9.0]   !enter the matrix elements row-wise
    xx=[1,1,0,1,0,1,0,1,1]
    !xx=[-1,-1,1,1,3,3,-1,-1,5,1,3,7]
    b=[3,4,6]

    a=reshape(xx,[m,n],order=[2,1]) !constructed the main matrix
    print*, "main matrix::"
    call write_matrix(a)

    call qr_decompose(m,n,a,q,r)
    print*, "orthogonal matrix::"
    call write_matrix(q)
    print*, "upper triangular matrix::"
    call write_matrix(r)

end program
