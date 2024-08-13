program printH
    implicit none
    character:: jobu*1, jobvt*1
    integer:: N, M, i, j, lwork, info, LDg, LDU, LDVT
    real*8, allocatable, dimension(:) :: work, S
    real*8, allocatable, dimension(:,:) :: g,U,VT,S2

    N = 4
    M = 4
    lwork = 40
    jobu = 'N'
    jobvt = 'N'
    LDg = N
    LDU = N
    LDVT = N

    allocate(g(LDg,M))
    allocate(work(lwork))
    allocate(S(M))
    allocate(S2(N,M))
    allocate(U(LDU,N))
    allocate(VT(LDVT,M))

    g = 1

    call sgesvd(jobu,jobvt,N,M,g,LDg,S,U,LDU,VT,LDVT,work,lwork,info)

    print *, S

    deallocate(g)
    deallocate(work)
    deallocate(S)
    deallocate(S2)
    deallocate(U)
    deallocate(VT)

    contains

    real function Hij(N,Nh,locs,i,j) result(storprod)
        implicit none 
        integer, intent(in) :: N, i, j
        logical, dimension (:), intent(in) :: locs
        integer, intent(in) :: Nh
        integer :: k = 0, d1 = 0, d2 = 0

        storprod = 1
        do k = 0,N-1
            d1 = mod(ishft(i-1, -k),2)
            d2 = mod(ishft(j-1, -k),2)
            if (locs(k+1)) then
                if (d1==1 .and. d2 ==1) then
                    storprod = -storprod
                end if
            end if
            if (.not. locs(k+1)) then
                if (.not. d1 == d2) then
                    storprod = 0
                end if
            end if
        end do

        storprod = storprod*1/sqrt(2.0)**Nh
    end function Hij

    function InvPerm(perm) result(inv_perm)
        implicit none
        integer, dimension (:), allocatable :: inv_perm
        integer, dimension (:), intent(in) :: perm
        integer :: i = 0
        allocate(inv_perm(size(perm)))
        do i = 1,size(perm)
            inv_perm(perm(i)) = i
        end do
    end function InvPerm

    function get_mat(N, Nh, locs, P, Q) result(C)
        implicit none
        logical, dimension(:), intent(in) :: locs
        integer, intent(in) :: N, Nh
        integer, dimension(:), intent(in) :: P, Q
        real, dimension(:,:), allocatable :: C
        integer :: i = 0, j = 0, k = 0, hadim, l=0
        real :: storsum = 0
        integer, dimension(N) :: Qinv, Pinv
        allocate(C(2**N, 2**N))
        Pinv = InvPerm(P)
        Qinv = InvPerm(Q)
        hadim = int(N/2)
        do i=1,2**N
            do j = 1,2**N
                storsum = 0
                do k = 1,2**N
                    l = Qinv(k)
                    l = l + (-1)**mod(ishft(l-1, -hadim),2)*2**hadim
                    storsum = storsum + Hij(N,Nh,locs,Pinv(i),k)*Hij(N,Nh,locs,Q(l),Pinv(j))
                end do 
                C(i, j) = storsum
            end do
        end do
    end function get_mat

    subroutine ptranspose(mat, hadim)
        implicit none
        integer :: in1 = 0, in2 = 0, out1 = 0, out2 = 0
        real, dimension(:, :), intent(inout) :: mat
        integer, intent(in) :: hadim
        real :: tmp = 0
        do i = 0,hadim**2-1
            do j = 0,hadim**2-1
                in2 = mod(i, hadim)
                in1 = i/hadim
                out2 = mod(j, hadim)
                out1 = j/hadim
                if (in2 < out1) then
                    tmp = mat(i+1, j+1)
                    mat(i+1,j+1) = mat(in1*hadim+out1+1, in2*hadim+out2+1)
                    mat(in1*hadim+out1+1, in2*hadim+out2+1) = tmp
                end if
            end do
        end do
    end subroutine ptranspose

end program printH