program printH
    implicit none
    real, dimension (4,4) :: H
    integer :: i = 0, j = 0
    do i = 1,4
        do j = 1,4
            write(*, fmt = "(F8.1)", advance = "no") Hij(2, 1, [.false., .true.], j, i)
        end do
        print *, ""
    end do

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
        integer, dimension (:), allocatable :: inv_perm
        integer, dimension (:), intent(in) :: perm
        integer :: i = 0
        allocate(inv_perm(size(perm)))
        do i = 1,size(perm)
            inv_perm(perm(i)) = i
        end do
    end function InvPerm

    function get_mat(N, Nh, locs, P, Q) result(C)
        logical, dimension(:), intent(in) :: locs
        integer, intent(in) :: N, Nh
        integer, dimension(:), intent(in) :: P, Q
        real, dimension(:,:), allocatable :: C
        integer :: i = 0, j = 0, k = 0, hadim
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
                    if (mod(ishft(-hadim, Qinv(k)), 2) == 0) then
                        storsum = storsum + Hij(N, Nh, locs, Pinv(i), k) * Hij(N, Nh, locs, Pinv(j),Q(Qinv(k)+2**hadim))
                    else
                        storsum = storsum + Hij(N, Nh, locs, Pinv(i), k) * Hij(N, Nh, locs, Pinv(j),Q(Qinv(k)-2**hadim))
                    end if
                end do 
            end do
        end do
    end function get_mat

end program printH