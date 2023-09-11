
function sum_of_ints(a,n) result(s)
    integer, intent(in) :: a(*), n
    integer :: s, i

    s = 0
    do i=1,n
        s = s + a(i)
    end do
    print*, "sum of ", n, " elements is ", s
end function sum_of_ints


function sum_of_doubles(a,n) result(s)
    double precision, intent(in) :: a(*)
    integer, intent(in) :: n
    double precision :: s
    integer :: i

    s = 0
    do i=1,n
        s = s + a(i)
    end do
    print*, "sum of ", n, " elements is ", s
end function sum_of_doubles
