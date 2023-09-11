
function sum_of_int(a,b) result(c)
    integer, intent(in) :: a, b
    integer :: c

    c = a + b
    print*,"sum of ", a, " and ", b, " is ", c
end function sum_of_int

function sum_of_double(a,b) result(c)
    double precision, intent(in) :: a, b
    double precision :: c

    c = a + b
    print*,"sum of ", a, " and ", b, " is ", c
end function sum_of_double
