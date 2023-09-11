

subroutine hello(name)
    character(len=*), intent(in) :: name 
    print* , "Hello, ", name, "!"
end subroutine hello

