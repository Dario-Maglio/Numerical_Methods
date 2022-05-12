program manyfiles

  implicit none
  character(len=20) :: filename="fdata", formato='(A5,I0,".dat")' 
  integer(4), parameter :: outunit=42,  N=5

  integer(4):: i

  do i=1, N
    write(filename, formato) filename, i
    open(unit=outunit,file=filename, form='formatted')
    
    
    write(outunit, *) i, N
    close(outunit)
  enddo
  
end program manyfiles
