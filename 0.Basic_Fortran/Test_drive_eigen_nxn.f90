Program drive_relazione
  
  Use library

  Implicit None
  integer(4), parameter:: n=4
  integer(4):: i
  Type(matrix) :: data
  Type(matrix_para) :: para
  Real(rk) :: lambda(1:n), vlamb(n,n)

  para%eps_x=1.E-6_rk
  para%eps_lambda=1.E-6_rk
  para%max_steps=99

  Allocate(data%A(1:n,1:n))
  data%A=Reshape([ 1._rk,5._rk,8._rk,5._rk,&
                   5._rk,-38._rk,2._rk,0._rk,&
                   8._rk,2._rk,50._rk,3._rk,&
                   5._rk,0._rk,3._rk,10._rk],[n,n])
  
  print *,''
  Print *,'            A'
  Do i=1,n
     print '(4(x,E13.6))',data%A(i,:)
  End do
  print *,'-----------------------------------'
  print *,''

  If (diagonalizationR(data, para, lambda, vlamb)/=0) then
    print*, 'Be careful with the Gershgorin circles.'
  end if
  
  Do i=1,n
     print *, ''
     print *, lambda(i)
     print '(4(x,E13.6))',vlamb(i,:)
     print *, ''
  End do
  

End Program drive_relazione
