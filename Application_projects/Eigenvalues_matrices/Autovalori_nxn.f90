Program drive_relazione
  !Per cambiare le dim di A basta cambiare n e il formato di stampa per ogni print
  
  Use LibraryMatr

  Implicit None
  integer(4), parameter:: n=6
  integer(4):: i
  Type(matrix) :: Data
  Type(matrix_para) :: para
  Real(rk) :: lambda(1:n),u(3:n),Uni(1:n,1:n),A(1:n,1:n),z,y
  Character(Len=1) :: n_l

  Write(n_l,'(I1)') n                                    !A cosa serve questo comando?

  para%eps_x=1.E-6_rk
  para%eps_lambda=1.E-6_rk
  para%max_steps=99

  A=Reshape([ 1._rk,1._rk,2._rk,5._rk,8._rk,0._rk,&
              1._rk,-39._rk,1._rk,0._rk,0._rk,2._rk,&
              2._rk,1._rk,-25._rk,2._rk,3._rk,0._rk,&
              5._rk,0._rk,2._rk,30._rk,0._rk,4._rk,&
              8._rk,0._rk,3._rk,0._rk,110._rk,3._rk,&
              0._rk,2._rk,0._rk,4._rk,3._rk,190._rk ],[n,n])
  
  Uni=0._rk
  Forall (i=1:n) Uni(i,i)=1._rk  
  Allocate(data%A(1:n,1:n))
  data%A=A
  
  Print *,''
  Print *,'            A'
  Do i=1,n
     Print '(6(x,E13.6))',data%A(i,:)
  End do
  Print *,'-----------------------------------'
  Print *,'-----------------------------------'
  Print *,''

  Print *,'Info matrice A'
  Print '(x,A,I1)','stat det = ',det(Data)
  Print '(x,A,x,E13.6)','det A= ',data%detA
  Data%TrA=0._rk
  Do i=1,n
     data%TrA=data%TrA + A(i,i)  
  End do
  Print '(x,A,x,E13.6)','Trac A= ',data%TrA
  Print '(x,A,I1)','stat inv = ',inv(Data)
  Print *,'             inv A'
  Do i=1,n
     Print '(6(x,E13.6))',data%iA(i,:)
  End do
  Print *,''
  Print *,'check A * inv A = I'
  Do i=1,n
     Print '(6(x,E13.6))',Matmul(data%A(i,:),data%iA)
  End do
  Print *,'-----------------------------------'
  Print *,'-----------------------------------'
  Print *, ''

  !Funzione che restituisce u e calcola centri e raggi di Gershgorin
  Print '(I4,x,A)',Gershgorin(Data,u),'-> stat Gershgorin'
  print *,''
  !u=[3._rk,-3._rk,11._rk,7._rk] !Vett coeff vicini agli autovalori (commentare Gershgorin per usare)

  Print *, 'Inizio ricerca di autovalori e autovettori:'

  Print *,'   1° autovalore (maggiore)'
  Print '(x,A,I1)','   stat eigen = ',eigen(Data,para)
  lambda(1)=data%lambda
  Print '(x,A,6(x,E13.6),x,A,(x,E13.6))','   1° x= (', data%x, ')       1° lambda=', lambda(1)
  Print '(x,A,x,I4)','   steps= ',para%steps
  Print *,''

  Print *, '   2° autovalore (minore)'
  data%A=data%iA
  Print '(x,A,I1)','   stat eigen = ',eigen(Data,para)
  lambda(2)=1._rk/data%lambda
  Print '(x,A,6(x,E13.6),x,A,(x,E13.6))','   2° x= (', data%x, ')       2° lambda=', lambda(2)
  Print '(x,A,x,I4)','   steps= ',para%steps
  Print *,''
  
  Do i=3,n-2
     Print '(x,I4,x,A)',i,'° autovalore (cerchi)'
     data%A=A-u(i)*Uni
     Print '(x,A,I1)','   stat inv = ',inv(Data)
     data%A=data%iA
     Print '(x,A,I1)','   stat eigen = ',eigen(Data,para)
     lambda(i)=(1._rk/data%lambda)+u(i)
     Print '(x,I4,x,A,6(x,E13.6),x,A,x,I4,x,A,(x,E13.6))',i,'° x = (', data%x, ') ',i,'° lambda =', lambda(i)
     Print '(x,A,x,I4)','   steps= ',para%steps
     Print *,''
  End do

  !Ripristino det(A)
  data%A=A
  y=det(data)
  
  Print *, 'Ultimi 2 autovalori (traccia e determinante)'
  print *,''
  z=data%TrA
  Do i=1,n-2
     z = z - lambda(i)
  end do
  y=data%detA
  Do i=1,n-2
     y=y/lambda(i)
  end do
  !lambda(n-1) + lambda (n) = z
  !lambda(n-1) * lambda (n) = y
  lambda(n)=(z+sqrt(z*z - 4._rk*y))/2._rk
  lambda(n-1)= z - lambda(n)
  !Calcolo autovettori
  !Penultimo autovettore
  u(n-1)=lambda(n-1)-1.E-2_rk
  data%A=A-u(n-1)*Uni
  Print '(x,A,I1)','   stat inv penultimo= ',inv(Data)
  data%A=data%iA
  Print '(x,A,I1)','   stat eigen penultimo= ',eigen(Data,para)
  If (abs((1._rk/data%lambda)+u(n-1)-lambda(n-1))>1.E1_rk*para%eps_lambda) Stop 'La ricerca del penultimo autovettore è fallita'
  Print '(x,I4,x,A,6(x,E13.6),x,A,x,I4,x,A,(x,E13.6))',n-1,'° x = (', data%x, ') ',&
       n-1,'° lambda =', lambda(n-1)
  Print '(x,A,x,I4)','   steps= ', para%steps
  Print *,''
  !Ultimo autovettore
  u(n)=lambda(n)-1.E-2_rk
  data%A=A-u(n)*Uni
  Print '(x,A,I1)','   stat inv ultimo= ',inv(Data)
  data%A=data%iA
  Print '(x,A,I1)','   stat eigen ultimo= ',eigen(Data,para)
  If (abs(((1._rk/data%lambda)+u(n))-lambda(n))>10._rk*para%eps_lambda) Stop 'La ricerca dell ultimo autovettore è fallita'
  Print '(x,I4,x,A,6(x,E13.6),x,A,x,I4,x,A,(x,E13.6))',n,'° x = (', data%x, ') ',&
       n,'° lambda =', lambda(n)
  Print '(x,A,x,I4)','   steps= ', para%steps
  Print *,''
  
  

End Program drive_relazione
