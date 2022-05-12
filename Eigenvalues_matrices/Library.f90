Module library

  Implicit None 
  Integer(4), Parameter :: rk=8
  Real(rk), Parameter :: pi=Acos(-1._rk)

  Abstract Interface                                        !Dichiarare come è fatta la funzione che verrà chiamata in una subroutine
     Function fun_R_R(x)                                    !Nome e variabili funzione
       Import :: rk                                         !Prende rk da sopra (necessario xk non è nel contains)
       Real(rk), Intent(in) :: x
       Real(rk) :: fun_R_R                                  !Anche nome della f è una variabile
     End Function fun_R_R
     Function fun_vR_RvR(y,x)
       import:: rk
       Real(rk), intent (in):: y(1:),x
       Real (rk)::fun_vR_RvR(1:size(y))
     end function fun_vR_RvR
     Function fun_mR_RvR (y,x)
       import:: rk
       Real(rk), intent (in):: y(1:),x
       Real (rk)::fun_mR_RvR(1:size(y),1:size(y))
     end function fun_mR_RvR
  End Interface

  Type::ode_data
     Real(rk), allocatable:: y0(:),y1(:),y(:),u(:,:)
     Real(rk)::x0,x1,dx
     Procedure(fun_vR_RvR), Nopass, pointer::fun            !Non passa alla funzione il tipo ode data come argomento
     Procedure(fun_vR_RvR), Nopass, pointer::dfun_t
     Procedure(fun_mR_RvR), Nopass, pointer::dfun_y  
  End Type ode_data

  Type::ode_para
     Real(rk)::eps
     integer(4)::steps,max_steps
     character(1)::Case
  end type ode_para

  Type :: matrix
     Real(rk), Allocatable :: A(:,:),At(:,:),iA(:,:),b(:),bt(:),x(:),Bm(:,:),Bmt(:,:)
     Integer(4), Allocatable :: v(:)
     Integer(4) :: p
     Real(rk) :: detA,TrA,lambda
  End Type matrix

  Type :: matrix_para
     Real(rk) :: eps_x,eps_lambda
     Integer(4) :: steps,max_steps
  End Type matrix_para
 
End Module library 
