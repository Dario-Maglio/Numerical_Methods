Module library


  Implicit None
  Integer(4), Parameter::rk=8
  Real(rk), Parameter::pi=Acos(-1._rk)

 

 Abstract Interface
     Function fun_R_R(x)
       Import :: rk
       Real(rk), Intent(in) :: x
       Real(rk) :: fun_R_R
     End Function fun_R_R

     Function fun_vR_vRR(y,x)
       Import::rk
       Real(rk),Intent(In)::y(1:),x
       Real(rk)::fun_vR_vRR(1:Size(y))
     End Function fun_vR_vRR

     Function fun_mR_vRR(y,x)
       Import::rk
       Real(rk),Intent(In)::y(1:),x
       Real(rk)::fun_mR_vRR(1:Size(y),1:Size(y))
       End Function fun_mR_vRR
    End Interface

    Type::ode_data
       Real(rk), Allocatable::y0(:),y1(:),y(:),u(:,:)
       Real(rk)::x0,x1,dx
       Procedure(fun_vR_vRR), Nopass, Pointer::fun
       Procedure(fun_vR_vRR), Nopass, Pointer::dfun_t
       Procedure(fun_mR_vRR), Nopass, Pointer::dfun_y
    End Type ode_data

    Type::ode_para
       Real(rk)::eps
       Integer(4)::steps,max_steps
       Character(len=1)::Case
    End Type ode_para

    Type :: matrix
     Real(rk), Allocatable :: A(:,:),At(:,:),iA(:,:),b(:),bt(:),x(:),Bm(:,:),Bmt(:,:)
     Integer(4), Allocatable :: v(:)
     Integer(4) :: p
     Real(rk) :: detA,lambda
  End Type matrix

  Type :: matrix_para
     Real(rk) :: eps_x,eps_lambda
     Integer(4) :: steps,max_steps
  End Type matrix_para

Contains

  
  
  Subroutine interp_lagrange(x,y,x0,y0)

    Implicit None
    Real(rk), Intent(In)::x(:), y(:), x0
    Real(rk), Intent(Out)::y0

    Integer(4)::l,m,N
    Real(rk)::prod
    
    N=size(x)
    if(Size(x)/=Size(y))Then
       Print *,'size x is different from size y!'
       stop
    End if
    Do l=1,N
       if(x(l)==x0) Then
          y0=y(l)
          return
       End if
    end do

    y0=0._rk
    Do l=1,N
       prod=1._rk
       Do m=1,N
          if(m==l) Cycle
          prod=prod*(x0-x(m))/(x(l)-x(m))
       End do
       y0=y0+prod*y(l)
    End do
  End Subroutine interp_lagrange


 Subroutine interp_aitken(x,y,x0,y0)
   

    Implicit None
    Real(rk), Intent(In) :: x(:), y(:), x0
    Real(rk), Intent(Out) :: y0

    Integer(4) :: k, i, N
    Real(Rk) :: y_(1:Size(x))

    N=Size(x)
    If (Size(x)/=Size(y)) Then
       Print *,'Size x is different from Size y!'
       Stop
    End If
    Do k=1,N
       If (x(k)==x0) Then
          y0=y(k)
          Return 
       End If
    End Do

    y_=y
    Do k=1,N-1
       Do i=1,N-k
          y_(i)=(x0-x(i+k))/((x(i)-x(i+k)))*y_(i)+(x0-x(i))/((x(i+k)-x(i)))*y_(i+1)
       End Do
    End Do
    y0=y_(1)
  End Subroutine interp_aitken



  
  Function diff(x,y,l,o,m,e,dy) Result(stat)
    Implicit None
    Real(rk), Intent(In)::x(0:),y(0:)
    Integer(4), Intent(In)::l,o,m,e
    Real(rk), Intent(Out)::dy
    Integer(4)::stat

    Integer(4)::k,N
    Real(rk)::xm,xp,ym,yp,bn(1:3,1:3)=0

    bn(1,1)=1._rk
    bn(1,2)=4._rk/3._rk
    bn(2,2)=-1._rk/3._rk
    bn(1,3)=1.5_rk
    bn(2,3)=-0.6_rk
    bn(3,3)=0.1_rk

    N=Size(x)-1
    If(Size(x)/=Size(y)) Then
       stat=1
       Return
    End If
    If(l<0.or.l>N) Then
       stat=2
       Return
    End If
    If(o<1.or.o>2) Then
       stat=3
       Return
    End If
    If(m>3) Then
       stat=4
       Return
    End If

    dy=0
    Do k=1,m
       If(l<k) Then
          xm=x(0)+(l-k)*(x(1)-x(0))
          xp=x(l+k)
          Call interp_aitken(x(0:Min(e-1,N)),y(0:Min(e-1,N)),xm,ym)
          yp=y(l+k)
       Else If(l>N-k) Then
          xm=x(l-k)
          xp=x(N)+(l-(N-k))*(x(N)-x(N-1))
          ym=y(l-k)
          Call interp_aitken(x(Max(0,N+1-e):N),y(Max(0,N+1-e):N),xp,yp)
       Else
          xm=x(l-k)
          xp=x(l+k)
          ym=y(l-k)
          yp=y(l+k)
       End If
       Select Case(o)
       Case(1)
          dy=dy+bn(k,m)*(yp-ym)/(xp-xm)
       Case(2)
          dy=dy+bn(k,m)*(yp-2._rk*y(l)+ym)/(0.25_rk*(xp-xm)**2)
       End Select
    End Do
    stat=0
  End Function diff
  
  Recursive function diff_fun(fun,x0,o,h_in,eps,max_steps,steps_,Dfun) Result(stat)

    Implicit None
    Procedure(fun_R_R):: fun
    Real(rk), intent(in)::x0,h_in,eps
    Integer(4), intent(in)::o,max_steps
    Integer(4), intent(out)::steps_
    Real(rk), intent (out):: Dfun
    Integer(4):: stat

    integer(4), Save::steps=1
    Real(rk):: dyh, dyh2
    Real(rk), Save::dyh_old

    If(steps==1) then
       Select case(o)
       case(1)
          dyh=(fun(x0+h_in)-fun(x0-h_in))/(2._rk*h_in)
       case(2)
          dyh=(fun(x0+h_in)-2._rk*fun(x0)+fun(x0-h_in))/(h_in**2)
       case default
          stat=2
          steps_=0
          steps=1
          Return
       End Select
    Else
       dyh=dyh_old
    End If
    Select Case(o)
    Case(1)
       dyh2=(fun(x0+h_in/2)-fun(x0-h_in/2))/h_in
    Case(2)
       dyh2=(fun(x0+h_in/2)-2._rk*fun(x0)+fun(x0-h_in/2))/((h_in**2)/4._rk)
    End Select
    dyh_old=dyh2
    If(Abs(dyh-dyh2)<=eps*Abs(dyh2)) Then
       Dfun=4._rk/3._rk*dyh2-dyh/3._rk
       stat=0
       steps_=steps
       steps=1
       Return
    End If
    steps=steps+1
    If(steps>max_steps) Then
       stat=1
       steps_=max_steps
       steps=1
       return
    End If
    stat=diff_fun(fun,x0,o,h_in/2,eps,max_steps,steps_,Dfun)
    Print *, h_in/2, Dfun, cos(x0),Abs(dfun-cos(x0))/abs(cos(x0)),eps
     If(stat/=0) Then
       steps_=steps
       steps=1
       End if
    End Function diff_fun
    
    Function integrate(x,y,m,integrale) Result(stat)

      Implicit None
      Real(rk), Intent(in)::x(0:), y(0:)
      Integer(4), Intent(in)::m
      Real(rk), Intent(out)::integrale
      Integer(4)::stat


      Real(rk)::w(0:Size(x)-1), h
      Integer(4):: i, N

      N=Size(x)-1
      Integrale=0
      If(Size(x)/=Size(y)) Then
         stat=1
         Return
      End If
         If(m<1.or.m>2) Then
            stat=2
            Return
         End If
      Select Case(m)
      Case(1)
         w=1._rk
         w(0)=0.5_rk
         w(N)=0.5_rk
      Case(2)
         w(0)=1._rk/3._rk
         w(N)=1._rk/3._rk
         Do i=1,N-1    
           If((i/2)*2==i) Then
            w(i)=2._rk/3._rk
           else
            w(i)=4._rk/3._rk
           end if
         End Do
         End Select
            h=x(1)-x(0)
            
            integrale=h*Sum(w*y)
        
   stat=0
   
        
            
 End Function Integrate


 Recursive Function int_fun_trap(fun,a,b,fa,fb,eps,min_steps,max_steps,steps_,Ifun) Result(stat)
   
    Implicit None
    Procedure(fun_R_R) :: fun
    Real(rk), Intent(In) :: a,b,fa,fb,eps
    Integer(4), Intent(In) :: min_steps,max_steps
    Integer(4), Intent(Out) :: steps_
    Real(rk), Intent(Out) :: Ifun
    Integer(4) :: stat

    Integer(4), Save :: steps=1
    Real(rk), Save :: a0,b0
    Real(rk) :: h,c,fc,S_I,S_II,Ifun_I,Ifun_II

    If (steps==1) Then
       a0=a
       b0=b
    Endif
    h=b-a
    c=0.5_rk*(a+b)
    S_I=0.5_rk*h*(fa+fb)
    fc=fun(c)
    S_II=0.5_rk*S_I+0.5_rk*h*fc
    If ((Abs(S_II-S_I)<=3._rk*eps*Abs(S_II)).And.(steps>=min_steps)) Then
       Ifun=S_II
       stat=0
    Else
       steps=steps+2
       If (steps>max_steps) Then
          Ifun=0._rk
          stat=1
       Else
          stat=Min(1,int_fun_trap(fun,a,c,fa,fc,eps,min_steps,max_steps,steps_,Ifun_I)&
               +int_fun_trap(fun,c,b,fc,fb,eps,min_steps,max_steps,steps_,Ifun_II))
          If (stat==0) Ifun=Ifun_I+Ifun_II
       End If
    Endif
    If (a==a0.And.b==b0) Then 
       steps_=steps
       steps=1
    Endif
  End Function int_fun_trap
  
  Recursive Function bisection(fun,A,B,epsx,epsy,max_s,zero,steps) Result (stat)

    Implicit None
    Procedure(fun_R_R)::fun
    Real(rk), Intent(in)::A,B,epsx,epsy
    Integer(4),Intent(in)::max_s
    Integer(4), Intent(Out)::steps
    Real(rk),Intent(Out)::zero
    Integer(4)::stat

    Real(rk)::A_,B_,C_

   steps=1
    
    If(fun(A)*fun(B)<0) Then
       A_=A
       B_=B
       C_=0.5_rk*(A_+B_)
       
       DO While(ABS(fun(C_))>epsy.AND.B_-A_>epsx)
          If(fun(A_)*fun(C_)<0) Then
             B_=C_
          else
             A_=C_
          End If
          C_=0.5_rk*(A_+B_)
          steps=steps+1

          If(steps>max_s) Then
             stat=2
             Return
          End If
       End Do
    Else
       stat=1
       Return
    End If

    stat=0
    zero=C_
    End Function bisection

    Recursive Function Newton(fun,dfun,A,epsy,max_s,zero,steps) Result(stat)

       Implicit None
       Procedure(fun_R_R)::fun, dfun
       Real(rk), Intent(in)::A,epsy
       Integer(4),Intent(in)::max_s
       Integer(4), Intent(Out)::steps
       Real(rk),Intent(Out)::zero
       Integer(4)::stat

       Real(rk)::X_
       X_=A
       steps=1

       Do While(ABS(fun(X_))>epsy)
          X_=X_-(fun(X_)/dfun(X_))
          steps=steps+1
       End Do

        If(steps>max_s) Then
             stat=2
             Return
          End If
       Stat=0
       zero=X_

     End Function Newton

     Recursive Function Secante(fun,A,B,epsy,max_s,zero,steps) Result(stat)
       
       Implicit None
       Procedure(fun_R_R)::fun
       Real(rk), Intent(in)::A,B,epsy
       Integer(4),Intent(in)::max_s
       Integer(4), Intent(Out)::steps
       Real(rk),Intent(Out)::zero
       Integer(4)::stat

       Real(rk):: x1,x2,x3
       x1=A
       x2=B
       steps=1

       Do While(ABS(fun(x2))>epsy)
          
          x3=x2-fun(x2)/((fun(x2)-fun(x1))/(x2-x1))
          x1=x2
          x2=x3
          steps=steps+1
       End Do
        If(steps>max_s) Then
             stat=2
             Return
          End If
          stat=0
          zero=x2
     End Function Secante
        
     Function euler(data) Result(stat)
       Implicit None
       Type(ode_data)::data
       Integer(4)::stat

       If(Size(data%y0)/=Size(data%y)) Then
          stat=1
          Return
        End If
    
      data%y=data%y0+data%dx*data%fun(data%y0,data%x0)
      stat=0
       End Function euler

       Function Picard(data) Result(stat)

          Implicit None
       Type(ode_data)::data
       Integer(4)::stat

       If(Size(data%y0)/=Size(data%y)) Then
          stat=1
          Return
        End If

       data%y=data%y0+0.5_rk*data%dx*(data%fun(data%y0,data%x0)+data%fun(data%y,data%x0+data%dx))
        stat=0
      End Function Picard

      Function Taylor(data) Result(stat)

          Implicit None
       Type(ode_data)::data
       Integer(4)::stat

       If(Size(data%y0)/=Size(data%y)) Then
          stat=1
          Return
        End If
    
        data%y=data%y0+data%dx*data%fun(data%y0,data%x0)+0.5_rk*(data%dx**2)*&
             (data%dfun_t(data%y0,data%x0)+Matmul(data%dfun_y(data%y0,data%x0),data%fun(data%y0,data%x0)))
        stat=0
      End Function Taylor

                                                        
      Function Picard_iterativo(data,para) Result(stat)

          Implicit None
          Type(ode_data)::data
          Type(ode_para)::para
          Real(rk)::yn(1:Size(data%y0)),yo(1:Size(data%y0))
       Integer(4)::stat

       If(Size(data%y0)/=Size(data%y)) Then
          stat=1
          Return
       End If

       yn=data%y
       para%steps=1

       Do While(ANY(ABS(yn-yo)>para%eps).OR.para%steps==1)
          yo=yn
          yn=data%y0+0.5_rk*data%dx*(data%fun(data%y0,data%x0)+data%fun(yo,data%x0+data%dx))
          para%steps=para%steps+1
          If(para%steps>para%max_steps)Then
             stat=2
             Return
          End If
       End Do
       
       data%y=yn
       stat=0
     End Function Picard_iterativo

     Function D3(Data) Result (stat)
 
    Implicit None
    Class(ode_data) :: Data
    Integer(4) :: stat

    If (Size(Data%y0)/=Size(Data%y)) Then
       stat=1
       Return
    End If
    If (Size(Data%y1)/=Size(Data%y)) Then
       stat=1
       Return
    End If
    Data%y=Data%y0+2._rk*Data%dx*Data%fun(Data%y1,Data%x0+Data%dx)
    stat=0
  End Function D3


  Function I3(Data) Result (stat)

    Implicit None
    Class(ode_data) :: Data
    Integer(4) :: stat

    If (Size(Data%y0)/=Size(Data%y)) Then
       stat=1
       Return
    End If
    If (Size(Data%y1)/=Size(Data%y)) Then
       stat=1
       Return
    End If
    Data%y=Data%y0+Data%dx/3._rk*(Data%fun(Data%y,Data%x0+2._rk*Data%dx)&
         +4._rk*Data%fun(Data%y1,Data%x0+Data%dx)+Data%fun(Data%y0,Data%x0))
    stat=0
  End Function I3


    Function Runge_Kutta(data) Result(stat)

    Implicit None
    Type(ode_data)::data
    Real(rk)::c1(1:Size(data%y)),c2(1:Size(data%y))
    Real(rk)::c3(1:Size(data%y)),c4(1:Size(data%y))
    Integer(4)::stat

    If(Size(data%y)/=Size(data%y0)) Then
       stat=1
       Return
    End If

    c1=data%dx*data%fun(data%y0,data%x0)
    c2=data%dx*data%fun(data%y0+0.5_rk*c1,data%x0+0.5_rk*data%dx)
    c3=data%dx*data%fun(data%y0+0.5_rk*c2,data%x0+0.5_rk*data%dx)
    c4=data%dx*data%fun(data%y0+c3,data%x0+data%dx)

    data%y=data%y0+(1._rk/6._rk)*(c1+2._rk*c2+2._rk*c3+c4)

    stat=0
    
  End Function Runge_Kutta

  Function gauss_jordan(data,N) Result(stat)
    Implicit None
    Type(matrix)::data
    Integer(4),Intent(In)::N
    Integer(4)::i,j,k
    Real(rk)::u_(N,N)
    Integer(4)::stat

     
     u_=data%A
    Do k=1,N-1
       data%At=u_
       Do i=k+1,N

            Do j=k,N

             u_(i,j)=data%At(i,j)-((data%At(i,k)/data%At(k,k))*data%At(k,j))

          End Do
       End Do
    End do
    data%At=u_
    
    stat=0
  End Function gauss_jordan
  
  Function equation_system(data) Result(stat)

    Implicit None
    Type(matrix)::data
    Integer(4)::i,j,k
    Real(rk)::u_(Size(data%b),Size(data%b)),s,b(Size(data%b))
    Integer(4)::stat

     If(Size(data%x)/=Size(data%b)) Then
       stat=1
       Return
    End If

    
    
   s=0._rk
   u_=data%A
   b=data%b
    Do k=1,Size(data%b)-1
       data%At=u_
       data%bt=b
       Do i=k+1,Size(data%b)

          b(i)=data%bt(i)-((data%At(i,k)/data%At(k,k))*data%bt(k))

            Do j=k,Size(data%b)

             u_(i,j)=data%At(i,j)-((data%At(i,k)/data%At(k,k))*data%At(k,j))

          End Do
          
       End Do
       
      
    End do
    data%At=u_
    data%bt=b
    
    
    data%x(Size(data%b))=data%bt(Size(data%b))/data%At(Size(data%b),Size(data%b))

    Do i=Size(data%b)-1,1,-1
       Do j=i+1,Size(data%b)

          s=s+(data%At(i,j)*data%x(j))
          
       End Do
       
      data%x(i)=1._rk/data%At(i,i)*(data%bt(i)-s)
      
    End Do
   
    stat=0
  End Function equation_system

  Function determinante(data,N) Result(stat)

    Implicit None
    Type(matrix)::data
    Integer(4)::i,j,k
    Integer(4), Intent(In)::N
    Real(rk)::u_(N,N)
    Integer(4)::stat


     data%detA=1._rk

     If(gauss_jordan(data,N)/=0) Then
        stat=1
        Return
     End If
     
       Do i=1,N
          data%detA=data%detA*data%A(i,i)
       End Do
      
       stat=0
     End Function determinante

     Function Inversa(A,N,inv) Result(stat)

    Implicit None
    Type(matrix)::data
    Integer(4),Intent(In)::N
    Real(rk),Intent(In)::A(:,:)
    Real(rk)::iden(N,N),c(N),b(N)
    Real(rk),Intent(Out)::inv(:,:)
    Integer(4)::i,j,k
    Integer(4)::stat

    Allocate(data%A(N,N),data%x(N),data%b(N))

    Do i=1,N
       Do j=1,N

          If(i==j)Then
             iden(i,j)=1._rk
          else
             iden(i,j)=0._rk
          End If
       End Do
    End Do

    Do j=1,N

       Do i=1,N

          b(i)=iden(i,j)
       End do

       data%A=A
       data%b=b
      
       If(equation_system(data)/=0) Then
          stat=1
          Return
       End If
        c=data%x
       inv(:,j) = c        
    End Do
    stat=0
  End Function Inversa

   Function tridiagsup(Data) Result (stat)

    Implicit None
    Type(matrix) :: Data
    Integer(4) :: stat

    Integer(4) :: n,i,j,k,q,r,bp,bm
    Real(rk) :: pivot

    n=Size(Data%A,1)

    If (n/=Size(Data%A,2)) Then
       stat=1
       Return
    End If

    If (Allocated(Data%v)) Deallocate(Data%v)
    Allocate(Data%v(1:n))
    Forall(i=1:n) Data%v(i)=i

    bp=0
    If (Allocated(Data%b)) Then
       If (Size(Data%b)/=n) Then
          stat=4
          Return
       End If
       bp=1
       If (Allocated(Data%bt)) Deallocate(Data%bt)
       Allocate(Data%bt(1:n))
       Data%bt=Data%b       
    Endif

    bm=0
    If (Allocated(Data%Bm)) Then
       If (Size(Data%Bm,1)/=n.Or.Size(Data%Bm,2)/=n) Then
          stat=5
          Return
       End If
       bm=1
       If (Allocated(Data%Bmt)) Deallocate(Data%Bmt)
       Allocate(Data%Bmt(1:n,1:n))
       Data%Bmt=Data%Bm
    Endif

    If (Allocated(Data%At)) Deallocate(Data%At)
    Allocate(Data%At(1:n,1:n))
    Data%At=Data%A

    Data%p=1
    Do k=1,n-1
       If (Any(Maxval(Abs(Data%At(Data%v(k:),:)),2)==0.D0)) Then
          stat=2
          Return
       Endif
       If (Maxval(Abs(Data%At(Data%v(k:),k)))==0.D0) Then
          stat=3
          Return
       Endif
       q=Maxloc(Abs(Data%At(Data%v(k:),k))/Maxval(Abs(Data%At(Data%v(k:),:)),2),1)+(k-1)
       If (q/=k) Then
          r=Data%v(q)
          Data%v(q)=Data%v(k)
          Data%v(k)=r
          Data%p=-Data%p
       Endif
       Do i=k+1,n
          pivot=Data%At(Data%v(i),k)
          Do j=1,n
             If (bm==1) Data%Bmt(Data%v(i),j)=Data%Bmt(Data%v(i),j)-pivot/Data%At(Data%v(k),k)*Data%Bmt(Data%v(k),j)
             If (j<k) Cycle
             Data%At(Data%v(i),j)=Data%At(Data%v(i),j)-pivot/Data%At(Data%v(k),k)*Data%At(Data%v(k),j)
          Enddo
          If (bp==1) Data%bt(Data%v(i))=Data%bt(Data%v(i))-pivot/Data%At(Data%v(k),k)*Data%bt(Data%v(k))
       Enddo
    Enddo
    stat=0

  End Function tridiagsup

  
  Function inv(Data) Result (stat)

    Implicit None
    Type(matrix) :: Data
    Integer(4) :: stat

    Integer(4) :: n,i,j,k
    Real(rk) :: sum

    n=Size(Data%A,1)

    If (n/=Size(Data%A,2)) Then
       stat=1
       Return
    End If

    If (Allocated(Data%Bm)) Deallocate(Data%Bm)
    Allocate(Data%Bm(1:n,1:n))
    Data%Bm=0._rk
    Do i=1,n
       Data%Bm(i,i)=1._rk
    Enddo

    stat=tridiagsup(Data)
    If (stat/=0) Return    

    If (Allocated(Data%iA)) Deallocate(Data%iA)
    Allocate(Data%iA(1:n,1:n))

    Do j=1,n
       Data%iA(n,j)=Data%Bmt(Data%v(n),j)/Data%At(Data%v(n),n)
    Enddo
    Do i=n-1,1,-1
       Do j=1,n
          sum=0
          Do k=i+1,n
             sum=sum+Data%At(Data%v(i),k)*Data%iA(k,j)
          Enddo
          Data%iA(i,j)=(Data%Bmt(Data%v(i),j)-sum)/Data%At(Data%v(i),i)
       Enddo
    Enddo
    stat=0

  End Function inv

  Function det(Data) Result (stat)

    Implicit None
    Type(matrix) :: Data
    Integer(4) :: stat

    Integer(4) :: i,n

    n=Size(Data%A,1)

    If (n/=Size(Data%A,2)) Then
       stat=1
       Return
    End If

    stat=tridiagsup(Data)
    If (stat/=0) Return    

    Data%detA=Data%p
    Do i=1,n
       Data%detA=Data%detA*Data%At(Data%v(i),i)
    Enddo
    stat=0

  End Function det


  Function eigen(data,para) Result(stat)

    Implicit None
    Type(matrix)::data
    Type(matrix_para)::para
    Real(rk)::rand_x(1:Size(data%A,1)),l
    Integer(4)::stat

    Call Random_seed()
    Call Random_Number(rand_x)

    para%steps=1
    data%x=rand_x
    Do While(para%steps==1.or.Sqrt(Dot_Product(data%x-rand_x,data%x-rand_x))>=para%eps_x.or.ABS(l-data%lambda)>=para%eps_lambda)

       rand_x=data%x
       data%x=(Matmul(data%A,rand_x))/(Sqrt(Dot_Product(Matmul(data%A,rand_x),Matmul(data%A,rand_x))))
        
       if(para%steps>=2) Then 
          l=data%lambda
       End If
       
       data%lambda=(Dot_Product(Matmul(data%A,data%x),data%x))/Dot_Product(data%x,data%x)

       para%steps=para%steps+1

       If(para%steps>para%max_steps) Then
          stat=1
          Return
       End If
    End Do
    stat=0
  End Function eigen
                 
End Module library
     
  
