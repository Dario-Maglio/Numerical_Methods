!!======================================================================================
!!== Module for numerical and data analysis in Fortran95. Author: Cafasso Dario.      ==
!!== Most of the subroutines from interpolation to ODE solution have been written     ==
!!== and optimized by professor Adolfo Avella of Salerno University, who shared       ==
!!== them with his students during the 'Fisica computazionale' course in UniSa.       ==
!!== Other contributions come from the Numerical methods course of Pisa University.   ==
!!== All of the contents have been strongly modified by the author, but also tested.  ==
!!======================================================================================

!!======================================================================================
!!Index
!!(the arguments followed by underscore are optional)
!!
!!General parameters and abstract interface
!!
!!Pseudorandom number generators:
!!-> ran2() (randomseed file is needed)
!!
!!Interpolation, derivation and integration algorithms:
!!-> interp(x,y,x0,y0,method_)
!!-> method: lagrange, aitken, updown
!!-> differentiate(x,y,l,dy,o_,m_,e_)
!!-> diff_fun(fun,x0,Dfun,o_,h_ini_,eps_,max_steps_,steps_)
!!-> integrate(x,y,integr,method_)
!!-> int_fun(fun,a,b,fa,fb,Ifun,eps_,min_steps_,max_steps_,steps_,fcc)
!!
!!Roots and ordinary differential equations:
!!-> bisection(fun,x,root,dx_,a_,b_,eps_fun_,eps_,max_steps_,steps_)
!!-> secant(fun,x,zero,dx_,eps_fun_,eps_,max_steps_,dxx_,a_,b_,steps_bis_,steps_)
!!-> newton(fun,dfun,x,zero,eps_fun_,eps_,max_steps_,dxx_,a_,b_,steps_bis_,steps_)
!!-> predictor(ode_data): Euler, Taylor3, D3
!!-> corrector(ode_data): I2 (or Picard), I3 (or Simpson)
!!-> corrector(ode_data, ode_para): IterPicard
!!-> Runge_Kutta(ode_data)
!!-> shooting(ode_data, ode_para)
!!
!!Linear algebra and eigenvalues problem:
!!-> tridiagsup(matrix)
!!-> linsys(matrix)
!!-> basic(matrix): det, trace, inv
!!-> eigen(matrix, matrix_para)
!!-> Gershgorin(matrix, out)
!!-> diagonalizationR(matrix, matrix_para, lambda, vlamb)
!!
!!======================================================================================




Module library

!!=================================================================
!!=================================================================
!!==	  General parameters and abstract interface            ====
!!=================================================================
!!=================================================================

  Implicit None 
  Integer(4), Parameter :: rk=8
  Real(rk), Parameter :: pi=Acos(-1._rk)
  
  Type::ode_data 
  !!Dati per sistema di equazioni differenziali della forma y = y0 + dx*fun(y0,x0)
  !!fun è una funzione vettoriale di un vettore e di un parametro reale (passo t)
  !!y1 e x1 sono un estremo finale mentre u conserva la traiettoria dello shooting
     Real(rk), allocatable:: y0(:),y1(:),y(:),u(:,:)
     Real(rk)::x0,x1,dx
     Procedure(fun_vR_RvR), Nopass, pointer::fun
     Procedure(fun_vR_RvR), Nopass, pointer::dfun_t
     Procedure(fun_mR_RvR), Nopass, pointer::dfun_y  
  End Type ode_data

  Type::ode_para
  !!Parametri di precisione delle equazioni differenziali
     Real(rk)::eps
     integer(4)::steps,max_steps
     character(1)::Case
  End type ode_para

  Type :: matrix
  !!Matrice, trasposta, inversa, coefficienti, matrici supporto, autovettore, autovalore ecc.
     Real(rk), Allocatable :: A(:,:),At(:,:),iA(:,:),b(:),bt(:),x(:),Bm(:,:),Bmt(:,:)
     Integer(4), Allocatable :: v(:)
     Integer(4) :: p
     Real(rk) :: detA,TrA,lambda
  End Type matrix

  Type :: matrix_para
  !!Parametri di precisione delle soluzioni di algebra lineare
     Real(rk) :: eps_x,eps_lambda
     Integer(4) :: steps,max_steps
  End Type matrix_para
  
!!*****************************************************************

  Abstract Interface                                        !Definire i tipi di funzioni che vengono utilizzati
     
     Function fun_R_R(x)
     !!Funzione reale di variabile reale
       Import :: rk                                         !Prende rk da sopra (necessario perché non è nel contains)
       Real(rk), Intent(in) :: x
       Real(rk) :: fun_R_R                                  !Anche nome della f è una variabile
     End Function fun_R_R
     
     Function fun_vR_RvR(y,x)
     !!Funzione vettoriale di variabile vettoriale y e reale x
       import:: rk
       Real(rk), intent (in):: y(1:),x
       Real (rk)::fun_vR_RvR(1:size(y))
     end function fun_vR_RvR
     
     Function fun_mR_RvR (y,x)
     !!Funzione matriciale di variabile vettoriale y e reale x
       import:: rk
       Real(rk), intent (in):: y(1:),x
       Real (rk)::fun_mR_RvR(1:size(y),1:size(y))
     end function fun_mR_RvR
  
  End Interface



Contains                                                    !Introduce procedure e funzioni della libreria o del sottoprogramma
!subroutine: funzione il cui return è un void               !subroutine si chiama, f si assegna
!funzione: restituisce valore con nome della funzione       !può restituire anche lo stato della funzione: 0 -> successo, altro -> errore

!!=================================================================
!!=================================================================
!!==             Pseudorandom number generators                ====
!!=================================================================
!!=================================================================

  Function ran2()
  !!Pseudorandom number generator with a long period (>2 x 10^18) before any aliasing occurs.
  !!The compiler generator period is 2^{256} - 1 (with multiple threads 2^{128}-1).
  
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,&
                ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,&
                ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,&
                rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
      !save iv,iy,idum2
      !data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
  End function ran2

!******************************************************************

  Subroutine ranstart
  !!Initializes the random numeber generator Ran2 reading the last saved state
  
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      read(23,*) idum
      read(23,*,end=117) idum2
      do i=1,32
         read(23,*) iv(i)
      enddo
      read(23,*) iy
      close(23)
      goto 118                          !!takes account of the first start
 117  if(idum.ge.0) idum = -idum -1
      close(23)
 118  continue

      return
  End subroutine ranstart

!******************************************************************
  
  Subroutine ranfinish
  !!Saves the last state of Ran2 generator
  
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      common /dasav/ idum,idum2,iv(32),iy

      open(unit=23, file='randomseed', status='unknown')
      write(23,*) idum
      write(23,*) idum2
      do i=1,32
         write(23,*) iv(i)
      enddo
      write(23,*) iy
      close(23)

      return
  End subroutine ranfinish

!!=================================================================
!!=================================================================
!!==   Interpolation, derivation and integration algorithms    ====
!!=================================================================
!!=================================================================

  Subroutine interp(x,y,x0,y0,method_)
  !!Compute the interpolated value y0 of the unknown function y=f(x) at x0
  !!given y(i) = f(x(i)) with i = 1...N by the selected method (aitken default).
  !!input: vettori punti e puntatori a x0, y0 in cui si desidera ottenere l'interpolazione.
  !!optional: method_ può essere lagrange, aitken o updown.
  !!note: consigliato l'utilizzo di 20-50 punti, dato che l'errore esplode con N.
    Implicit None
    Real(rk), Intent(In) :: x(:), y(:), x0                  !def come variabili in input con Intent(in) warning se modificate
    Character(6), optional :: method_
    Real(rk), Intent(Out) :: y0                             !qui analogamente warning se non viene assegnata y0 alla fine
    
    If (present(method_)) Then                              !definizione di argomento opzionale
    	Select Case(trim(method_))                          !Invece di un if annidato (è meglio così)
   	 Case('lagran')                                     !Trim taglia spazi vuoti in stringhe
   	   Call interp_lagrange(x,y,x0,y0)
   	 Case('aitken')
    	   Call interp_aitken(x,y,x0,y0)
         Case('updown')
   	   Call interp_updown(x,y,x0,y0)
   	 Case Default
    	   Print *,'Sorry, this method had not been implemented yet!'
    	End Select
    Else
    	Call interp_aitken(x,y,x0,y0)
    End If
  End Subroutine interp
  
!******************************************************************

  Subroutine interp_lagrange(x,y,x0,y0)
  !!Compute the interpolated value y0 of the unknown function y=f(x) at x0
  !!given y(i)=f(x(i)) with i = 1...N by means of the Lagrange method.
    Implicit None
    Real(rk), Intent(In) :: x(:), y(:), x0
    Real(rk), Intent(Out) :: y0

    Integer(4) :: l, m, N
    Real(rk) :: prod

    N=Size(x)                                      !Size prende la misura dell'array
    If (Size(x)/=Size(y)) Then
       Print *,'Size x is different from Size y!'
       Stop                                        !Ferma il programma
    End If
    Do l=1,N
       If (x(l)==x0) Then
          y0=y(l)
          Return                                   !Ferma l'esecuzione del sottoprogramma
       End If
    End Do

    y0=0._rk
    Do l=1,N
       prod=1._rk
       Do m=1,N
          If(m==l) Cycle                           !Salta il passo del ciclo 
          prod=prod*(x0-x(m))/(x(l)-x(m))
       End Do
       y0=y0+prod*y(l)
    End Do
  End Subroutine interp_lagrange

!******************************************************************

  Subroutine interp_aitken(x,y,x0,y0)
  !!Compute the interpolated value y0 of the unknown function y=f(x) at x0
  !!given y(i)=f(x(i)) with i = 1...N by means of the Aitken method.
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

!******************************************************************

  Function nearest_i(x,xs)
  !!Search the closest point x(i) to xs and return i, i = 1...N
  !!Adolfo Avella
  !!ver:1.0 Salerno, 2013.04.07
    
    Implicit None
    Integer(4) :: nearest_i
    Real(rk), Intent(In) :: x(:),xs

    nearest_i = Minloc(Abs(x-xs),1)

  End Function nearest_i

  Subroutine interp_updown(x,y,x0,y0)
  !!Compute the interpolated value y0 of the unknown function y=f(x) at x0
  !!given y(i)=f(x(i)) with i = 1...N by means of the plus-minus differences.
  !!Adolfo Avella
  !!ver:1.0 Salerno, 2013.04.07
  
    Implicit None
    Real(rk), Intent(In) :: x(:), y(:), x0
    Real(rk), Intent(Out) :: y0

    Real(rk) :: delta_plus(Size(x),Size(x))  ! Auxiliary matrices delta+
    Real(rk) :: delta_minus(Size(x),Size(x)) ! Auxiliary matrices delta-
    Real(rk) :: aux,dx,delta(Size(x))        ! Auxiliary quantity, difference between xs and closest point, delta to be summed up
    Integer(4) :: N,i,j,near,near_p,near_m   ! N, running integers, closest-point index, running indexes in re-summation

    delta_plus(:,1) = y  ! Loading the auxiliary matrices delta+(:,0) with values of y
    delta_minus(:,1) = y ! Loading the auxiliary matrices delta-(:,0) with values of y
    N = Size(x)          ! Compute the value of N
    Do j = 1, N-1        ! up-down recursion
       Do i = 1, N-j
          aux = (delta_plus(i,j) - delta_minus(i+1,j)) / (x(i) - x(i+j))
          delta_plus(i,j+1) = (x0 - x(i+j)) * aux
          delta_minus(i,j+1) = (x0 - x(i)) * aux
       Enddo
    Enddo
    near=nearest_i(x,x0)
    near_m=near
    near_p=near
    dx=x0-x(near)
    delta(1) = delta_plus(near,1)
    Do i= 2, N
       If ((dx<0.Or.near_p==N).And.near_m>1) Then
          near_m = near_m - 1
          delta(i) = delta_plus(near_m,i)
          dx = -dx
       Else
          delta(i) = delta_minus(near_m,i)
          near_p = near_p + 1
          dx = -dx
       Endif
    Enddo
    y0 = delta(N)
    Do i= N-1, 1, -1
       y0 = y0 + delta(i)
    Enddo
  End Subroutine interp_updown

!******************************************************************
!------------------------------------------------------------------
!******************************************************************

  Function differentiate(x,y,l,dy,o_,m_,e_) result (stat)          !result è var da utilizzare al posto di Der
  !!Compute the o-th derivative dy of the unknown function y=f(x) at x(l)
  !!given y(i)=f(x(i)) with i = 0...N by means of the (2m+1)-point formula.
  !!input: vettori punti, l indice di x0 in cui calcolare la derivata, dy puntatore al risultato.
  !!opzionali: o=1,2 ordine della derivata, m=1,2,3 per la derivata a 2m+1 punti, 
  !!e=10 punti da usare per estrapolazione ai bordi.
  !!note: l'errore è O(h^2m).
  !stat=0 successo
  !stat=1 Size(x)/=Size(y)
  !stat=2 l<0 or l>N
  !stat=3 o<1 or o>2                                
  !stat=4 m<1 or m>3

    Implicit none
    real(rk), Intent (In):: x(0:), y(0:)
    integer(4), Intent (In) :: l
    real(rk), Intent (Out):: dy
    integer(4), optional :: o_,m_,e_
    integer(4) :: o=1, m=3, e = 10                                              
    integer(4):: stat
    
    integer(4)::k,N
    real(rk):: xm,xp,ym,yp,bn(1:3,1:3)=0                 !bn: matrice dei coefficienti della comb lin
    if(present(o_)) o = o_                               !p: plus, m: minus
    if(present(m_)) m = m_
    if(present(e_)) e = e_

    bn(1,1)=1._rk
    bn(1,2)=4._rk/3._rk
    bn(2,2)=-1._rk/3._rk
    bn(1,3)=1.5_rk
    bn(2,3)=-0.6_rk
    bn(3,3)=0.1_rk
    
    N=size(x)-1
    if (Size(x)/=Size(y)) then                           !Controllo preliminare
       stat=1
       return
    end if
     if (l<0.or.l>N) then
       stat=2
       return
    end if
     if (o<1.or.o>2) then
       stat=3
       return
    end if
     if (m<1.or.m>3) then
       stat=4
       return
    end if

    dy=0                                                 !Algoritmo di derivazione
    Do k=1,m
       if(l<k) then                                      !Estrapola punti a sinistra
          xm=x(0) + (l-k)*(x(1)-x(0))
          call interp_aitken(x(0:Min(e-1,N)),y(0:Min(e-1,N)),xm,ym)
          xp=x(l+k)
          yp=y(l+k)
       else if (l>N-k)then                               !Estrapola punti a destra
          xm=x(l-k)
          ym=y(l-k)
          xp=x(N)+ (l-(N-k))*(x(N)-x(N-1))
          call interp_aitken(x(0:Max(0,N+1-e)),y(0:Max(0,N+1-e)),xp,yp)
       else
          xm=x(l-k)
          xp=x(l+k)
          ym=y(l-k)
          yp=y(l+k)
       end if
       
       Select case(o)                                    !Formula a seconda dell'ordine 'o' della derivata
       Case(1)
          dy=dy+bn(k,m)*(yp-ym)/(xp-xm)
       Case(2)
          dy=dy+bn(k,m)*(yp-2._rk*y(l)+ym)/(0.25_rk*(xp-xm)**2)
       End select
    End do
    stat=0
  End Function differentiate
  
!******************************************************************

  Recursive Function diff_fun(fun,x0,Dfun,o_,h_ini_,eps_,max_steps_,steps_) Result(stat)
  !!Compute the o-th derivative Dfun of the R->R function fun at x0 in steps_ steps
  !!starting with a h_ini increment and asking for a relative precision eps
  !!without overcoming max_steps iterative steps.
  !!input: funzione analitica, punto in cui calcolare la derivata e puntatore risultato Dfun.
  !!optional: o=1,2 ordine della derivata, h_ini=1 passo iniziale, eps=1.E-6_rk precisione,
  !!max_steps=1000 e steps puntatore al numero di step impiegati per raggiungere eps.
  !stat=0 success
  !stat=1 max_steps overcome
  !stat=2 diff. order not implemented

    Implicit None
    Procedure(fun_R_R) :: fun
    Real(rk), Intent(in) :: x0
    Real(rk), optional :: h_ini_, eps_
    Integer(4), optional :: o_, max_steps_, steps_
    Real(rk) :: h_ini=0.5_rk, eps=1.E-6_rk
    Integer(4) :: o=1, max_steps=1000
    Real(rk), Intent(out) :: Dfun
    Integer(4) :: stat

    Integer(4), Save :: steps=1
    Real(rk) :: Dyh,Dyh2
    Real(rk), Save :: Dyh_old
    
    if(present(o_)) o = o_                               
    if(present(max_steps_)) max_steps = max_steps_
    if(present(eps_)) eps = eps_
    if(present(h_ini_)) h_ini = h_ini_
    
    If (steps==1) Then                                                  !Inizializzazione 
       Select Case(o)
       Case(1)
          Dyh=(fun(x0+h_ini)-fun(x0-h_ini))/(2._rk*h_ini)
       Case(2)
          Dyh=(fun(x0+h_ini)-2._rk*fun(x0)+fun(x0-h_ini))/(h_ini**2)
       Case Default
          stat=2
          steps_=0
          steps=1
          Return
       End Select
    Else 
       Dyh=Dyh_old                                                      !Slittamento del valore di Dyh2 precedente
    End If
    
    Select Case(o)
    Case(1)
       Dyh2=(fun(x0+h_ini/2)-fun(x0-h_ini/2))/h_ini
    Case(2)
       Dyh2=(fun(x0+h_ini/2)-2._rk*fun(x0)+fun(x0-h_ini/2))/((h_ini**2)/4._rk)
    End Select
    Dyh_old=Dyh2                                                        !Mette nella variabile salvata il risultato

    If (Abs(Dyh-Dyh2)<=eps*Abs(Dyh2)) Then                              !Controllo per l'interruzione del programma
       Dfun=4._rk/3._rk*Dyh2-Dyh/3._rk
       stat=0
       If(present(steps_)) steps_=steps
       steps=1
       Return
    End If
    
    steps=steps+1
    If (steps>max_steps) Then
       stat=1
       If(present(steps_)) steps_=max_steps
       steps=1
       Return
    End If
    
    stat=diff_fun(fun,x0,Dfun,o,h_ini/2,eps,max_steps,steps_)           !Chiamata ricorsiva con la distanza dimezzata 
    
  End Function diff_fun

!******************************************************************
!-----------------------------------------------------------------
!******************************************************************

  Function integrate(x,y,integr,method_) result(stat)
  !!Compute the integral of the unknown function y=f(x) given 
  !!y(i)=f(x(i)) with i = 0...N by means of trapezoid or simpson method.
  !!if method is not present simpson is used by default.
  !!if N is odd (total number of points (N+1) is even) trapezoid is used.
  !!input: vettori punti e puntatore a var risultato.
  !!optional: method può essere simpson o trapezi.
  !!note: fornire un numero di punti dispari per simpson. 
  !!l'errore dei trapezi è O(h^2), mentre per simpson O(h^4).
  !stat = 0 success
  !stat = 1 Size(x)/=Size(y)
  !stat = 2 x points not equidistant
  !stat = 3 requested simpson for even total number of points, changed to trapezoid
  !stat = 4 method not implemented

    Implicit None
    Real(rk), Intent(In) :: x(0:), y(0:)
    Character(7), optional :: method_
    Real(rk) :: integr
    Integer(4):: stat

    Integer(4) ::N,k
    Character(7) :: method
    Real(rk) :: w(0:Size(x)-1), Delta
    delta=x(1)-x(0)
    N=Size(x) - 1

    If (Size(x)/=Size(y)) Then
       stat=1
       Return
    End if
    
    If (Any((x(N:1)-x(N-1:0))/=delta)) Then
       stat=2
       Return
    End if
    
    If (present(method_)) Then
       method=method_
    Else
       If (mod(N,2)==0) Then
          method='simpson'
       Else
          method='trapezi'
       End If
    End if
    
    If (mod(N,2)==1.And.method=='simpson') Then
       method='trapezi'
       stat=3
    End if

    Select Case(trim(method))                                          !Coefficienti del metodo d'integrazione
    Case('simpson')
       Forall(k=0:N,mod(k,2)==0) w(k)=2._rk/3._rk
       Forall(k=0:N,mod(k,2)==1) w(k)=4._rk/3._rk
       w(0)=1._rk/3._rk
       w(N)=1._rk/3._rk
    Case('trapezi')
       w=1._rk
       w(0)=0.5_rk
       w(N)=w(0)                                           
    Case default
       stat = 4
       return
    End Select
    
    integr = delta*dot_product(w,y)
    If (stat/=3) stat=0

  End function integrate

!**********************************************************************

   Recursive Function int_fun(fun,a,b,fa,fb,Ifun,eps_,min_steps_,max_steps_,steps_,fcc) Result(stat)
   !!Compute the integral Ifun of the R->R function fun between a and b in steps_ steps in minimum 
   !!min_steps and asking for a relative precision eps without overcoming max_steps iterative steps.
   !!fa=fun(a) and fb=fun(b) are given to avoid recomputing and fcc = fun(c) is optional and is given 
   !!if simpson method should be used instead of trapezoid method.
   !!input: funzione analitica, estremi intervallo, funzione ivi valutata e puntatore risultato.
   !!optional: eps=1.E-6_rk precisione, min_steps=10, max_steps=1000, steps puntatore al numero di step 
   !!impiegati per raggiungere eps e fcc = fun(c).
   !stat=0 success
   !stat=1 max_steps overcome
    
    Implicit None
    Procedure(fun_R_R) :: fun
    Real(rk), Intent(In) :: a,b,fa,fb
    Real(rk), Intent(In), optional :: fcc,eps_
    Integer(4), optional :: min_steps_,max_steps_, steps_
    Real(rk), Intent(Out) :: Ifun
    Integer(4) :: stat

    Integer(4), Save :: steps=1
    Real(rk), Save :: a0,b0
    Integer(4) :: min_steps=10, max_steps=10000
    Real(rk) :: d,h,c,e,fc,fd,fe,S_I,S_II,Ifun_I,Ifun_II,alpha,eps=1.E-6_rk
    
    If(present(eps_)) eps = eps_
    If(present(min_steps_)) min_steps = min_steps_
    If(present(max_steps_)) max_steps = max_steps_

    If (steps==1) Then                                    !Inizializzazione con assegnazione degli estremi 
       a0=a
       b0=b
    End if
    
    h=b-a                                                 !Salviamo la distanza
    c=0.5_rk*(a+b)                                        !Il punto medio
    
    If (.not.present(fcc)) Then
       S_I=0.5_rk*h*(fa+fb)
       fc=fun(c)
       S_II=0.5_rk*S_I+0.5_rk*h*fc
       alpha=3._rk
    Else
       fc=fcc
       S_I=h*(fa+4._rk*fc+fb)/6._rk
       d=0.5_rk*(a+c)
       fd=fun(d)
       e=0.5_rk*(c+b)
       fe=fun(e)
       S_II=h*(fa+4._rk*fd+fc)/12._rk+h*(fc+4._rk*fe+fb)/12._rk
       alpha=15._rk
    Endif
    
    If ((Abs(S_II-S_I)<=alpha*eps*Abs(S_II)).And.(steps>=min_steps)) Then
       Ifun=S_II                                          !Se c'è la precisione richiesta salva e stat=0
       stat=0
    Else
       steps=steps+2
       If (steps>max_steps) Then                          !Controllo sul numero di step
          Ifun=0._rk
          stat=1
       Else                                               !chiamata iterativa ai due intervalli dx e sx
          If (.not.present(fcc)) Then
             stat=Min(1,int_fun(fun,a,c,fa,fc,Ifun_I,eps,min_steps,max_steps,steps_)&
                  +int_fun(fun,c,b,fc,fb,Ifun_II,eps,min_steps,max_steps,steps_))
          Else
             stat=Min(1,int_fun(fun,a,c,fa,fc,Ifun_I,eps,min_steps,max_steps,steps_,fd)&
                  +int_fun(fun,c,b,fc,fb,Ifun_II,eps,min_steps,max_steps,steps_,fe))
          Endif
          If (stat==0) Ifun=Ifun_I+Ifun_II
       End If
    Endif    
    
    If (a==a0.And.b==b0) Then                             !Il calcolo non è neanche iniziato 
       If(present(steps_)) steps_=steps
       steps=1
    Endif
    
  End Function int_fun
  
!!=================================================================
!!=================================================================
!!==       Roots and ordinary differential equations           ====
!!=================================================================
!!=================================================================

  Function braket(fun,x,dx,a,b) Result (stat)
  !!Compute the minimal interval [a,b] symmetric around the point x
  !!that contains a zero of the function fun without exceeding the
  !!original interval [a,b] and using dx has incremental step.
  !stat = 0 success
  !stat = 1 no zero between a and b

    Implicit None
    Procedure(fun_R_R) :: fun
    Real(rk), Intent(In) :: x,dx
    Real(rk), Intent(InOut) :: a,b
    Integer(4) :: stat

    Real(rk) :: a0,b0

    a0=x-dx
    b0=x+dx
    Do While (a0>=a.And.b0<=b)
       If (fun(a0)*fun(b0)<0) Then
          a=a0
          b=b0
          stat=0
          Return
       Endif
       a0=a0-dx
       b0=b0+dx
    Enddo
    stat=1
  End Function braket

!**********************************************************************

  Function bisection(fun,x,root,dx_,a_,b_,eps_fun_,eps_,max_steps_,steps_) Result (stat)
  !!Compute the zero of the function fun closer to x by means of bisection in steps_ steps without 
  !!exceeding max_steps. The zero is returned when fun(zero)<=eps_fun or |a_n-b_n|<=eps for n-step.
  !!input: puntatore alla funzione, x vicino al presunto zero, root puntatore al risultato.
  !!optional: eps_fun=1.E-6_rk, eps=1.E-6_rk, max_steps=100, steps puntatore al numero di step,  
  !!dx=0.001 incremento e massimo intervallo per la ricerca [a=-10,b=10].
  !stat = 0 success
  !stat = 1 no zero between a and b
  !stat = 2 max_steps exceeded

    Implicit None
    Procedure(fun_R_R) :: fun
    Real(rk), Intent(In) :: x
    Real(rk), Intent(Out) :: root
    Integer(4), optional :: max_steps_,steps_
    Real(rk), Intent(In), Optional :: dx_,a_,b_,eps_fun_,eps_
    Integer(4) :: stat

    Real(rk) :: a,b,aa,bb,c,fa,fc,dx=0.001_rk,eps_fun=1.E-6_rk,eps=1.E-6_rk
    Integer(4) :: braket_stat,steps,max_steps=100
    
    If(present(eps_)) eps = eps_
    If(present(eps_fun_)) eps_fun = eps_fun_
    If(present(dx_)) dx = dx_
    If(present(max_steps_)) max_steps = max_steps_

    steps=1
    If (present(a_).And.present(b_)) Then
       a=a_
       b=b_
    Else
       a=x-10._rk
       b=x+10._rk
    Endif
    
    braket_stat=braket(fun,x,dx,a,b)
    
    If (braket_stat/=0) Then
       stat=1
       If(present(steps_)) steps_=steps
       return
    Endif
    
    aa=a
    bb=b
    fa=fun(aa)
    Do While (steps<=max_steps)
       c=(aa+bb)/2._rk
       fc=fun(c)
       If ((Abs(fc)<=eps_fun).Or.(Abs(aa-bb)<=eps)) Then
          root=c
          stat=0
          If(present(steps_)) steps_=steps
          return
       Endif
       If (fa*fc<0) Then
          bb=c
       Else
          aa=c
          fa=fc
       Endif
       steps=steps+1
    Enddo
    
    root=c
    If(present(steps_)) steps_=steps
    stat=2
    
  End Function bisection

!**********************************************************************

  Function secant(fun,x,zero,dx_,eps_fun_,eps_,max_steps_,dxx_,a_,b_,steps_bis_,steps_) Result (stat)
  !!Compute the zero of the function fun closer to x by means of secant in steps_ steps without 
  !!exceeding max_steps. The zero is returned when fun(zero)<=eps_fun or |a_n-b_n|<=eps for n-step.
  !!input: puntatore alla funzione, x vicino al presunto zero, zero puntatore al risultato.
  !!optional: eps_fun=1.E-6_rk, eps=1.E-6_rk, max_steps=100, steps puntatore al numero di step,  
  !!dx=0.001 incremento del secondo punto.
  !!note: one can start with maximum steps_bis of bisection in the interval [a,b] with pass dxx.
  !stat = 0 success
  !stat = 1 no zero between a and b
  !stat = 2 max_steps exceeded

    Implicit None
    Procedure(fun_R_R) :: fun
    Real(rk), Intent(In) :: x
    Real(rk), Intent(Out) :: zero
    Integer(4), optional :: max_steps_,steps_,steps_bis_
    Real(rk), Intent(In), optional :: dx_,dxx_,a_,b_,eps_fun_,eps_
    Integer(4) :: stat

    Real(rk) :: x_,x__,x___,fx_,dfx__,fx__,fx___,zero_bis,dx=0.001,eps_fun=1.E-6_rk,eps=1.E-6_rk
    Integer(4) :: bisection_stat,steps,max_steps=100
    
    If(present(eps_)) eps = eps_
    If(present(eps_fun_)) eps_fun = eps_fun_
    If(present(dx_)) dx = dx_
    If(present(max_steps_)) max_steps = max_steps_

    steps=1
    If(present(steps_bis_)) Then
       If (present(dxx_).and.present(a_).and.present(b_)) Then
          bisection_stat=bisection(fun,x,zero_bis,dxx_,a_,b_,eps_fun,eps,steps_bis_,steps)
       Else
          bisection_stat=bisection(fun,x,zero_bis,eps_fun_=eps_fun,eps_=eps,max_steps_=steps_bis_,steps_=steps)
       End if
       
       If(bisection_stat==0) Then
          If(present(steps_)) steps_=0
          zero=zero_bis
          stat=0
          return
       Else If(bisection_stat==1) Then
          If(present(steps_)) steps_=0
          zero=zero_bis
          stat=1
          return
       End if
       x_=zero_bis
    Else
       x_=x
    Endif
    
    x__=x_+dx
    fx_=fun(x_)
    fx__=fun(x__)
    Do While (steps<=max_steps)
       dfx__=(fx__-fx_)/(x__-x_)
       x___=x__-fx__/dfx__
       fx___=fun(x___)
       If ((Abs(fx___)<=eps_fun).Or.(Abs(x___-x__)<=eps)) Then
          zero=x___
          stat=0
          If(present(steps_)) steps_=steps
          Return
       Endif
       x_=x__
       x__=x___
       fx_=fx__
       fx__=fx___
       steps=steps+1
    Enddo
    
    zero=x___
    If(present(steps_)) steps_=steps
    stat=2
    
  End Function secant

!**********************************************************************

  Function newton(fun,dfun,x,zero,eps_fun_,eps_,max_steps_,dxx_,a_,b_,steps_bis_,steps_) Result (stat)
  !!Compute the zero of the function fun closer to x by means of secant in steps_ steps without 
  !!exceeding max_steps. The zero is returned when fun(zero)<=eps_fun or |a_n-b_n|<=eps for n-step.
  !!dfun is the derivative function of fun.
  !!input: puntatore alla funzione e alla sua derivata, x vicino al presunto zero, zero puntatore al risultato.
  !!optional: eps_fun=1.E-6_rk, eps=1.E-6_rk, max_steps=100, steps puntatore al numero di step,  
  !!note: one can start with maximum steps_bis of bisection in the interval [a,b] with pass dxx.
  !stat = 0 success
  !stat = 1 no zero between a and b
  !stat = 2 max_steps exceeded

    Implicit None
    Procedure(fun_R_R) :: fun,dfun
    Real(rk), Intent(In) :: x
    Real(rk), Intent(Out) :: zero
    Integer(4), optional :: max_steps_,steps_,steps_bis_
    Real(rk), Intent(In), optional :: dxx_,a_,b_,eps_fun_,eps_
    Integer(4) :: stat

    Real(rk) :: x_,x__,fx_,dfx_,zero_bis,dx=0.001,eps_fun=1.E-6_rk,eps=1.E-6_rk
    Integer(4) :: bisection_stat,steps,max_steps=100
    
    If(present(eps_)) eps = eps_
    If(present(eps_fun_)) eps_fun = eps_fun_
    If(present(max_steps_)) max_steps = max_steps_

    steps=1
    If(present(steps_bis_)) Then
       If (present(dxx_).and.present(a_).and.present(b_)) Then
          bisection_stat=bisection(fun,x,zero_bis,dxx_,a_,b_,eps_fun,eps,steps_bis_,steps)
       Else
          bisection_stat=bisection(fun,x,zero_bis,eps_fun_=eps_fun,eps_=eps,max_steps_=steps_bis_,steps_=steps)
       End if
       
       If(bisection_stat==0) Then
          If(present(steps_)) steps_=0
          zero=zero_bis
          stat=0
          return
       Else If(bisection_stat==1) Then
          If(present(steps_)) steps_=0
          zero=zero_bis
          stat=1
          return
       End if
       x_=zero_bis
    Else
       x_=x
    Endif
    
    fx_=fun(x_)
    dfx_=dfun(x_)
    Do While (steps<=max_steps)
       x__=x_-fx_/dfx_
       fx_=fun(x__)
       If ((Abs(fx_)<=eps_fun).Or.(Abs(x__-x_)<=eps)) Then
          zero=x__
          stat=0
          If(present(steps_)) steps_=steps
          return
       Endif
       x_=x__
       dfx_=dfun(x_)
       steps=steps+1
    Enddo
    
    zero=x__
    If(present(steps_)) steps_=steps
    stat=2

  End Function newton
  
!******************************************************************
!-----------------------------------------------------------------
!******************************************************************
  
  Function Euler(Data) Result (stat) !Predictor
  !!Compute data%y at point data%x0+data%dx starting from data%y0
  !!at point data%x0 given the current data%fun by Euler (Taylor-O(tau^2)) method
  !stat = 0 success
  !stat = 1 data%y0 and data%y are not congruent

    Implicit None
    Class(ode_data) :: Data
    Integer(4) :: stat

    If (Size(Data%y0)/=Size(Data%y)) Then
       stat=1
       return
    End If
    
    Data%y=Data%y0+Data%dx*Data%fun(Data%y0,Data%x0)
    stat=0
    
  End Function euler

!**********************************************************************

  Function Taylor3(Data) Result (stat) !Predictor
  !!Compute data%y at point data%x0+data%dx starting from data%y0
  !!at point data%x0 given the current data%fun by Taylor-O(tau^3) method
  !stat = 0 success
  !stat = 1 data%y0 and data%y are not congruent

    Implicit None
    Class(ode_data) :: Data
    Integer(4) :: stat

    If (Size(Data%y0)/=Size(Data%y)) Then
       stat=1
       return
    End If
    
    Data%y=Data%y0+Data%dx*Data%fun(Data%y0,Data%x0)+0.5_rk*((Data%dx)**2)*(Data%dfun_t(Data%y0,Data%x0)&
         +Matmul(Data%dfun_y(Data%y0,Data%x0),Data%fun(Data%y0,Data%x0)))
    stat=0
    
  End Function taylor3

!**********************************************************************

  Function I2(Data) Result (stat) !Corrector
  !!Compute data%y at point data%x0+data%dx starting from data%y0
  !!at point data%x0 given the current data%fun by I2 (non-iterative Picard) method
  !!To be used as corrector as it uses data%y too
  !stat = 0 success
  !stat = 1 data%y0 and data%y are not congruent

    Implicit None
    Class(ode_data) :: Data
    Integer(4) :: stat

    If (Size(Data%y0)/=Size(Data%y)) Then
       stat=1
       Return
    End If
    
    Data%y=Data%y0+0.5_rk*Data%dx*(Data%fun(Data%y0,Data%x0)+Data%fun(Data%y,Data%x0+Data%dx))
    stat=0
    
  End Function I2

!**********************************************************************

  Function D3(Data) Result (stat)
  !!Compute data%y at point data%x0+2*data%dx starting from data%y0
  !!at point data%x0 given the current data%fun by D3 method
  !!It needs data%y1 at point data%x0+data%dx too
  !stat = 0 success
  !stat = 1 data%y0 and data%y are not congruent
  !stat = 2 data%y1 and data%y are not congruent

    Implicit None
    Class(ode_data) :: Data
    Integer(4) :: stat

    If (Size(Data%y0)/=Size(Data%y)) Then
       stat=1
       Return
    End If
    
    If (Size(Data%y1)/=Size(Data%y)) Then
       stat=2
       Return
    End If
    
    Data%y=Data%y0+2._rk*Data%dx*Data%fun(Data%y1,Data%x0+Data%dx)
    stat=0
    
  End Function D3
  
!**********************************************************************


  Function I3(Data) Result (stat) !Corrector
  !!Compute data%y at point data%x0+2*data%dx starting from data%y0
  !!at point data%x0 given the current data%fun by I3 (Simpson) method
  !!It needs data%y1 at point data%x0+data%dx too
  !!To be used as corrector as it uses data%y too
  !stat = 0 success
  !stat = 1 data%y0 and data%y are not congruent
  !stat = 2 data%y1 and data%y are not congruent

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
  
!**********************************************************************
  
  Function IterPicard (data,para) result (stat) !Corrector
  !!Compute data%y by iterative Picard method
  !!To be used as corrector as it uses data%y too
  !stat = 0 success
  !stat = 1 max_steps exceeded
  
    Implicit None
    Type(ode_data):: data
    Type(ode_para):: para
    Integer(4):: stat

    real(rk)::y_n(1:size(data%y)),y_o(1:size(data%y))

    y_n=data%y
    para%steps=0
    
    do while (any(abs(y_n-y_o)>para%eps).or.(para%steps==0))
       y_o=y_n
       para%steps=para%steps+1
       
       if (para%steps>para%max_steps) then
         stat=1
         return
       end if
       
       y_n=data%y0+0.5_rk*data%dx*(data%fun(data%y0,data%x0)+data%fun(y_o,data%x0+data%dx))
    end do
    
    data%y=y_n
    stat=0
    
  End Function IterPicard

!**********************************************************************

  Function Runge_Kutta(Data) Result (stat)
  !!Compute data%y at point data%x0+data%dx starting from data%y0
  !!at point data%x0 given the current data%fun by Runge_Kutta method
  !stat = 0 success
  !stat = 1 data%y0 and data%y are not congruent


    Implicit None
    Class(ode_data) :: Data
    Integer(4) :: i, stat

    Real(rk) :: c(0:4,1:size(data%y)), a(1:4)=[1,2,2,1]

    If (size(Data%y0)/=size(Data%y)) then
       stat=1
       return
    End If

    c(0,:)=0
    do i=1,4
       c(i,:)=data%dx*data%fun(data%y0 + c(i-1,:)/a(i),data%x0 + data%dx/a(i))
    end do
    
    data%y=data%y0+matmul(a,c(1:,:))/6._rk
    stat=0
    
  End Function Runge_Kutta
  
!**********************************************************************
  
   Function shooting(Data,para) Result (stat)
   !!Compute u(x0->x1)[=data%u(1,x0->x1)] and du(x0->x1)[=data%u(2,x0->x1)]
   !!with the shooting method (with secant and Runge-Kutta)
   !!for the three cases below where (x1>x0 and dx>0) or  (x1<x0 and dx<0)
   !!para%case = 'A' known  u(x0)[=y0(1)] e  u(x1)[=y1(1)] - 
   !!note posizione iniziale e finale (determina velocità iniziale)
   !!para%case = 'B' known  u(x0)[=y0(1)] e du(x1)[=y1(2)] - 
   !!note posizione iniziale e velocità finale (determina velocità iniziale)
   !!para%case = 'C' known du(x0)[=y0(2)] e du(x1)[=y1(2)] - 
   !!note velocità iniziale e velocità finale (determina posizione iniziale)
   !stat = 0 success
   !stat = 1 no zero between ini-10 and ini+10
   !stat = 2 max_steps exceeded in secant

    Implicit None
    Type(ode_data), Target :: Data
    Type(ode_para) :: para
    Integer(4) :: stat

    Real(rk) :: y0_(1:Size(Data%y0)),x0_
    Real(rk), Pointer :: var,com
    Real(rk) :: ini,exa,zero
    Integer(4) :: i,max_i

    max_i=Ceiling(Abs((Data%x1-Data%x0)/Data%dx)) + 1     !Ceiling prende l'intero più piccolo
    Allocate(Data%u(1:Size(Data%y0),1:max_i))

    y0_=Data%y0
    x0_=Data%x0

    Select Case(para%Case)
    Case('A')
       ini = Data%y0(2)
       var => Data%y0(2)
       exa = Data%y1(1)
       com => Data%y(1)
    Case('B')
       ini = Data%y0(2)
       var => Data%y0(2)
       exa = Data%y1(2)
       com => Data%y(2)
    Case('C')
       ini = Data%y0(1)
       var => Data%y0(1)
       exa = Data%y1(2)
       com => Data%y(2)
    End Select

    stat=secant(aux,ini,zero,dx_=Data%dx,eps_fun_=para%eps,eps_=para%eps,max_steps_=para%max_steps,steps_bis_=0,steps_=para%steps)
    
    If (stat==1) Then
       Print *,'No zero between:',ini-10,' and ',ini+10
       Return
    End If
    If (stat==2) Then
       Print *,'Exceeded max_steps in Secant'
       Return
    End If

    Data%y0=y0_
    Data%x0=x0_
    var=zero
    Data%u(:,1)=Data%y0
    Do i=2,max_i
       stat=Runge_Kutta(Data)
       Data%u(:,i)=Data%y
       Data%y0=Data%y
       Data%x0=Data%x0+Data%dx
    Enddo
    Data%y0=y0_
    Data%x0=x0_
    var=zero

  Contains

    Function aux(x)
      Implicit None
      Real(rk), Intent(in) :: x
      Real(rk) :: aux

      Data%y0=y0_
      Data%x0=x0_
      var=x       
      Do i=2,max_i
         stat=Runge_Kutta(Data)
         Data%y0=Data%y
         Data%x0=Data%x0+Data%dx
      Enddo
      aux=com-exa
    End Function aux

  End Function shooting
  
!!=================================================================
!!=================================================================
!!==         Linear algebra and eigenvalues problem            ====
!!=================================================================
!!=================================================================

   Function tridiagsup(Data) Result (stat)
   !!Expresse in tridiagonal form the matrix A taking care of moving 
   !!the coefficinet vector b and matrix Bm too (if allocated).
   !stat = 0 success
   !stat = 1 matrix A is not square
   !stat = 2 proportional rows
   !stat = 3 more than one diagonal element is equal to zero
   !stat = 4 b is incompatible with A
   !stat = 5 Bm is incompatible with A

    Implicit None
    Type(matrix) :: Data
    Integer(4) :: stat

    Integer(4) :: n,i,j,k,q,r,bp,bm
    Real(rk) :: pivot

    n=Size(Data%A,1)

    If (n/=Size(Data%A,2)) Then                            !Controllo sul fatto che la matrice sia quadrata
       stat=1
       Return
    End If

    If (Allocated(Data%At)) Deallocate(Data%At)
    Allocate(Data%At(1:n,1:n))
    Data%At=Data%A

    If (Allocated(Data%v)) Deallocate(Data%v)
    Allocate(Data%v(1:n))
    Forall(i=1:n) Data%v(i)=i

    bp=0                                                   !bp fa da switch per l'uso del vettore termini noti
    If (Allocated(Data%b)) Then                            
       If (Size(Data%b)/=n) Then                           !Controllo dimensione di b
          stat=4
          Return
       End If
       bp=1
       If (Allocated(Data%bt)) Deallocate(Data%bt)
       Allocate(Data%bt(1:n))
       Data%bt=Data%b       
    Endif

    bm=0                                                   !bm fa da switch per l'uso della matrice Bm
    If (Allocated(Data%Bm)) Then
       If (Size(Data%Bm,1)/=n.Or.Size(Data%Bm,2)/=n) Then  !Controllo nel primo e nel secondo indice di Bm Size(cosa, indice)
          stat=5                                           !per verificare se è quadrata
          Return
       End If
       bm=1
       If (Allocated(Data%Bmt)) Deallocate(Data%Bmt)
       Allocate(Data%Bmt(1:n,1:n))
       Data%Bmt=Data%Bm
    Endif

    Data%p=1                                               !p tiene conto del segno delle permutazioni
    Do k=1,n-1
       If (Any(Maxval(Abs(Data%At(Data%v(k:),:)),2)<1.D-12)) Then
          stat=2                                           !Controlla se il maxval dalla k riga in poi al variare del secondo indice
          Return                                           !(scorrendo le colonne della riga) è zero, e quindi se ci sono righe proporz.
       Endif
       If (Maxval(Abs(Data%At(Data%v(k:),k)))<1.D-12) Then
          stat=3                                           !Controlla se un elem in diagonale è nullo
          Return
       Endif
       q=Maxloc(Abs(Data%At(Data%v(k:),k))/Maxval(Abs(Data%At(Data%v(k:),:)),2),1)+(k-1) !+(k-1) perché maxloc è applicato a vettore ridotto
       If (q/=k) Then                                      !q indica la riga della colonna col massimo valore ????????
          r=Data%v(q)
          Data%v(q)=Data%v(k)                              !scambio di righe tramite v
          Data%v(k)=r
          Data%p=-Data%p                                   !cambio segno
       Endif

       
       Do i=k+1,n
          pivot=Data%At(Data%v(i),k)                       
          Do j=1,n
             If (bm==1) Data%Bmt(Data%v(i),j)=Data%Bmt(Data%v(i),j)-pivot/Data%At(Data%v(k),k)*Data%Bmt(Data%v(k),j)
             !Se lo switch di Bm è attivo, allora triangola anche Bm
             If (j<k) Cycle
             Data%At(Data%v(i),j)=Data%At(Data%v(i),j)-pivot/Data%At(Data%v(k),k)*Data%At(Data%v(k),j)
          Enddo
          If (bp==1) Data%bt(Data%v(i))=Data%bt(Data%v(i))-pivot/Data%At(Data%v(k),k)*Data%bt(Data%v(k))
          !Se lo switch di b è attivo, allora triangola anche b
       Enddo
    Enddo
    stat=0

  End Function tridiagsup

!**********************************************************************

  Function linsys(Data) Result (stat)
  !!Compute the solution x of linear system equations given by A*x=b
  !!in which A's rows are linearly independents.
  !stat = 0 success
  !stat = 1 matrix A is not square
  !stat = 2 proportional rows
  !stat = 3 more than one diagonal element is equal to zero
  !stat = 4 b is incompatible with A
  !stat = 5 b is not allocated or incompatible with A
  !stat = 6 data%At(data%v(n),n) = 0

    Implicit None
    Type(matrix) :: Data
    Integer(4) :: stat

    Integer(4) :: n,i,j
    Real(rk) :: sum

    n=Size(Data%A,1)

    If (n/=Size(Data%A,2)) Then                             !Controlla se la matrice è quadrata
       stat=1
       Return
    End If

    If (.Not.Allocated(Data%b).Or.Size(Data%b)/=n) Then     !Verifica presenza dei termini noti
       stat=5
       Return
    End If

    stat=tridiagsup(Data)                                   !Effettua la riduzione a scalini
    If (stat/=0) Return    

    If (Data%At(Data%v(n),n)==0) Then                       !Se l'ultimo termine è nullo teniamo una serie di problemi
       stat=6                                               !O sistema incompatibile o infinite soluzioni
       Return
    Endif

    If (Allocated(Data%x)) Deallocate(Data%x)               !Vettore delle incognite
    Allocate(Data%x(1:n))

    Data%x(n)=Data%bt(Data%v(n))/Data%At(Data%v(n),n)       !Risolve l'ultima riga del sistema
    Do i=n-1,1,-1                                           !Itera dall'ultima alla prima riga
       sum=0
       Do j=i+1,n
          sum=sum+Data%At(Data%v(i),j)*Data%x(j)            !Calcola il valore della riga a meno dell'incognita
       Enddo
       Data%x(i)=(Data%bt(Data%v(i))-sum)/Data%At(Data%v(i),i) !Trova l'incognita
    Enddo
    stat=0

  End Function linsys

!**********************************************************************

  Function det(Data) result(stat)
  !!Compute the determinant of the matrix A.
  !stat = 0 success
  !stat = 1 matrix A is not square
  !stat = 2 proportional rows
  !stat = 3 more than one diagonal element is equal to zero

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
    If (stat/=0) return    

    Data%detA=Data%p
    Do i=1,n
       Data%detA=Data%detA*Data%At(Data%v(i),i)
    Enddo
    stat=0

  End Function det
  
!**********************************************************************

  Function trace(Data) result(stat)
  !!Compute the determinant of the matrix A.
  !stat = 0 success
  !stat = 1 matrix A is not square
  
  Implicit None
    Type(matrix) :: Data
    Integer(4) :: stat

    Integer(4) :: i,n

    n=Size(Data%A,1)

    If (n/=Size(Data%A,2)) Then
       stat=1
       return
    End If
    
    Data%TrA=0._rk
    Do i=1,n
       Data%TrA=Data%TrA + Data%A(i,i)  
    End do
    
    stat = 0
  
  End Function trace

!**********************************************************************

  Function inv(Data) Result (stat)
  !!Compute the inverse of the matrix A.
  !stat = 0 success
  !stat = 1 matrix A is not square
  !stat = 2 proportional rows
  !stat = 3 more than one diagonal element is equal to zero

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

    stat=det(data)
    if (abs(Data%detA)<1.E-12_rk) stop 'Determinante nullo' 

    If (Allocated(Data%Bm)) Deallocate(Data%Bm)                  !Genera la matrice identità
    Allocate(Data%Bm(1:n,1:n))
    Data%Bm=0._rk
    Do i=1,n
       Data%Bm(i,i)=1._rk
    Enddo

    stat=tridiagsup(Data)                                        !Riduce a scalini considerando anche Bm
    If (stat/=0) Return    

    If (Allocated(Data%iA)) Deallocate(Data%iA)                  !Alloca la matrice inversa
    Allocate(Data%iA(1:n,1:n))

    Do j=1,n
       Data%iA(n,j)=Data%Bmt(Data%v(n),j)/Data%At(Data%v(n),n)   !Risolve l'ultima riga del sistema per ogni colonna di Bmt
    Enddo
    Do i=n-1,1,-1                                                !Itera al contrario come linsys
       Do j=1,n                                                  !Varia la colonna che stiamo risolvendo
          sum=0
          Do k=i+1,n
             sum=sum+Data%At(Data%v(i),k)*Data%iA(k,j)           !Come linsys 
          Enddo
          Data%iA(i,j)=(Data%Bmt(Data%v(i),j)-sum)/Data%At(Data%v(i),i)
       Enddo
    Enddo
    stat=0

  End Function inv

!**********************************************************************

  Function eigen(data, para) result (stat)
  !!Compute the highest eigenvalue (i.e. eigenvector)
  !stat = 0 success
  !stat = 1 matrix A is not square
  !stat = 2 proportional rows
  !stat = 3 more than one diagonal element is equal to zero
  !stat = 4 superato max steps
    
    Implicit none
    Type(matrix):: data
    Type(matrix_para) ::para

    Real (rk)::lambda,c(1:size(data%A,1)),v(1:size(data%A,1))
    Integer(4):: stat

    if (allocated(data%x)) deallocate(data%x)
    allocate(data%x(1:size(data%A,1)))
    call random_seed()                                           !Generazione casuale di num tra 0 e 1
    call random_number(data%x)
    data%x=data%x/sqrt(dot_product(data%x,data%x))
    data%lambda=0._rk
    para%steps=0

    Do while (para%steps==0.or.(Any(abs(data%x-v)>para%eps_x).and.Any(abs(data%x+v)>para%eps_x))&
         .or.(data%lambda-lambda)>para%eps_lambda)               !Controllo sulla precisione di auto-vett/val
      v=data%x
      lambda=data%lambda 
      c=matmul(data%A,v)
      data%x=c/sqrt(dot_product(c,c))
      data%lambda=dot_product(data%x,matmul(data%A,data%x))/&
           sqrt(dot_product(data%x,data%x))
      
      para%steps=para%steps+1
      if (para%steps>para%max_steps)then
         stat=4
         return
      end if
    End do
    stat=0
    
  End Function eigen
  
!**********************************************************************

   Function Gershgorin(data,u) result (stat)
   !!Compute the center and the radius of the Gershgorin circles, verifies if there's
   !!an intersection and then it returns the vector of the n-2 intermediate circle's centers.
   !!It works only for more than 3 circles. 
   !stat = 0 success
   !stat = 1 two circles intersect
   !stat = 4 size < 4
    
   Implicit none
   Type(matrix):: data
   Type(matrix_para) ::para
   Real(rk):: u(3:size(data%A,1))

   Real(rk):: v(1:size(data%A,1),1:2),a,b,c(1:2)
   Integer(4):: stat,i,j,n

   stat = 0
   n = size(data%A,1)
   If (n < 4) then
      stat = 4
      return
    end if

   Do i=1,n                       !Elementi sulla diagonale(centro) e distanza(raggio) dell'autovalore associato
       v(i,1)=data%A(i,i)
       a=0._rk
       b=a
       do j=1,n
          a=a+abs(data%A(i,j))
          b=b+abs(data%A(j,i))
       end do
       v(i,2)=Minval([a,b])-abs(v(i,1))
   End do

   Do j=n-1,1,-1                   !Dispone i centri e i raggi in ordine decrescente (dei centri)
     Do i=1,j
      if (v(i,1)<v(i+1,1)) then
         c=v(i,:)
         v(i,:)=v(i+1,:)
         v(i+1,:)=c
      end if
     End do
   End do

   Do i=1,n-1                      !Verifica se i cerchi si intersecano
       if ((v(i,1)-v(i,2))-(v(i+1,1)+v(i+1,2))<0._rk) then
         Print '(x,A,x,(x,E13.6),x,A,x,(x,E13.6),x,A)', 'ATTENTION: The circles in ',v(i,1),' and ',v(i+1,1), ' intersect.'
         stat = 1
       End if
   End do
   
   Print *,''
   Print *, 'Gershgorin'
   Do i=1,n
       Print '(2(x,E13.6))', v(i,:)
   End do
   Print *,''

   Do j=n-1,1,-1                   !Dispone i centri e i raggi in ordine decrescente in modulo(dei centri)
     Do i=1,j
      if (abs(v(i,1)) < abs(v(i+1,1))) then
         c=v(i,:)
         v(i,:)=v(i+1,:)
         v(i+1,:)=c
      end if
     End do
   End do
   
   u=v(2:n-1,1)                    !Consideriamo i coefficienti di quelli compresi tra il min e il max autoval in modulo
   
  End Function Gershgorin
  
!**********************************************************************
  
  Function diagonalizationR(data, para, lambda, vlamb) result (stat)
  !!Compute eigenvectors and eigenvalues of a real matrix with non-intersecting Gershgorin circles.
  !!The allowed dimension is nxn with n>3.
  !stat = 0 success
  !stat = 1 determinant or inverse matrix are not defined
  !stat = 2 some circles intersect
  !stat = 3 eigen failed
  !stat = 4 size < 4
  
    Implicit none
    Type(matrix) :: data
    Type(matrix_para) :: para
    Real(rk) :: lambda(1:size(data%A,1)), vlamb(1:size(data%A,1),1:size(data%A,1))
    
    Real(rk) :: A(1:size(data%A,1),1:size(data%A,1)), Uni(1:size(data%A,1),1:size(data%A,1)), u(3:size(data%A,1)), z, y
    Integer(4):: stat, i, j, n
    
    n = size(data%A,1)
    If (n<4) then
      stat = 4
      return
    end if
    If(det(data)/=0) then
      stat = 1
      return
    end if
    If(inv(data)/=0) then
      stat = 1
      return
    end if
    If(Gershgorin(Data,u)/=0) then
      print *, 'The intersection may affect the research giving several times the same eigenvalue.'
      stat = 2
    end if
  
    A = data%A                                                   !Save the value of data%A
  
    If(eigen(Data,para)/=0) then                                 !Compute eigenvector and eigenvalue (highest)
      print *, 'Highest eigen failed!'
      stat = 3
      return
    end if
    lambda(1) = data%lambda
    vlamb(1,:) = data%x

    data%A=data%iA
    If(eigen(Data,para)/=0) then                                 !Compute eigenvector and eigenvalue (ground)
      print*, 'Lowest eigenvalue failed!'     
      stat = 3
      return
    end if
    lambda(2)=1._rk/data%lambda
    vlamb(2,:) = data%x
  
    Forall(i=1:n, j=1:n, i/=j) Uni(i,j)=0._rk
    Forall(i=1:n) Uni(i,i)=1._rk
    Do i=3,n-2                                                   !Compute n-4 eigenvectors and eigenvalues
       data%A=A-u(i)*Uni
       If(inv(Data)/=0) then
         print*, 'Error inverting the matrix associated to the eigenvalue!'     
         stat = 3
         return
       end if
       data%A=data%iA
       If(eigen(Data,para)/=0) then
         print*, 'Eigenvalue failed!'    
         stat = 3
         return
       end if
       lambda(i)=(1._rk/data%lambda)+u(i)
       vlamb(i,:) = data%x
    End do

    data%A=A                                                    !Start to compute the last two eigenvalues
    y=det(data)
    y=data%detA
    Do i=1,n-2
       y=y/lambda(i)
    end do
    
    z=trace(data)
    z=data%TrA
    Do i=1,n-2
       z = z - lambda(i)
    end do
  
    !We have to work on the following relations:
    !lambda(n-1) + lambda (n) = z
    !lambda(n-1) * lambda (n) = y
    
    lambda(n)=(z+sqrt(z*z - 4._rk*y))/2._rk
    lambda(n-1)= z - lambda(n)
  
    u(n-1)=lambda(n-1)-1.E-2_rk                                   !Start to compute the last two eigenvectors
    data%A=A-u(n-1)*Uni
    If(inv(Data)/=0) then
      print*, 'Error inverting the matrix associated to the eigenvalue!'     
      stat = 3
      return
    end if
    data%A=data%iA
    If ((eigen(Data,para)/=0).or.(abs((1._rk/data%lambda)+u(n-1)-lambda(n-1))>1.E1_rk*para%eps_lambda)) then
      print*, 'Eigenvalue n-1 failed!'    
      stat = 3
      return
    end if
    vlamb(n-1,:) = data%x

    u(n)=lambda(n)-1.E-2_rk
    data%A=A-u(n)*Uni
    If(inv(Data)/=0) then
      print*, 'Error inverting the matrix associated to the eigenvalue!'     
      stat = 3
      return
    end if
    data%A=data%iA
    If ((eigen(Data,para)/=0).or.(abs(((1._rk/data%lambda)+u(n))-lambda(n))>10._rk*para%eps_lambda)) then
      print*, 'Eigenvalue n failed!'    
      stat = 3
      return
    end if
    vlamb(n,:) = data%x
    
    data%A = A                                                    !Reset data%A          
    stat = 0
  
  End function diagonalizationR

End Module library 
