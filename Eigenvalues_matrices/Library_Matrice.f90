Module LibraryMatr

  Use Library

Contains
!*********************************************************************************************************************************************
   Function tridiagsup(Data) Result (stat)
   ! stat = 0 success
   ! stat = 1 matrix A is not square
   ! stat = 2 proportional rows
   ! stat = 3 more than one diagonal element is equal to zero
   ! stat = 4 b is incompatible with A
   ! stat = 5 Bm is incompatible with A

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

!*********************************************************************************************************************************************
  Function linsys(Data) Result (stat) 
  ! stat = 0 success
  ! stat = 1 matrix A is not square
  ! stat = 2 proportional rows
  ! stat = 3 more than one diagonal element is equal to zero
  ! stat = 4 b is incompatible with A
  ! stat = 5 b is not allocated or incompatible with A
  ! stat = 6 data%At(data%v(n),n) = 0

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

    If (.Not.Allocated(Data%b).Or.Size(Data%b)/=n) Then     !Teniamo problemi con b
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

!*********************************************************************************************************************************************
  Function det(Data) Result (stat)
  ! stat = 0 success
  ! stat = 1 matrix A is not square
  ! stat = 2 proportional rows
  ! stat = 3 more than one diagonal element is equal to zero

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

!*********************************************************************************************************************************************
  Function inv(Data) Result (stat)
  ! stat = 0 success
  ! stat = 1 matrix A is not square
  ! stat = 2 proportional rows
  ! stat = 3 more than one diagonal element is equal to zero

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

!*********************************************************************************************************************************************
  Function eigen(data, para) result (stat)
  ! stat = 0 success
  ! stat = 1 matrix A is not square
  ! stat = 2 proportional rows
  ! stat = 3 more than one diagonal element is equal to zero
  ! stat = 4 superato max steps
    
    Implicit none
    Type(matrix):: data
    Type(matrix_para) ::para

    Real (rk)::lambda,c(1:size(data%A,1)),v(1:size(data%A,1))
    Integer(4):: stat

    if (allocated(data%x)) deallocate(data%x)
    allocate(data%x(1:size(data%A,1)))
    call random_seed()                                                          !Generazione casuale di num tra 0 e 1
    call random_number(data%x)
    data%x=data%x/sqrt(dot_product(data%x,data%x))
    data%lambda=0._rk
    para%steps=0

    Do while (para%steps==0.or.(Any(abs(data%x-v)>para%eps_x).and.Any(abs(data%x+v)>para%eps_x))&
         .or.(data%lambda-lambda)>para%eps_lambda)                                  !Controllo sulla precisione di auto-vett/val
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
  
!*********************************************************************************************************************************************
   Function Gershgorin(data,u) result (stat)
    
   Implicit none
   Type(matrix):: data
   Type(matrix_para) ::para
   Real (rk)::u(3:size(data%A,1))

   Real (rk)::v(1:size(data%A,1),1:2),a,b,c(1:2)
   Integer(4):: stat,i,j,n

   n=(size(data%A,1))

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

   Do j=n-1,1,-1                   !Dispone i centri e i raggi in ordine decrescente(dei centri)
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
         Print '(x,A,x,(x,E13.6),x,A,x,(x,E13.6),x,A)', 'ATTENZIONE: I cerchi di centri',v(i,1),'e',v(i+1,1), 'si intersecano.'
       End if
   End do
   Print *, '   Centro e raggio del cerchio'
   Do i=1,n
       Print '(2(x,E13.6))', v(i,:)
   End do

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
   
   stat=0

  End Function Gershgorin
  
!*********************************************************************************************************************************************
  
End Module LibraryMatr
