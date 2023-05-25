      subroutine tridia(n,dia,dup,dlo,b,x)
      implicit real*8 (a-h,o-z)
      dimension dia(n),dup(n),dlo(n),b(n),x(n)

! It solves the tridiagonal problem
!
!   dia  dup    0     0     0
!   dlo  dia   dup    0     0
!    0   dlo   dia   dup    0       X   =   b
!    0    0    dlo   dia   dup
!    0    0     0    dlo   dia
! 
! note  dlo(1) = dlo_{21}
!       dia(1) = dia_{11}
!       dup(1) = dup_{12}
! First, go down, use dia to remove dlo

      do i = 2,n
         if(dia(i-1).eq.0.d0) print *,"warning, I cannot compute the inverse!"
         factor = dlo(i-1)/dia(i-1)
!         dlo(i-1) = 0.d0
         dia(i)   = dia(i)-factor*dup(i-1)
         b(i)     = b(i) - factor*b(i-1)
      enddo

! now, ready to go up
      if(dia(n).eq.0.d0) print *,"warning, I cannot compute the inverse!"
      x(n) = b(n)/dia(n)
      do i=n-1,1,-1
         b(i) = b(i)-dup(i)*x(i+1)
         if(dia(i).eq.0.d0) print *,"warning, I cannot compute the inverse!"
         x(i) = b(i)/dia(i)
      enddo

      return
      end  
