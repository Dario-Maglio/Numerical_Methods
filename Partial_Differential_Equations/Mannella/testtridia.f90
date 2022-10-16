      implicit real*8(a-h,o-z)
      parameter (n=10)
      dimension dia(n),dup(n),dlo(n),x(n),b(n)

      do i=1,n
         dia(i) = dcos(1.d0*i)
         dup(i) = dsin(2.d0*i)
         dlo(i) = dsin(3.d0*i)
         x(i) = dcos(2.d0*i)
      enddo
      b(1) = dia(1)*x(1)+dup(1)*x(2)
      do i=2,n-1
         b(i) = dlo(i-1)*x(i-1)+dia(i)*x(i)+dup(i)*x(i+1)
      enddo
      b(n) = dlo(n-1)*x(n-1)+dia(n)*x(n)
      print *,x
      print *,b
      x=0.d0
      call tridia(n,dia,dup,dlo,b,x)
      print *,x
      print *,b
      stop
      end
      
