      implicit real*8 (a-h,p-z)

      parameter (n=100)
      dimension uold(-n:n),unew(-n:n)
      dimension wold(-n:n),wnew(-n:n),wint(-n:n)
      f(u) = v*u
      dt = 1.d-3
      dx = 1.d-1
      v=1.d0
      alpha = v*dt/dx
      alp2 = dt*v/dx*0.5
      alp1 = dt*v/dx
      ninner=100
      coef = dt/dx
      
! initial state
      do i=-n,n
         x=i*dx
         uold(i)=dexp(-x*x/2.d0)*dcos(x*6.3d0)
         wold(i) = uold(i)
      enddo
      do itime = 1,100
         do inner = 1,ninner
            do i=-n+1,n-1
! lax
               unew(i) = 0.5d0*(uold(i+1)+uold(i-1))+
     1                   0.5d0*coef*(f(uold(i+1))-f(uold(i-1)))
            enddo
            do i=-n,n
               uold(i)=unew(i)
            enddo
!lw
            do i=-n,n-1
               wint(i)=0.5d0*(wold(i+1)+wold(i))+
     1                 0.5d0*coef*(f(wold(i+1))-f(wold(i)))
            enddo
            do i=-n+1,n-1
               wold(i) = wold(i)+coef*(f(wint(i))-f(wint(i-1)))
            enddo
            
         enddo
         do i=-n,n
            write(10,*) i*dx,itime*dt*ninner,uold(i),wold(i)
            write(11,*) i*dx,itime*dt*ninner,wold(i)
         enddo
         write(10,*)
         write(11,*)
      enddo
      stop
      end
