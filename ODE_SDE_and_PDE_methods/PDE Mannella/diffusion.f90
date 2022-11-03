      implicit real*8 (a-h,o-z)
      parameter (n=100)
      dimension uold(-n:n),unew(-n:n)
      dimension dia(-n:n),dlo(-n:n),dup(-n:n)
      dimension bnoto(-n:n),xsol(2*n+1)
      f(x) = dcos(6.3*x) + dcos(time)! -ak*(x-1)

      diff=1.d-1
      dx = 1.d-1
      dt = 3.d-3
      ninner = 1000
      ntime = 20
      alpha = diff*dt/(dx*dx)
      ak = 1.d0
      beta = dt/(2.d0*dx)

!! integrates u_t = D u_xx
      print *,dt,ninner
!! Explicit
!     initial state
      do i = -n,n
         x=i*dx
         uold(i) = dexp(-x*x/(2.d0*diff))
      enddo
      unew=0.d0
      itime = 0
      do i=-n,n
         write(20,*) i*dx,itime*dt*ninner,uold(i)
      enddo
      write(20,*)

      do itime = 1,ntime
         do inner = 1,ninner
            do i=-n+1,n-1
               time = dt*(inner+itime*ninner)
               x=i*dx
               unew(i) = uold(i)+alpha*(uold(i+1)-2*uold(i)+uold(i-1)) &
                      -beta*(f(x+dx)*uold(i+1)-f(x-dx)*uold(i-1))

            enddo
            uold = unew
         enddo
         do i=-n,n
            write(20,*) i*dx,itime*dt*ninner,uold(i)
         enddo
         write(20,*)
      enddo
      
!! implicit
!     initial state
      do i = -n,n
         x=i*dx
         uold(i) = dexp(-x*x/(2.d0*diff))
      enddo
      unew=0.d0
      itime = 0
      do i=-n,n
         write(21,*) i*dx,itime*dt*ninner,uold(i)
      enddo
      write(21,*)

      do itime = 1,ntime

         do inner = 1,ninner
               time = dt*(inner+itime*ninner)
! necessario solo se diffusione e drift dipendono da tempo e spazio
            do i=-n,n
               x=i*dx
               dia(i) = (1+2*alpha)
               dlo(i) = -alpha-beta*f(x-dx)
               dup(i) = -alpha+beta*f(x+dx)
            enddo
            call tridia(2*n-1,dia,dup,dlo(-n+1),uold(-n+1),unew(-n+1))
            uold = unew
         enddo
         do i=-n,n
            write(21,*) i*dx,itime*dt*ninner,uold(i)
         enddo
         write(21,*)
      enddo

!! KN
!     initial state
      do i = -n,n
         x=i*dx
         uold(i) = dexp(-x*x/(2.d0*diff))
      enddo
      unew=0.d0
      itime = 0
      do i=-n,n
         write(22,*) i*dx,itime*dt*ninner,uold(i)
      enddo
      write(22,*)

      do itime = 1,ntime

         do inner = 1,ninner
               time = dt*(inner+itime*ninner)
! necessario solo se diffusione e drift dipendono da tempo e spazio
            do i=-n,n
               x=i*dx
               dia(i) = (1+alpha)
               dlo(i) = -alpha*0.5d0-beta*f(x-dx)*0.5d0
               dup(i) = -alpha*0.5d0+beta*f(x+dx)*0.5d0
            enddo
            do i=-n+1,n-1
               x=dx*i
               bnoto(i) = uold(i)*(1-alpha)+ &
                      alpha*0.5d0*(uold(i-1)+uold(i+1)) + &
                      beta*0.5d0*(f(x-dx)*uold(i-1)-f(x+dx)*uold(i+1))
            enddo
            call tridia(2*n-1,dia,dup,dlo(-n+1),bnoto(-n+1),unew(-n+1))
            uold = unew
         enddo
         do i=-n,n
            write(22,*) i*dx,itime*dt*ninner,uold(i)
         enddo
         write(22,*)
      enddo
      print *,dt,ninner
      stop
      end
