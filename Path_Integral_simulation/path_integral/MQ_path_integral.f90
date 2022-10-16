Program MQ_path_integral

   Use flibrary

   Implicit none
   Integer(4), parameter:: measures=1500, i_decorrel=5, i_termal=100, d_metrop=1
   Integer(4), parameter:: iflag=1, d_flag=1, side_lattice=10
   Real(rk):: eta=0.1
   
   Type(lattice_data) :: data
   Integer(4):: iter, idec
   
   data%d_flag=d_flag
   data%side_latt = side_lattice
   data%i_flag = iflag
   
   data%measures = measures
   data%i_decorrel = i_decorrel
   data%i_termal = i_termal
   data%d_metrop = d_metrop
   data%eta = eta
   
   call gen_lattice(data)
   
   call ranstart()
   
   !Print to test config
   do iter=1, data%vol_latt
      if (modulo(i-1, data%side_latt)==0) print *,'---'
      print '(E15.5)', data%field(i)
   end do
   print *,'---'
   !Print for test nearest
   do iter=1,data%vol_latt
      print *, data%nn(i,:)
   end do
   
   !Termalizzazione metropolis
   do iter = 1, i_termal
      call update_metropolis(data)
   enddo
   
   !Raccolta dati
   do iter = 1, measures
      do idec = 1, i_decorrel
         call update_metropolis(data)
      enddo
      call measure(data)
   enddo
   
   
   
   call ranfinish()
   
   open(1,file='lattice',status='unknown')
   write(1,*) data%field
   close(1)

End program MQ_path_integral
