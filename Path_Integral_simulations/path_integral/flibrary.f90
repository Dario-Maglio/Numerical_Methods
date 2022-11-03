Module flibrary

!!=================================================================
!!=================================================================
!!==	  General parameters and abstract interface            ====
!!=================================================================
!!=================================================================

  Implicit None 
  Integer(4), Parameter :: rk=8
  Real(rk), Parameter :: pi=Acos(-1._rk)
  
  Type :: lattice_data
  !!Features and structure of the n-dimensional lattice with side side_latt
      Integer(4):: measures, i_decorrel, i_termal, d_metrop
      Integer(4):: i_flag, d_flag, side_latt, vol_latt
      Procedure(fun_R_R), Nopass, pointer:: force
      Integer(4), allocatable:: nn(:,:)
      Real(rk), allocatable:: field(:)
      Real(rk):: eta
   End type

  Abstract Interface
     
     Function fun_R_R(x)
     !!Funzione reale di variabile reale
       Import :: rk
       Real(rk), Intent(in) :: x
       Real(rk) :: fun_R_R 
     End Function fun_R_R

  End Interface



Contains

!!=================================================================
!!=================================================================
!!==           Lattice operations and MC method                ====
!!=================================================================
!!=================================================================

   Subroutine gen_lattice(data)
   !!Generate the lattice starting from the parameters in ising_data with side side_latt.
   !!Allocate also the matrix nn with the indexes of the nearest neighbors for each point. 
   !!The d_flag indicates the dimensions and the geometry of the lattice.
   !!Start d_flag= 1 -> 1D chain, 2 -> 2D square
   !!Start i_flag= 0 -> cold, 1 -> hot, else -> the last saved configuration.
   
      Type(lattice_data) :: data
      Integer(4):: i, val
      
      !Allocate field and nn of the correct dimension
      select case(data%d_flag)
         case(1)
            data%vol_latt = data%side_latt
            allocate(data%nn(data%vol_latt, 2))
         case(2)
            data%vol_latt = data%side_latt**2
            allocate(data%nn(data%vol_latt, 4))
         case default
            print *, 'Conditions not implemented yet!'
            return
      end select
      allocate(data%field(data%vol_latt))
      
      !Initialize the lattice with the initial thermal condition
      select case(data%i_flag)
         case(0)                 
            forall (i=1:data%vol_latt) data%field(i) = 0._rk
         case(1)
            call random_number(data%field)
            forall (i=1:data%vol_latt) data%field(i) = real(1. - 2.*data%field(i), rk)
         case default
            open(1,file='lattice',status='old')
            read(1,*) data%field
            close(1)
      end select
      
      
      !Initialize the nearest neighbors matrix with periodic boundary conditions
      do i=1, data%vol_latt
         val = modulo(i, data%side_latt)
         data%nn(i,1) = i + 1
         data%nn(i,2) = i - 1
         if (val == 0) then
            data%nn(i,1) = data%nn(i,1) - data%side_latt
         else if (val == 1) then
            data%nn(i,2) = data%nn(i,2) + data%side_latt
         end if
         
         if (data%d_flag>1) then
            data%nn(i,3) = i - data%side_latt
            data%nn(i,4) = i + data%side_latt
            val = (i - 1) / data%side_latt
            if (val == 0) then
               data%nn(i,3) = data%vol_latt + data%nn(i,3)
            else if (val == (data%side_latt-1)) then
               data%nn(i,4) = data%nn(i,4) - data%vol_latt
            end if
         end if
      end do 
   
   End subroutine gen_lattice
   
   Subroutine update_metropolis(data)
   
      Type(lattice_data) :: data
   
      Integer(4):: i, j, sizenn
      Real(rk):: c1, c2, x, phi, phi_prova
   
      c1 = 1._rk/data%eta
      c2 = (1._rk/data%eta + data%eta/2._rk)
      sizenn = size(data%nn(i:))
      
      do i = 1, nlatt                     
         
         force = 0.
         do j=1, sizenn
            force = force + field(data%nn(i:j))
         end do
         
         phi =  field(i)
         
         phi_prova = phi + 2.*d_metro*(0.5-ran2())                  

         p_rat = c1 * phi_prova * force - c2 * phi_prova**2 - (c1 * phi * force - c2 * phi**2)
         
         x = log(ran2())                      
                                              
         if (x.lt.p_rat) data%field(i) = phi_prova
         
      end do                     
   
   End subroutine update_metropolis

   Subroutine measure()

      parameter (nlatt = 10)
      common/lattice/field(nlatt)
      common/move/npp(nlatt),nmm(nlatt)

      obs1 = 0.0
      obs2 = 0.0
      do i = 1,nlatt        
        obs1 = obs1 + field(i)**2
        obs2 = obs2 + (field(i)-field(npp(i)))**2
      enddo                             

      obs1 = obs1/float(nlatt) !! media sul singolo path di y^2 
      obs2 = obs2/float(nlatt) !! media sul singolo path di Delta y^2 
      !write(2,*) obs1,obs2

      return
   End

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
  
!******************************************************************
  
End module flibrary
