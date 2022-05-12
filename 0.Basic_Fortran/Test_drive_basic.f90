!! Drive per testare le subroutine della libreria fino alla ricerca degli zeri

Program TestDrive

  Use library

  Implicit none
  Real(rk), parameter :: a=-1._rk, b=1._rk
  Integer(4), parameter :: N = 1000, o = 2
  Real(rk) :: delta, dy0, dy1, x(0:N), y(0:N), integ, zero
  Integer(4) :: i
  
  delta = (b-a)/N
  
  do i=0,N
    x(i) = a + i*delta
    y(i) = fun(x(i))
  end do
  
  open(1, file ='function.dat')
  
  do i=0,N
    if(differentiate(x, y, i, dy0, o_=o, m_=3)/=0) stop 'Errore funzione discreta'
    if(diff_fun(fun, x(i), dy1, o_=o)/=0) stop 'Errore funzione analitica'
    write(1,*) x(i), y(i), dfun(x(i)), dy0, dy1
  end do
    
    write(*,*) 'discrete int', integrate(x, y, integ), integ
    write(*,*) 'analitic int', int_fun(fun, x(0), x(N), y(0), y(N), integ), integ
    write(*,*) 'bisect stat', bisection(fun, x(2), zero), zero
    write(*,*) 'Secant stat', secant(fun, x(2), zero), zero
    write(*,*) 'Newton stat', newton(fun,dfun,x(2),zero), zero
    
    
  close(1)
  
  
  
Contains
  
  Function fun(x)
    Implicit none
    Real(rk), parameter :: xm=0.5_rk
    Real(rk), Intent(in) :: x
    Real(rk) :: fun
    
    fun = (x + xm)*exp(x**2)
  
  End function
  
  Function dfun(x)
    Implicit none
    Real(rk), parameter :: xm=0.5_rk
    Real(rk), Intent(in) :: x
    Real(rk) :: dfun
    
    if(diff_fun(fun, x, dfun, o_=1)/=0) stop 'Errore analitica'
  
  End function

End program TestDrive
