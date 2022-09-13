Program Drive_ode

  Use library

  Implicit None
  Real(rk), Parameter :: g=9.81_rk,l=2._rk, errore=1.E-12_rk
  Real(rk), Parameter :: theta0=pi*3/2._rk,omega0=0._rk,t0=0._rk,dt=5.E-3_rk,tfin=30._rk,b=0.5_rk

  Real(rk) :: theta(2)
  Type(ode_data) :: data
  Type(ode_para) :: para

  Allocate(data%y0(2),data%y(2))
  data%y0=[theta0,omega0]
  data%x0=t0
  data%dx=dt
  data%fun => pendulum_f
  data%dfun_t => pendulum_t
  data%dfun_y => pendulum_y
  para%eps=errore
  para%max_steps=1000
  theta=data%y0
  
  Open(1,file='OdeRunge.dat')
  Write(*,'(E19.12,3(x,E19.12))') data%x0,theta(1),-cos(theta(1)),sin(theta(1))
  Write(1,'(E19.12,3(x,E19.12))') data%x0,theta(1),-cos(theta(1)),sin(theta(1))
  
  Do While (data%x0<=tfin)
     data%y0=theta
     if (Runge_Kutta(data)/=0) Stop 'Runge-Kutta failed'
     theta=data%y
     data%x0=data%x0+data%dx
     Write(*,'(E19.12,3(x,E19.12))') data%x0,theta(1),-cos(theta(1)),sin(theta(1))
     Write(1,'(E19.12,3(x,E19.12))') data%x0,theta(1),-cos(theta(1)),sin(theta(1))
  EndDo
  
  Close(1)
  
  data%y0=[theta0,omega0]
  data%x0=t0
  theta=data%y0
  
  Open(2,file='OdePicard.dat')
  Write(*,'(E19.12,3(x,E19.12))') data%x0,theta(1),-cos(theta(1)),sin(theta(1))
  Write(2,'(E19.12,3(x,E19.12))') data%x0,theta(1),-cos(theta(1)),sin(theta(1))
  
  Do While (data%x0<=tfin)
     data%y0=theta
     if (Taylor3(data)/=0) stop 'taylor failed'
     if (IterPicard(data,para)/=0) stop 'iterPicard failed'
     theta=data%y
     data%x0=data%x0+data%dx
     Write(*,'(E19.12,3(x,E19.12))') data%x0,theta(1),-cos(theta(1)),sin(theta(1))
     Write(2,'(E19.12,3(x,E19.12))') data%x0,theta(1),-cos(theta(1)),sin(theta(1))
  EndDo
  
  Close(2)
  
  data%y0=[theta0,omega0]
  data%x0=t0
  theta=data%y0
  
  Open(3,file='OdeBasic.dat')
  Write(*,'(E19.12,3(x,E19.12))') data%x0,theta(1),-cos(theta(1)),sin(theta(1))
  Write(3,'(E19.12,3(x,E19.12))') data%x0,theta(1),-cos(theta(1)),sin(theta(1))
  
  Do While (data%x0<=tfin)
     data%y0=theta
     if (euler(data)/=0) stop 'euler failed'
     theta=data%y
     data%x0=data%x0+data%dx
     Write(*,'(E19.12,3(x,E19.12))') data%x0,theta(1),-cos(theta(1)),sin(theta(1))
     Write(3,'(E19.12,3(x,E19.12))') data%x0,theta(1),-cos(theta(1)),sin(theta(1))
  EndDo
  
  Close(3)

Contains

  Function pendulum_lin(theta,t)
    Implicit None
    Real(rk), Intent(in) :: theta(1:),t
    Real(rk) :: pendulum_lin(1:Size(theta))

    pendulum_lin(1)=theta(2)
    pendulum_lin(2)=-g/l*theta(1)
  End Function pendulum_lin

  Function pendulum(theta,t)
    Implicit None
    Real(rk), Intent(in) :: theta(1:),t
    Real(rk) :: pendulum(1:Size(theta))

    pendulum(1)=theta(2)
    pendulum(2)=-g/l*Sin(theta(1))
  End Function pendulum

  Function pendulum_f(theta,t)
    Implicit None
    Real(rk), Intent(in) :: theta(1:),t
    Real(rk) :: pendulum_f(1:Size(theta))

    pendulum_f(1)=theta(2)
    pendulum_f(2)=-g/l*Sin(theta(1))-b*theta(2)+b*sin(t)
  End Function pendulum_f

    Function pendulum_t(theta,t)
    Implicit None
    Real(rk), Intent(in) :: theta(1:),t
    Real(rk) :: pendulum_t(1:Size(theta))

    pendulum_t(1)=0._rk
    pendulum_t(2)=b*cos(t)
  End Function pendulum_t

  Function pendulum_y(theta,t)
    Implicit None
    Real(rk), Intent(in) :: theta(1:),t
    Real(rk) :: pendulum_y(1:Size(theta),1:Size(theta))

    pendulum_y(1,1)=0._rk
    pendulum_y(1,2)=1._rk
    pendulum_y(2,1)=-g/l*cos(theta(1))
    pendulum_y(2,2)=-b
  End Function pendulum_y

End Program Drive_ode
