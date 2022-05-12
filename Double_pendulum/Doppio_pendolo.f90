Program drive_doppio_pendolo

  Use library

  Implicit None
  Real(rk), Parameter::g=9.81_rk,l=0.5_rk,m1=2._rk,m2=3.5_rk,b=0.5_rk
  Real(rk), Parameter::theta0=pi/4._rk,phi0=-pi/3._rk,omega1_0=0._rk,omega2_0=0._rk,t0=0._rk,dt=5.E-5_rk,tfin=30._rk
  Real(rk)::theta_phi(4)
  Type(ode_data)::data
  

  Allocate(data%y0(4),data%y(4))
 
  data%y0=[theta0,omega1_0,phi0,omega2_0]
  data%x0=t0
  data%dx=dt
  data%fun=>pendolo
  theta_phi=data%y0
  steps=0
  Open(1,file='doppio_pendolo.dat')
  Write(*,'(E19.12,2(x,E19.12))') data%x0,theta_phi(1),theta_phi(3),&
             l*Sin(theta_phi(1)),-l*Cos(theta_phi(1)),l*(Sin(theta_phi(1))&
             +Sin(theta_phi(3))),&
             -l*(Cos(theta_phi(1))+Cos(theta_phi(3)))
  Write(1,'(E19.12,2(x,E19.12))') data%x0,theta_phi(1),theta_phi(3),&
             l*Sin(theta_phi(1)),-l*Cos(theta_phi(1)),l*(Sin(theta_phi(1))&
             +Sin(theta_phi(3))),&
             -l*(Cos(theta_phi(1))+Cos(theta_phi(3)))
 
  Do While(data%x0<=tfin)
    data%y0=theta_phi
    If(Runge_Kutta(data)/=0) Stop 'Runge-Kutta failed'
    data%x0=data%x0+data%dx
    theta_phi=data%y
     
    Do while(theta_phi(1)>=2._rk*pi)
      theta_phi(1)=theta_phi(1)-2._rk*pi
    End Do
    Do while(theta_phi(3)>=2._rk*pi)
      theta_phi(3)=theta_phi(3)-2._rk*pi
    End Do
        
    Write(*,'(E19.12,6(x,E19.12))') data%x0,theta_phi(1),theta_phi(3),&
         l*Sin(theta_phi(1)),-l*Cos(theta_phi(1)),l*(Sin(theta_phi(1))&
         +Sin(theta_phi(3))),-l*(Cos(theta_phi(1))+Cos(theta_phi(3)))
    Write(1,'(E19.12,6(x,E19.12))') data%x0,theta_phi(1),theta_phi(3),&
             l*Sin(theta_phi(1)),-l*Cos(theta_phi(1)),l*(Sin(theta_phi(1))&
             +Sin(theta_phi(3))),-l*(Cos(theta_phi(1))+Cos(theta_phi(3)))
             
   End Do
  
   Close(1)
   
   
 Contains

   Function pendolo(theta,t)
     Implicit None
     Real(rk),Intent(In)::theta(1:),t
     Real(rk)::pendolo(1:Size(theta))
     
     pendolo(1)=theta(2)
     pendolo(2)=-b/(m1+m2)*(2._rk*theta(2)+Cos(theta(1)-theta(3))*theta(4))&
           -m2/(m1+m2)*(theta(4)**2-theta(2)*theta(4))*Sin(theta(1)-theta(3))&
           -g/l*Sin(theta(1))-m2/(m1+m2)*Cos(theta(1)-theta(3))&
           *((-b/m2*(theta(4)+Cos(theta(1)-theta(3))*theta(2))&
          -(theta(2)*theta(4)-theta(2)**2)*Sin(theta(1)-theta(3))&
          -g/l*Sin(theta(3))-Cos(theta(1)-theta(3))*(-b/(m1+m2)&
          *(2._rk*theta(2)+Cos(theta(1)-theta(3))*theta(4))&
          -m2/(m1+m2)*(theta(4)**2-theta(2)*theta(4))*Sin(theta(1)-theta(3))&
          -g/l*Sin(theta(1))))/(1._rk-m2/(m1+m2)*(Cos(theta(1)-theta(3)))**2))
        
     pendolo(3)=theta(4)
     pendolo(4)=(-b/m2*(theta(4)+Cos(theta(1)-theta(3))*theta(2))&
          -(theta(2)*theta(4)-theta(2)**2)*Sin(theta(1)-theta(3))&
          -g/l*Sin(theta(3))-Cos(theta(1)-theta(3))*(-b/(m1+m2)&
          *(2._rk*theta(2)+Cos(theta(1)-theta(3))*theta(4))&
          -m2/(m1+m2)*(theta(4)**2-theta(2)*theta(4))*Sin(theta(1)-theta(3))&
          -g/l*Sin(theta(1))))/(1._rk-m2/(m1+m2)*(Cos(theta(1)-theta(3)))**2)
   End Function pendolo

 End Program drive_doppio_pendolo
 
