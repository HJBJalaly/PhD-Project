home
clear
close all

% parameters of robot
L=1;       % Link Lentgh
NoiseSTD=0.05;

% EF motion: Circle
f=1;
A=1;
phi=pi;
Tres=0.0025;
time=0:Tres:1*f;
Xef=A*cos(2*pi/f*time)+1.5;
Yef=A*sin(2*pi/f*time)+1;
DXef=-(2*pi/f)*A*sin(2*pi/f*time);
DYef= (2*pi/f)*A*cos(2*pi/f*time);

% Initial Point
q1=deg2rad(0);
q2=deg2rad(8.031);
q3=deg2rad(51.317);


for i=1:length(time)-1
   
    XPn= [Xef(i+1),Yef(i+1)]';
    XPo= [Xef(i),Yef(i)]';
    XPo= L*[cos(q1(i))+cos(q1(i)+q2(i))+cos(q1(i)+q2(i)+q3(i));
           sin(q1(i))+sin(q1(i)+q2(i))+sin(q1(i)+q2(i)+q3(i))];
    
    Qo =[q1(i) , q2(i) , q3(i)]';
    
    DeltaX =XPn- XPo;
    
    Jac= L*[ -sin(q1(i))-sin(q1(i)+q2(i))-sin(q1(i)+q2(i)+q3(i)) ,  -sin(q1(i)+q2(i))-sin(q1(i)+q2(i)+q3(i)) ,-sin(q1(i)+q2(i)+q3(i)) ;
             +cos(q1(i))+cos(q1(i)+q2(i))+cos(q1(i)+q2(i)+q3(i)) ,  cos(q1(i)+q2(i))+cos(q1(i)+q2(i)+q3(i)) ,cos(q1(i)+q2(i)+q3(i))];
    
    JacSharp= Jac'*(Jac*Jac')^-1;
         
    DeltaQ= JacSharp*DeltaX;
    
    Qn= Qo+DeltaQ+(NoiseSTD)^2*randn(3,1);
    q1(i+1)=Qn(1);
    q2(i+1)=Qn(2);
    q3(i+1)=Qn(3);
    
end

h1=figure;

for i=1:length(time)-1
    plot(Xef,Yef,'linewidth',2,'linestyle','-.')
    axis equal
    hold on
    
    Link1PosF =[0,0];
    Link1PosE =[(L)*cos(q1(i))   , (L)*sin(q1(i))];

    Link2PosF = Link1PosE;
    Link2PosE =Link1PosE+[(L)*cos(q1(i)+q2(i))  , (L)*sin(q1(i)+q2(i))];

    Link3PosE = Link2PosE;
    Link3PosF =Link2PosE+[(L)*cos(q1(i)+q2(i)+q3(i))  , (L)*sin(q1(i)+q2(i)+q3(i))];


    plot([Link1PosF(1) Link1PosE(1)], [Link1PosF(2) Link1PosE(2)],'linewidth',4,'color','r')
    plot([Link2PosF(1) Link2PosE(1)], [Link2PosF(2) Link2PosE(2)],'linewidth',4,'color','g')
    plot([Link3PosF(1) Link3PosE(1)], [Link3PosF(2) Link3PosE(2)],'linewidth',4,'color','m')

    Xdesin= L*[cos(q1(1:i))+cos(q1(1:i)+q2(1:i))+cos(q1(1:i)+q2(1:i)+q3(1:i));
               sin(q1(1:i))+sin(q1(1:i)+q2(1:i))+sin(q1(1:i)+q2(1:i)+q3(1:i))];
       
    plot(Xdesin(1,:),Xdesin(2,:),'linewidth',2)

    drawnow;
    pause(.2)
    hold off
end



%%

clear all;clc
load Data
s1=1;s2=201;
t=X(s1:s2,1);f=Y(s1:s2,1);

XX=[ones(s2-s1+1,1) t t.^2 t.^3 t.^4 t.^5];
T=(XX'*XX)\(XX'*f);
hold on
plot(t,f);grid on
eig(XX'*XX)
plot(t,T'*XX','r')
hold off
