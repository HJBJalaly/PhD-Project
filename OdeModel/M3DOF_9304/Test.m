home
clear
close all

% envirment
% parameters of robot
m=1;
L=1;
% EF motion
f=1;
A=1;
phi=pi;
Tres=0.005;
time=0:Tres:1*f;

% Circle motion
Xef=A*cos(2*pi/f*time)+1.5;
Yef=A*sin(2*pi/f*time)+1;
DXef=-(2*pi/f)*A*sin(2*pi/f*time);
DYef= (2*pi/f)*A*cos(2*pi/f*time);


q1=deg2rad(0);
q2=deg2rad(8.031);
q3=deg2rad(51.317);

h1=figure;
plot(Xef,Yef)
axis equal
figure
plot(time,Yef)
hold all
plot(time,Xef)



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
    
    Qn= Qo+DeltaQ;
    q1(i+1)=Qn(1);
    q2(i+1)=Qn(2);
    q3(i+1)=Qn(3);
    
end

figure()
subplot(3,1,1)
plot(time,q1)
subplot(3,1,2)
plot(time,q2)
subplot(3,1,3)
plot(time,q3)


Xdesin= L*[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
           sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];
       
figure(h1)
hold all
plot(Xdesin(1,:),Xdesin(2,:))




% pos=[xef;yef];
% Dpos=[Dxef;Dyef];
% D2pos=[D2xef;D2yef];
% 
% 
% % Joint Trajectories
% 
% for tt=1:length(time)-1
%     lPos=L*[cos(q1(tt))+cos(q1(tt)+q2(tt))+cos(q1(tt)+q2(tt)+q3(tt));
%             sin(q1(tt))+sin(q1(tt)+q2(tt))+sin(q1(tt)+q2(tt)+q3(tt))];
% 
%     cPos=pos(:,tt);
%     dx=cPos-lPos;
% 
%     JJ=L*[-sin(q1(tt))-sin(q1(tt)+q2(tt))-sin(q1(tt)+q2(tt)+q3(tt)), -sin(q1(tt)+q2(tt))-sin(q1(tt)+q2(tt)+q3(tt)), -sin(q1(tt)+q2(tt)+q3(tt));
%            cos(q1(tt))+cos(q1(tt)+q2(tt))+cos(q1(tt)+q2(tt)+q3(tt)),  cos(q1(tt)+q2(tt))+cos(q1(tt)+q2(tt)+q3(tt)),  cos(q1(tt)+q2(tt)+q3(tt))]; 
% 
%     dq=JJ'*(JJ*JJ')^-1*dx;
% 
%     q1(tt+1)=q1(tt)+dq(1);
%     q2(tt+1)=q2(tt)+dq(2);
%     q3(tt+1)=q3(tt)+dq(3);
%     
% end
% 
% RPosVal=L*[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
%         sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];
% 
%     
% % Required torque for desired path
% q=[q1;q2;q3];
% Dq=differential(q,time,Tres);
% D2q=differential(Dq,time,Tres);
% torqueDesire=TorqueCalculator(D2q,Dq,q,g,m,m,m,L,L,L);
% 
% Middle=ceil(length(time)/2);
% 
% % show    
% close all
% figure('name','Path (3DoF)')
%     plot(xef,yef,'linewidth',2.5,'linestyle','-','color','g')
%     hold on 
%     plot(RPosVal(1,Middle:end),RPosVal(2,Middle:end),'linewidth',2,'linestyle','-.','color','b')
%     legend('Desired','Static Path Planing')
%     hold off
%     axis equal
% 
% figure('name','Joint Trajectories')
%     subplot(3,1,1)
%     plot(time(Middle:end),q1(Middle:end))
%     grid on
%     subplot(3,1,2)
%     plot(time(Middle:end),q2(Middle:end))
%     grid on
%     subplot(3,1,3)
%     plot(time(Middle:end),q3(Middle:end))
%     grid on
%     
% figure('name','Desired Torque vs Time')
%     plot(time(Middle:end),torqueDesire(:,Middle:end))
%     legend('\tau_1','\tau_2','\tau_3')
%     xlabel('time (s)')
%     ylabel('\tau')
%     grid on
%     
% 
% 
% figure('name','Desired Torque vs Angle')
%     subplot(3,1,1)
%     plot(q1(Middle:end),torqueDesire(1,Middle:end))
%     grid on
%     subplot(3,1,2)
%     plot(q2(Middle:end),torqueDesire(2,Middle:end))
%     grid on
%     subplot(3,1,3)
%     plot(q3(Middle:end),torqueDesire(3,Middle:end))
%     grid on
%     
% Y=[q1(1:end)',q2(1:end)',q3(1:end)'];
% AnimBot3DOF(time(1:end),Y,L);
