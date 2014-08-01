function [Torque,Q,D1Q,D2Q,IntU2,IntUdq,IntAbsUdq,Error]=ShowTime(Coef,time,Tres,Degree,Xef,Yef,m,L,g,ShowFlag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% close all

%% Regerate Trjectory from Coef
CoefP_q1=Coef(1:(Degree+1));
CoefP_q2=Coef((Degree+1)+1:2*(Degree+1));
CoefP_q3=Coef(2*(Degree+1)+1:3*(Degree+1));
D1CoefP_q1=CoefP_q1(1:end-1).*(Degree:-1:1);
D1CoefP_q2=CoefP_q2(1:end-1).*(Degree:-1:1);
D1CoefP_q3=CoefP_q3(1:end-1).*(Degree:-1:1);
D2CoefP_q1=D1CoefP_q1(1:end-1).*(Degree-1:-1:1);
D2CoefP_q2=D1CoefP_q2(1:end-1).*(Degree-1:-1:1);
D2CoefP_q3=D1CoefP_q3(1:end-1).*(Degree-1:-1:1);

Q1=polyval(CoefP_q1,time);
Q2=polyval(CoefP_q2,time);
Q3=polyval(CoefP_q3,time);
D1Q1=polyval(D1CoefP_q1,time);
D1Q2=polyval(D1CoefP_q2,time);
D1Q3=polyval(D1CoefP_q3,time);
D2Q1=polyval(D2CoefP_q1,time);
D2Q2=polyval(D2CoefP_q2,time);
D2Q3=polyval(D2CoefP_q3,time);


%% EF

% % For Circle
% Xef=A*cos(2*pi/f*time)+1.5;
% Yef=A*sin(2*pi/f*time)+1;

% % For Line
% Xef=A*cos(2*pi/(1*f)*time);
% Yef=2*ones(size(time));

Pos=[Xef;Yef];

RPos=L*[cos(Q1)+cos(Q1+Q2)+cos(Q1+Q2+Q3);
        sin(Q1)+sin(Q1+Q2)+sin(Q1+Q2+Q3)];
Error=sum(sum((RPos-Pos).^2))*Tres;

if(ShowFlag)
    
    figure('name','WorkSapce')
        plot(Pos(1,:),Pos(2,:),'linewidth',2,'linestyle','-','color','b')
        title('WorkSpace','FontWeight','bold')
        hold on
        plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','-','color','r')
        hold off
        xlabel('x (m)')
        ylabel('y (m)')

    
    
    Y=[Q1',Q2',Q3'];
    AnimBot3DOF(time,Y,L);
end
%% Trajectory
if(ShowFlag)
    
    figure('name','compare trajectory')
        subplot(3,1,1)
        plot(time,Q1)
        title('Jonits Trajectory','FontWeight','bold')
        grid on
        xlabel('time (s)')
        ylabel('q_1 (rad)')

        subplot(3,1,2)
        plot(time,Q2)
        grid on
        xlabel('time (s)')
        ylabel('q_2 (rad)')

        subplot(3,1,3)
        plot(time,Q3)
        grid on
        xlabel('time (s)')
        ylabel('q_3 (rad)')

end
%% Torque

mL1=m;
mL2=m;
mL3=m;
LL1=L;
LL2=L;
LL3=L;

Torque=zeros(3,length(time));
IntU2=0;
IntAbsUdq=0;
IntUdq=0;

for i=1:length(time)
    q1=Q1(i);
    q2=Q2(i);
    q3=Q3(i);
    D1q1=D1Q1(i);
    D1q2=D1Q2(i);
    D1q3=D1Q3(i);
    D2q1=D2Q1(i);
    D2q2=D2Q2(i);
    D2q3=D2Q3(i);

    
    MM=[mL3*LL1*LL3*cos(q2+q3)+LL1*LL2*cos(q2)*mL2+2*LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL3*LL1^2+LL1^2*mL2+mL1*LL1^2/3+mL2*LL2^2/3,mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6;
        mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6;
        mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6,mL3*LL3^2/3];

    CC= [-mL3*LL1*LL3*sin(q2+q3)*(D1q2+D1q3)/2-LL2*(LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL1*LL3*(D1q1+D1q2+D1q3)*sin(q2+q3)/2-LL2*(sin(q2)*D1q1*LL1*(mL2+2*mL3)+LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL3*(LL1*sin(q2+q3)+LL2*sin(q3))*(D1q1+D1q2+D1q3)/2;
          mL3*LL1*D1q1*LL3*sin(q2+q3)/2+(sin(q2)*D1q1*LL1*(mL2+2*mL3)-mL3*LL3*sin(q3)*D1q3)*LL2/2,-mL3*LL2*LL3*sin(q3)*D1q3/2,-mL3*LL3*sin(q3)*(D1q1+D1q2+D1q3)*LL2/2;
          mL3*(LL1*D1q1*sin(q2+q3)+LL2*sin(q3)*(D1q1+D1q2))*LL3/2,LL3*LL2*sin(q3)*(D1q1+D1q2)*mL3/2,0];

    GG= [g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2+2*cos(q1)*LL1*mL3+2*cos(q1)*LL1*mL2+cos(q1)*LL1*mL1)/2;
         g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2)/2;
         g*mL3*LL3*cos(q1+q2+q3)/2];


    Torque(:,i) = MM*[D2q1;D2q2;D2q3] + CC*[D1q1;D1q2;D1q3] + GG;
%     IntU2=IntU2+Torque(:,i)'*Torque(:,i)*Tres;
%     IntAbsUdq=IntAbsUdq+abs(Torque(:,i)'*[D1q1;D1q2;D1q3])*Tres;
%     IntUdq=IntUdq+(Torque(:,i)'*[D1q1;D1q2;D1q3])*Tres;
end

IntU2=sum(sum(Torque.^2))*Tres;
IntAbsUdq=sum(sum(abs(Torque.*[D1Q1;D1Q2;D1Q3])))*Tres;
IntUdq=sum(abs(sum((Torque.*[D1Q1;D1Q2;D1Q3]),2)))*Tres;

if(ShowFlag)
    
    figure('name','Torque vs Time')
        plot(time,Torque)
        title('Torque','FontWeight','bold')
        legend('\tau_1','\tau_2','\tau_3')
        xlabel('time (s)')
        ylabel('\tau')
        grid on

    figure('name','Torque vs Angle')
        subplot(3,1,1)
        plot(Q1,Torque(1,:))
        title('Torque Angle Profile','FontWeight','bold')
        hold on
        plot(Q1(1),Torque(1,1),'linewidth',2,'linestyle','none','marker','*','markersize',6)
        xlabel('q_1 (rad)')
        ylabel('\tau_1')
        hold off
        grid on

        subplot(3,1,2)
        plot(Q2,Torque(2,:))
        hold on
        plot(Q2(1),Torque(2,1),'linewidth',2,'linestyle','none','marker','*','markersize',6)
        xlabel('q_2 (rad)')
        ylabel('\tau_2')
        hold off
        grid on

        subplot(3,1,3)
        plot(Q3,Torque(3,:))
        hold on
        plot(Q3(1),Torque(3,1),'linewidth',2,'linestyle','none','marker','*','markersize',6)
        xlabel('q_3 (rad)')
        ylabel('\tau_3')
        hold off
        grid on

end


Q=[Q1;Q2;Q3];
D1Q=[D1Q1;D1Q2;D1Q3];
D2Q=[D2Q1;D2Q2;D2Q3];
end

