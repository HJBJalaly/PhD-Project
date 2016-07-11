function Main()



%% Initialization

clear
close all
home
rand('twister', sum(100*clock));


%% create a active sample motion with 3 active joints as initial input of optimization
home
clear
clc
close all
% envirment
g=9.81*1;
% parameters of robot
m=1;
L=1;
% EF motion
A=.75;
Tres=0.005;
f=.5;
phi=pi/2;

%%
% Circle motion
Name='Test';
f=0.5;
phi=pi/2;
Tres=0.005;
time=0:Tres:40/f;
Middle=20;%ceil(1*length(time)/40);
Start=20;
xef=A*cos(2*pi*f*time)+0;
yef=A*sin(2*pi*f*time)+2;
Dxef=-(2*pi*f)*A*sin(2*pi*f*time);
Dyef= (2*pi*f)*A*cos(2*pi*f*time);
D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time);
D2yef=-(2*pi*f)^2*A*sin(2*pi*f*time);

 
% Ellipose motion
Name='Ellipose';
f=0.5;
phi=pi/2;
Tres=0.005;
time=0:Tres:100/f;
Middle=ceil(2*length(time)/100);
Start=1;
xef=A*cos(2*pi*f*time)+0;
yef=A/2*sin(2*pi*f*time)+2;
Dxef=-(2*pi*f)*A*sin(2*pi*f*time);
Dyef= (2*pi*f)*A/2*cos(2*pi*f*time);
D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time);
D2yef=-(2*pi*f)^2*A/2*sin(2*pi*f*time);
% 
 
% % Rose motion
% Name='Rose';
% f=0.5;
% phi=pi/2;
% Tres=0.005;
% time=0:Tres:60/f;
% Middle=ceil(58*length(time)/60);
% Start=20;
% xef=A*cos(2*2*pi*f*time).*cos(2*pi*f*time)+0;
% yef=A*cos(2*2*pi*f*time).*sin(2*pi*f*time)+2;
% Dxef=- 2*A*pi*f*sin(2*pi*f*time).*cos(4*pi*f*time) - 4*A*pi*f*sin(4*pi*f*time).*cos(2*pi*f*time);
% Dyef= 2*A*pi*f*cos(2*pi*f*time).*cos(4*pi*f*time) - 4*A*pi*f*sin(2*pi*f*time).*sin(4*pi*f*time);
% D2xef=16*A*pi^2*f^2*sin(2*pi*f*time).*sin(4*pi*f*time) - 20*A*pi^2*f^2*cos(2*pi*f*time).*cos(4*pi*f*1*time);
% D2yef=- 20*A*pi^2*f^2*sin(2*pi*f*time).*cos(4*pi*f*time) - 16*A*pi^2*f^2*sin(4*pi*f*time).*cos(2*pi*f*time);


% % Line motion: Horizontal
% time=0:Tres:40/f;
% Middle=ceil(38*length(time)/40);
% Start=20;
% xef=A*cos(2*pi*f*time+phi);
% yef=2.5*ones(size(time));
% Dxef=-(2*pi*f)*A*sin(2*pi*f*time+phi);
% Dyef= 2*zeros(size(time));
% D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time+phi);
% D2yef= 2*zeros(size(time));

pos=[xef;yef];
Dpos=[Dxef;Dyef];
D2pos=[D2xef;D2yef];
hTPath=figure('name','TPath');
    plot(xef(Start:end),yef(Start:end),'linewidth',2.5,'linestyle','-','color','g')
    hold off
    axis equal
    xlabel('X (m)')
    ylabel('Y (m)')



%
% hTPath=figure('name','TPath');
hJointTrajectory=figure('name','Joint Trajectories');
hPowerTime=figure('name','Desired Power vs Time');
hTorqueTime=figure('name','Desired Torque vs Time');
hTorqueAngle=figure('name','Desired Torque vs Angle');
hPoincare=figure('name','Poincare');

%%
for ii=1:1
    
    ii

    q1=rand*2*pi;
    q2=rand*2*pi;
    q3=rand*2*pi;
    q1=deg2rad(-60);
    q2=deg2rad(48.031);
    q3=deg2rad(-48.031);

   
   

    for tt=1:length(time)-1
        lPos=L*[cos(q1(tt))+cos(q1(tt)+q2(tt))+cos(q1(tt)+q2(tt)+q3(tt));
                sin(q1(tt))+sin(q1(tt)+q2(tt))+sin(q1(tt)+q2(tt)+q3(tt))];

        cPos=pos(:,tt);
        dx=cPos-lPos;

        JJ=L*[-sin(q1(tt))-sin(q1(tt)+q2(tt))-sin(q1(tt)+q2(tt)+q3(tt)), -sin(q1(tt)+q2(tt))-sin(q1(tt)+q2(tt)+q3(tt)), -sin(q1(tt)+q2(tt)+q3(tt));
               cos(q1(tt))+cos(q1(tt)+q2(tt))+cos(q1(tt)+q2(tt)+q3(tt)),  cos(q1(tt)+q2(tt))+cos(q1(tt)+q2(tt)+q3(tt)),  cos(q1(tt)+q2(tt)+q3(tt))]; 

        dq=JJ'*(JJ*JJ')^-1*dx;

        q1(tt+1)=q1(tt)+dq(1);
        q2(tt+1)=q2(tt)+dq(2);
        q3(tt+1)=q3(tt)+dq(3);

    end
    
    q1=InRangeShifter(q1);
    q2=InRangeShifter(q2);
    q3=InRangeShifter(q3);
    RPosVal=L*[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
            sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];


    % Required torque for desired path
    q=[q1;q2;q3];
    Dq=differential(q,time,Tres);
    D2q=differential(Dq,time,Tres);
    torqueDesire=TorqueCalculator3(D2q,Dq,q,g,m,m,m,L,L,L);



    % show    
%     figure(hTPath)
% %         plot(xef(Start:end),yef(Start:end),'linewidth',2.5,'linestyle','-','color','g')
% %         hold on 
%         plot(RPosVal(1,Start:end),RPosVal(2,Start:end),'linewidth',2,'linestyle','-.','color','b')
%         legend('Desired','Static Path Planing')
%         hold off
%         axis equal

    figure(hJointTrajectory)
        subplot(3,1,1)
        plot(time(Start:end),q1(Start:end))
        grid on
        hold all
        ylabel('q_1 (deg)')
        subplot(3,1,2)
        plot(time(Start:end),q2(Start:end))
        grid on
        ylabel('q_2 (deg)')
        hold all
        subplot(3,1,3)
        plot(time(Start:end),q3(Start:end))
        grid on
        ylabel('q_3 (deg)')
        xlabel('time (s)')
        hold all
        


    figure(hPowerTime)
        subplot(3,1,1)
        plot(time(Start:end),torqueDesire(1,Start:end).*Dq(1,Start:end))
        grid on
        ylabel('power_1 (w)')
        hold all
        subplot(3,1,2)
        plot(time(Start:end),torqueDesire(2,Start:end).*Dq(2,Start:end))
        ylabel('power_2 (w)')
        hold all
        grid on
        subplot(3,1,3)
        plot(time(Start:end),torqueDesire(3,Start:end).*Dq(3,Start:end))
        hold all
        ylabel('power_3 (w)')
        xlabel('time (s)')
        grid on


    figure(hTorqueTime)
        subplot(3,1,1)
        plot(time(Start:end),torqueDesire(1,Start:end))
        grid on
        ylabel('\tau_1 (N.m)')
        hold all
        subplot(3,1,2)
        plot(time(Start:end),torqueDesire(2,Start:end))
        grid on
        ylabel('\tau_2 (N.m)')
        hold all
        subplot(3,1,3)
        plot(time(Start:end),torqueDesire(3,Start:end))
        grid on
        ylabel('\tau_3 (N.m)')
        hold all
        xlabel('time (s)')
        

    figure(hTorqueAngle)
        subplot(3,1,1)
        hp=plot(rad2deg(q(1,Middle:end)),torqueDesire(1,Middle:end),'linewidth',2);
        hold all
        plot(rad2deg(q(1,Middle)),torqueDesire(1,Middle),'linestyle','none','marker','*',...
            'markersize',10,'markerfacecolor',get(hp,'color'),'markeredgecolor',get(hp,'color'))
        title('Initial Desired Torque-Angle Profile','FontWeight','normal','FontSize',16,'FontName','Times');
        xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('\tau_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        set(gca,'YMinorGrid','on')
        grid on

        subplot(3,1,2)
        hp=plot(rad2deg(q(2,Middle:end)),torqueDesire(2,Middle:end),'linewidth',2);
        hold all
        plot(rad2deg(q(2,Middle)),torqueDesire(2,Middle),'linestyle','none','marker','*',...
            'markersize',10,'markerfacecolor',get(hp,'color'),'markeredgecolor',get(hp,'color'))
        xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('\tau_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        set(gca,'YMinorGrid','on')
        grid on

        subplot(3,1,3)
        hp=plot(rad2deg(q(3,Middle:end)),torqueDesire(3,Middle:end),'linewidth',2);
        hold all
        plot(rad2deg(q(3,Middle)),torqueDesire(3,Middle),'linestyle','none','marker','*',...
            'markersize',10,'markerfacecolor',get(hp,'color'),'markeredgecolor',get(hp,'color'))
        xlabel('q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('\tau_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        set(gca,'YMinorGrid','on')
        grid on


    figure(hPoincare)
        subplot(3,1,1)
        hp=plot(rad2deg(q(1,Middle:10:end)),rad2deg(Dq(1,Middle:10:end)),'linewidth',2);
        hold all
        plot(rad2deg(q(1,Middle)),rad2deg(Dq(1,Middle)),'linestyle','none','marker','*',...
            'markersize',10,'markerfacecolor',get(hp,'color'),'markeredgecolor','k');%get(hp,'color')),get(hp,'color'))
        xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('Dq_1 (deg/s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        set(gca,'YMinorGrid','on')
        grid on

        subplot(3,1,2)
        hp=plot(rad2deg(q(2,Middle:10:end)),rad2deg(Dq(2,Middle:10:end)),'linewidth',2);
        hold all
        plot(rad2deg(q(2,Middle)),rad2deg(Dq(2,Middle)),'linestyle','none','marker','*',...
            'markersize',10,'markerfacecolor',get(hp,'color'),'markeredgecolor','k');%get(hp,'color')),get(hp,'color'))
        xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('Dq_2 (deg/s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        set(gca,'YMinorGrid','on')
        grid on

        subplot(3,1,3)
        hp=plot(rad2deg(q(3,Middle:10:end)),rad2deg(Dq(3,Middle:10:end)),'linewidth',2);
        hold all
        plot(rad2deg(q(3,Middle)),rad2deg(Dq(3,Middle)),'linestyle','none','marker','*', ...
            'markersize',10,'markerfacecolor',get(hp,'color'),'markeredgecolor','k');%get(hp,'color'))
        xlabel('q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('Dq_3 (deg/s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        set(gca,'YMinorGrid','on')
        hold all
        grid on
    
end
% Y=[q1(1:end)',q2(1:end)',q3(1:end)'];
%%
figure(hTPath);
SavePdf([Name,'_Task'])
figur(hJointTrajectory);
SavePdf([Name,'_TrajectoryTime'])
figure(hPowerTime);
SavePdf([Name,'_PowerTime'])
figure(hTorqueTime);
SavePdf([Name,'_TorqueTime'])
figure(hTorqueAngle);
SavePdf([Name,'_TorqueAngle'])
figure(hPoincare);
SavePdf([Name,'_Poincare'])


% AnimBot3DOF(time(1:end),Y,L);

