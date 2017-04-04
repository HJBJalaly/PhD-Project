function Main()



%% Initialization

clear
close all
home
rand('twister', sum(100*clock));

%% create a active sample motion with 2 active joints as initial input of optimization
home
% clear
clc
% envirment
g=9.81*1;
% parameters of robot
m=.5;
L=.5;
% EF motion
A=.3;

% Circle motion
f=0.5;
phi=pi/2;
Tres=0.01;
time=0:Tres:40/f;
Middle =ceil(39*length(time)/40);
Middle2=ceil(38*length(time)/40);
xef=A*cos(2*pi*f*time)+0;
yef=A*sin(2*pi*f*time)+.65;
Dxef=-(2*pi*f)*A*sin(2*pi*f*time);
Dyef= (2*pi*f)*A*cos(2*pi*f*time);
D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time);
D2yef=-(2*pi*f)^2*A*sin(2*pi*f*time);
q1=deg2rad(-60);
q2=deg2rad(48.031);

% % Line motion: Horizontal
% f=1;
% phi=0;
% Tres=0.002;
% time=0:Tres:20/f;
% Middle=ceil(19*length(time)/20);
% Middle2=ceil(38*length(time)/40);
% xef=A*cos(2*pi*f*time+phi);
% yef=.25*ones(size(time));
% Dxef=-(2*pi*f)*A*sin(2*pi*f*time+phi);
% Dyef= 2*zeros(size(time));
% D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time+phi);
% D2yef= 2*zeros(size(time));
% q1=deg2rad( 22.7023); % 41.7023 deg for y=2.5m  ,  22.2 deg for y=2.0m
% q2=deg2rad(29.2663); % 18.2663 deg for y=2.5m  ,  29.468 deg for y=2.0m
% q3=deg2rad(80.3377); % 44.3377 deg for y=2.5m  ,  71.431 deg for y=2.0m
% % 

% % Line motion: Vertical
% phi=0;
% yef=A*cos(2*pi*f*time+phi);
% xef=2.5*ones(size(time));
% Dyef=-(2*pi*f)*A*sin(2*pi*f*time+phi);
% Dxef= 2*zeros(size(time));
% D2yef=-(2*pi*f)^2*A*cos(2*pi*f*time+phi);
% D2xef= 2*zeros(size(time));
% q1=deg2rad( 41.7023); % 41.7023 deg for y=2.5m  ,  22.2 deg for y=2.0m
% q2=deg2rad(18.2663); % 18.2663 deg for y=2.5m  ,  29.468 deg for y=2.0m
% q3=deg2rad(44.3377); % 44.3377 deg for y=2.5m  ,  71.431 deg for y=2.0m


% Line motion
% time=0:Tres:4/f;
% phi=pi;
% xef=A*(1-cos(2*pi*f/2*time+phi))-1;
% yef=-A/4*sin(2*pi*f*time+phi)+2.5;
% 
% Dyef=(2*pi*f/2)*A*sin(2*pi*f/2*time+phi);
% Dxef= -(2*pi*f)*A/4*cos(2*pi*f*time+phi);
% 
% D2yef=(2*pi*f/2)^2*A*cos(2*pi*f/2*time+phi);
% D2xef=(2*pi*f)^2*A/4*sin(2*pi*f*time+phi);
% q1=deg2rad( 41.7023); % 41.7023 deg for y=2.5m  ,  22.2 deg for y=2.0m
% q2=deg2rad(18.2663); % 18.2663 deg for y=2.5m  ,  29.468 deg for y=2.0m
% q3=deg2rad(44.3377); % 44.3377 deg for y=2.5m  ,  71.431 deg for y=2.0m
% 


pos=[xef;yef];
Dpos=[Dxef;Dyef];
D2pos=[D2xef;D2yef];


% Joint Trajectories

for tt=1:length(time)-1
    lPos=L*[cos(q1(tt))+cos(q1(tt)+q2(tt));
            sin(q1(tt))+sin(q1(tt)+q2(tt))];

    cPos=pos(:,tt);
    dx=cPos-lPos;

    JJ=L*[-sin(q1(tt))-sin(q1(tt)+q2(tt)), -sin(q1(tt)+q2(tt));
           cos(q1(tt))+cos(q1(tt)+q2(tt)),  cos(q1(tt)+q2(tt))]; 

    dq=JJ'*(JJ*JJ')^-1*dx;

    q1(tt+1)=q1(tt)+dq(1);
    q2(tt+1)=q2(tt)+dq(2);
    
end

rPosVal=L*[cos(q1)+cos(q1+q2);
        sin(q1)+sin(q1+q2)];

    
% Required torque for desired path
q=[q1;q2];
Dq=differential(q,time,Tres);
D2q=differential(Dq,time,Tres);
torqueDesire=TorqueCalculator2(D2q,Dq,q,g,m,m,L,L);

% close all
figure('name','Path (3DoF)')
    plot(xef(Middle:end),yef(Middle:end),'linewidth',2.5,'linestyle','-','color','g')
    hold on 
    plot(rPosVal(1,Middle:end),rPosVal(2,Middle:end),'linewidth',2,'linestyle','-.','color','b')
    legend('Desired','Static Path Planing')
    hold off
    axis equal
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('x (m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('y (m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    grid on

figure('name','Joint Trajectories')
    subplot(2,1,1)
    plot(time(Middle:end),q(1,Middle:end))
    ylabel('q_1 (rad)')
    grid on
    subplot(2,1,2)
    plot(time(Middle:end),q(2,Middle:end))
    grid on
    ylabel('q_2 (rad)')
    xlabel('time (s)')
    

figure('name','Desired Power vs Time')
    subplot(2,1,1)
    plot(time(Middle:end),torqueDesire(1,Middle:end).*Dq(1,Middle:end))
    ylabel('Power_1')   
    grid on
    subplot(2,1,2)
    plot(time(Middle:end),torqueDesire(2,Middle:end).*Dq(2,Middle:end))
    grid on
    ylabel('Power_2')   
    xlabel('time (s)')   
    
    
    
figure('name','Desired Torque vs Time')
    plot(time(Middle:end),torqueDesire(:,Middle:end))
    legend('u_1','u_2')
    xlabel('time (s)')
    ylabel('u')
    grid on
    


figure('name','Desired Torque vs Angle')
    subplot(2,2,1)
    plot(rad2deg(q(1,Middle:end)),torqueDesire(1,Middle:end),'linewidth',2,'linestyle','none','marker','*')
    title('Desired Torque-Angle Profile','FontSize',16);
    set(gca,'YMinorGrid','on')
    xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    grid on
    
    subplot(2,2,3)
    plot((rad2deg(q(2,Middle:end))),torqueDesire(2,Middle:end),'linewidth',2)
    xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
    subplot(2,2,2)
    plot(rad2deg(Dq(1,Middle:end)),torqueDesire(1,Middle:end),'linewidth',2)
    title('Desired Torque-Velocity Profile','FontSize',16);
    xlabel('Dq_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
    subplot(2,2,4)
    plot(rad2deg(Dq(2,Middle:end)),torqueDesire(2,Middle:end),'linewidth',2)
    xlabel('Dq_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
Y=[q1(Middle:end)',q2(Middle:end)'];
AnimBot3DOF(time(Middle:end),Y,L);

IntAbsUdqDesire_Opt=sum(sum(abs(torqueDesire(:,Middle:end).*Dq(:,Middle:end)),2))*Tres

%%  Generate Initial value for Optimization
% close all

% DoF system
nn=2; % number of joints
% DoF of Optimization 
% rQ=20; % Degree of joint trajectory
rP=10; % Degree of Compliance polunomial torque
rD=0; % Degree of Damper polynomial torque
% WeightMatrix
% Weight=[1 1]';

% % Landa for [DQ  D2q ]
% SeletcStr={'DQ','D2Q'};
% SelectLanda=[0 1];
% Landa=[1e-3 1e-3];


Time=time(Middle-15:end-15)-time(Middle-15);
Q1=q1(Middle-15:end-15);
Q2=q2(Middle-15:end-15);
QJ=[Q1;Q2];
DQ1=Dq(1,Middle-15:end-15);
DQ2=Dq(2,Middle-15:end-15);
DQJ=[DQ1;DQ2];


XEF=xef(Middle-15:end-15);
YEF=yef(Middle-15:end-15);
DXEF=Dxef(Middle-15:end-15);
DYEF=Dyef(Middle-15:end-15);
D2XEF=D2xef(Middle-15:end-15);
D2YEF=D2yef(Middle-15:end-15);
TorqueDesQ1=torqueDesire(1,Middle-15:end-15);
TorqueDesQ2=torqueDesire(2,Middle-15:end-15);
TorqueDesire=[TorqueDesQ1;TorqueDesQ2];


% [Alpha_Q1,BezireCoef_q1]= BezierCoeffinet(Time,Q1,rQ);
% [Alpha_Q2,BezireCoef_q2]= BezierCoeffinet(Time,Q2,rQ);
% 
% Q1val=polyval(Alpha_Q1,Time);
% Q2val=polyval(Alpha_Q2,Time);
% QVal=[Q1val;Q2val];

% RPosVal=L*[cos(Q1val)+cos(Q1val+Q2val);
%         sin(Q1val)+sin(Q1val+Q2val)];

    
[BetaOp1,ThetaOp1] = LSParamPoly(Q1',DQ1,TorqueDesQ1',rP,rD,0);
[BetaOp2,ThetaOp2] = LSParamPoly(Q2',DQ2,TorqueDesQ2',rP,rD,0);
BetaOptimal=[BetaOp1;BetaOp2];
ThetaOptimal=[ThetaOp1;ThetaOp2];


% TorquePassiveQ1val=polyval(CoefBLS_UPassive1,Q1);
% TorquePassiveQ2val=polyval(CoefBLS_UPassive2,Q2);
% TorquePassiveVal=[TorquePassiveQ1val; TorquePassiveQ2val];
    

% Integral Matrix
Cost=0;
for i=1:nn
    QQ=[];
    DQ=[];
    for j=1:length(QJ(i,:))
         QQ(j,:) = QJ(i,j).^(rP:-1:0)';
         DQ(j,:) = DQJ(i,j).^(rD:-1:0)';
    end
    
    Cost=Cost + ...
          ( 1/2*(TorqueDesire(i,:)' - QQ*BetaOptimal((i-1)*(rP+1)+1:(i)*(rP+1))-DQ*ThetaOptimal((i-1)*(rD+1)+1:(i)*(rD+1)) )'* (TorqueDesire(i,:)' - QQ*BetaOptimal((i-1)*(rP+1)+1:(i)*(rP+1))-DQ*ThetaOptimal((i-1)*(rD+1)+1:(i)*(rD+1)) ));                      
end


TorquePassiveOptimalQ1=polyval(BetaOptimal(0*(rP+1)+1:1*(rP+1)),Q1);
TorquePassiveOptimalQ2=polyval(BetaOptimal(1*(rP+1)+1:2*(rP+1)),Q2);
TorquePassiveOptimal=[TorquePassiveOptimalQ1; TorquePassiveOptimalQ2];

TorqueDamperOptimalQ1=polyval(ThetaOptimal(0*(rD+1)+1:1*(rD+1)),DQ1);
TorqueDamperOptimalQ2=polyval(ThetaOptimal(1*(rD+1)+1:2*(rD+1)),DQ2);
TorqueDamperOptimal=[TorqueDamperOptimalQ1; TorqueDamperOptimalQ2];

TorqueActive=TorqueDesire-TorquePassiveOptimal-TorqueDamperOptimal;
WorkActuator    =sum(sum(abs(TorqueActive.*Dq(:,Middle:end)),2))*Tres;
WorkCompliance  =sum(sum(abs(TorquePassiveOptimal.*Dq(:,Middle:end)),2))*Tres;
WorkDamper      =sum(sum(abs(TorqueDamperOptimal.*Dq(:,Middle:end)),2))*Tres;
WorkDesire      =sum(sum(abs(TorqueDesire.*Dq(:,Middle:end)),2))*Tres;
Title=sprintf('%10s  %15s %15s %15s %15s'  ,   'Cost','Req. Work','Cmpl. Work','Dmpl. Work','Act. Work');
Result=sprintf('% 11.2f  % 11.2f % 15.2f % 15.2f % 15.2f'  ,   Cost,WorkDesire,WorkCompliance,WorkDamper,WorkActuator);
disp('-------------------------------------------------------')
disp(Title)
disp(Result)


% show time
% close all

% figure('name','compare trajectory')
%     subplot(2,1,1)
%     plot(Time,Q1)
%     hold all
%     plot(Time,Q1val,'-.')
%     hold off
%     grid on
%     subplot(2,1,2)
%     plot(Time,Q2)
%     hold all
%     plot(Time,Q2val,'-.')
%     hold off
%     grid on
    
% figure('name','Path (3DoF)')
%     plot(xef,yef,'linewidth',2.5,'linestyle','-','color','g')
%     hold on 
%     plot(RPosVal(1,:),RPosVal(2,:),'linewidth',2,'linestyle','-.','color','b')
%     legend('Desired','Static Path Planing')
%     hold off
%     axis equal

figure('name','Passive Torques')
    subplot(2,2,1)
    plot(rad2deg(Q1),TorqueDesQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(rad2deg(Q1),TorquePassiveOptimalQ1,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    title('Parallel Compliance','FontWeight','bold','FontSize',15,'FontName','mwa_cmb10');
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Spring Torque')
    subplot(2,2,2)
    plot(rad2deg(DQ1),TorqueDesQ1-0*TorquePassiveOptimalQ1,'linewidth',2)
    hold on
    plot(rad2deg(DQ1),TorqueDamperOptimalQ1,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    hold off
    grid on
    title('Parallel Damper','FontWeight','bold','FontSize',15,'FontName','mwa_cmb10');
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('Dq_1 (deg/s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Damper Torque')
    
    
    subplot(2,2,3)
    plot(rad2deg(Q2),TorqueDesQ2,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(rad2deg(Q2),TorquePassiveOptimalQ2,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Spring Torque')
    subplot(2,2,4)
    plot(rad2deg(DQ2),TorqueDesQ2,'linewidth',2)
    hold on
    plot(rad2deg(DQ2),TorqueDamperOptimalQ2,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('Dq_2 (deg/s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Damper Torque')

figure('name','Passive Torques')
    subplot(2,1,1)
    plot(Time,TorqueDesQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorquePassiveOptimalQ1,'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueDamperOptimalQ1,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    plot(Time,TorqueActive(1,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    L=legend('Required Torque','Spring Torque','Damper Torque','Actuator Torque');
    set(L,'orientation','horizontal')
    subplot(2,1,2)
    plot(Time,TorqueDesQ2,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorquePassiveOptimalQ2,'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueDamperOptimalQ2,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    plot(Time,TorqueActive(2,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('Time (s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    L=legend('Required Torque','Spring Torque','Damper Torque','Actuator Torque');
    set(L,'orientation','horizontal')
    
    
figure('name','Substitude Toruq')
    subplot(2,2,1)
    plot(rad2deg(Q1),TorqueDesQ1-TorqueDamperOptimalQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(rad2deg(Q1),TorquePassiveOptimalQ1,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    title('Parallel Compliance','FontWeight','bold','FontSize',15,'FontName','mwa_cmb10');
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_1-ud_1^* (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Spring Torque')
    subplot(2,2,2)
    plot(rad2deg(DQ1),TorqueDesQ1-TorquePassiveOptimalQ1,'linewidth',2)
    hold on
    plot(rad2deg(DQ1),TorqueDamperOptimalQ1,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    hold off
    grid on
    title('Parallel Damper','FontWeight','bold','FontSize',15,'FontName','mwa_cmb10');
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('Dq_1 (deg/s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_1-up_1^* (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Damper Torque')
    
    
    subplot(2,2,3)
    plot(rad2deg(Q2),TorqueDesQ2-TorqueDamperOptimalQ2,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(rad2deg(Q2),TorquePassiveOptimalQ2,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_2-ud_2^* (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Spring Torque')
    
    legend('Required Torque','PassiveTorque')
    subplot(2,2,4)
    plot(rad2deg(DQ2),TorqueDesQ2-TorquePassiveOptimalQ2,'linewidth',2)
    hold on
    plot(rad2deg(DQ2),TorqueDamperOptimalQ2,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('Dq_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_2-up_2^* (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Damper Torque')
    
%%

figure('name','Passive Torques')
    subplot(2,2,1)
        plot(rad2deg(Q1),TorqueDesQ1,'linewidth',2,'linestyle','-','color','b')
        hold on
        plot(rad2deg(Q1),TorquePassiveOptimalQ1,'linewidth',2,'linestyle','-.','color','g')
        hold off
        grid on
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        legend('u_r_1','u_p_1')
    subplot(2,2,2)
        plot(rad2deg(Q2),TorqueDesQ2,'linewidth',2,'linestyle','-','color','b')
        hold on
        plot(rad2deg(Q2),TorquePassiveOptimalQ2,'linewidth',2,'linestyle','-.','color','g')
        hold off
        grid on
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        legend('u_r_2','u_p_2')
    subplot(2,2,3)
        plot(rad2deg(Q1),TorqueActive(1,:),'linewidth',2,'linestyle','-','color','r')
        hold off
        grid on
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('u_a_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        legend('u_a_1')
    subplot(2,2,4)
        plot(rad2deg(Q2),TorqueActive(2,:),'linewidth',2,'linestyle','-','color','r')
        hold off
        grid on
        set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('u_a_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        legend('u_a_2')
        
figure('name','Passive Torques')
    subplot(2,1,1)
    plot(Time,TorqueDesQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorquePassiveOptimalQ1,'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueDamperOptimalQ1,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    plot(Time,TorqueActive(1,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    L=legend('Required Torque','Spring Torque','Damper Torque','Actuator Torque');
    set(L,'orientation','horizontal')
    subplot(2,1,2)
    plot(Time,TorqueDesQ2,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorquePassiveOptimalQ2,'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueDamperOptimalQ2,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    plot(Time,TorqueActive(2,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('Time (s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    L=legend('Required Torque','Spring Torque','Damper Torque','Actuator Torque');
    set(L,'orientation','horizontal')
    
    
figure('name','Substitude Toruq')
    subplot(2,2,1)
    plot(rad2deg(Q1),TorqueDesQ1-TorqueDamperOptimalQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(rad2deg(Q1),TorquePassiveOptimalQ1,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    title('Parallel Compliance','FontWeight','bold','FontSize',15,'FontName','mwa_cmb10');
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_1-ud_1^* (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Spring Torque')
    subplot(2,2,2)
    plot(rad2deg(DQ1),TorqueDesQ1-TorquePassiveOptimalQ1,'linewidth',2)
    hold on
    plot(rad2deg(DQ1),TorqueDamperOptimalQ1,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    hold off
    grid on
    title('Parallel Damper','FontWeight','bold','FontSize',15,'FontName','mwa_cmb10');
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('Dq_1 (deg/s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_1-up_1^* (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Damper Torque')
    
    
    subplot(2,2,3)
    plot(rad2deg(Q2),TorqueDesQ2-TorqueDamperOptimalQ2,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(rad2deg(Q2),TorquePassiveOptimalQ2,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_2-ud_2^* (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Spring Torque')
    
    legend('Required Torque','PassiveTorque')
    subplot(2,2,4)
    plot(rad2deg(DQ2),TorqueDesQ2-TorquePassiveOptimalQ2,'linewidth',2)
    hold on
    plot(rad2deg(DQ2),TorqueDamperOptimalQ2,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('Dq_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_2-up_2^* (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Damper Torque')
    
    
 