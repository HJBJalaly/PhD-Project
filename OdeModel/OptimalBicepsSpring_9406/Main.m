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

% 
% % Circle motion
% f=0.5;
% phi=pi/2;
% Tres=0.01;
% time=0:Tres:40/f;
% Middle=ceil(39*length(time)/40);
% xef=A*cos(2*pi*f*time)+0;
% yef=A*sin(2*pi*f*time)+.5;
% Dxef=-(2*pi*f)*A*sin(2*pi*f*time);
% Dyef= (2*pi*f)*A*cos(2*pi*f*time);
% D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time);
% D2yef=-(2*pi*f)^2*A*sin(2*pi*f*time);
% q1=deg2rad(-60);
% q2=deg2rad(48.031);

% Ellipose motion
Name='Ellipose';
f=0.5;
phi=pi/2;
Tres=0.01;
time=0:Tres:40/f;
Middle=ceil(39*length(time)/40);
% Start=20;
xef=A*cos(2*pi*f*time)+0;
yef=A/2*sin(2*pi*f*time)+.5;
Dxef=-(2*pi*f)*A*sin(2*pi*f*time);
Dyef= (2*pi*f)*A/2*cos(2*pi*f*time);
D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time);
D2yef=-(2*pi*f)^2*A/2*sin(2*pi*f*time);
q1=deg2rad(-60);
q2=deg2rad(48.031);
 

% % Line motion: Horizontal
% f=1;
% phi=0;
% Tres=0.002;
% time=0:Tres:20/f;
% Middle=ceil(19*length(time)/20);
% xef=A*cos(2*pi*f*time+phi);
% yef=2.5*ones(size(time));
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
    
Y=[q1(1:end)',q2(1:end)'];
% AnimBot3DOF(time(1:end),Y,L);

IntAbsUdqDesire_Opt=sum(sum(abs(torqueDesire(:,Middle:end).*Dq(:,Middle:end)),2))*Tres

%%  Generate Initial value for Optimization
% close all

% DoF system
nn=2; % number of joints
% DoF of Optimization 
% rQ=20; % Degree of joint trajectory
rU=2; % Degree of passive torque
rB=5; % Degree of passive torque
SampleRate=1;
% % Landa for [D2Q  D2BQ ]
% SeletcStr={'DQ','D2Q'};
% SelectLanda=[0 1];
Landa=[0 .000001 .1 0.1];


Time=time(Middle:end)-time(Middle);
Q1=InRangeShifter(q1(Middle:end));
Q2=InRangeShifter(q2(Middle:end));
QJ=[Q1;Q2];
Qhat1=Q1+Q2;
% /DQ2=Dq(2,Middle:end);
QhatJ=[Qhat1];


XEF=xef(Middle:end);
YEF=yef(Middle:end);
DXEF=Dxef(Middle:end);
DYEF=Dyef(Middle:end);
D2XEF=D2xef(Middle:end);
D2YEF=D2yef(Middle:end);
TorqueDesQ1=torqueDesire(1,Middle-1:end-1);
TorqueDesQ2=torqueDesire(2,Middle-1:end-1);
TorqueDesire=[TorqueDesQ1;TorqueDesQ2];


% RPosVal=L*[cos(Q1val)+cos(Q1val+Q2val);
%         sin(Q1val)+sin(Q1val+Q2val)];

    
[BetaOptimal,ThetaOptimal,Cost,TorquePassiveOptimal,TorqueBicepsOptimal]= ...
            LSParamPolyBiceps(QJ,QhatJ,TorqueDesire,nn,rU,rB,Landa,SampleRate);

% BetaOp1=OptimalSpring((1-1)*(rU+1)+1:(1)*(rU+1));
% BetaOp2=OptimalSpring((2-1)*(rU+1)+1:(2)*(rU+1));
% ThetaOp1=OptimalSpring(nn*(rU+1)+(1-1)*(rB+1)+1:nn*(rU+1)+(1)*(rB+1));
% BetaOptimal=[BetaOp1;BetaOp2];
% ThetaOptimal=[ThetaOp1];


TorquePassiveOptimalQ1=TorquePassiveOptimal(1,:);
%polyval(BetaOptimal(0*(rU+1)+1:1*(rU+1)),Q1);
TorquePassiveOptimalQ2=TorquePassiveOptimal(2,:);
%polyval(BetaOptimal(1*(rU+1)+1:2*(rU+1)),Q2);
% TorquePassiveOptimal=[TorquePassiveOptimalQ1; TorquePassiveOptimalQ2];

TorqueBicepsOptimalQ1=TorqueBicepsOptimal(1,:);
%polyval(ThetaOptimal(0*(rB+1)+1:1*(rB+1)),Qhat1);
% TorqueBicepsOptimal=[TorqueBicepsOptimalQ1];

TorqueActive=TorqueDesire-TorquePassiveOptimal-[TorqueBicepsOptimal;zeros(size(TorquePassiveOptimalQ1))]-[zeros(size(TorquePassiveOptimalQ1));TorqueBicepsOptimal];
WorkActuator    =sum(sum(abs(TorqueActive.*Dq(:,Middle:end)),2))*Tres;
WorkDesire      =sum(sum(abs(TorqueDesire.*Dq(:,Middle:end)),2))*Tres;
Title=sprintf('%10s  %15s %15s'  ,   'Cost','Req. Work','Act. Work');
Result=sprintf('% 11.2f  % 11.2f % 15.2f'  ,   Cost,WorkDesire,WorkActuator);
disp('-------------------------------------------------------')
disp(Title)
disp(Result)

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
    legend('Required Torque','Mono Spring Torque')
    subplot(2,2,2)
    plot(rad2deg(Qhat1),TorqueDesQ1,'linewidth',2)
    hold on
    plot(rad2deg(Qhat1),TorqueBicepsOptimalQ1,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    hold off
    grid on
    title('Parallel Biceps Compliance','FontWeight','bold','FontSize',15,'FontName','mwa_cmb10');
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('q_1+q_2 (deg/s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Biceps Spring Torque')
    
    
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
%     subplot(2,2,4)
%     plot(rad2deg(DQ2),TorqueDesQ2,'linewidth',2)
%     hold on
%     plot(rad2deg(DQ2),TorqueDamperOptimalQ2,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
%     hold off
%     grid on
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     xlabel('Dq_2 (deg/s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
%     ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
%     legend('Required Torque','Damper Torque')

figure('name','Time Torques')
    subplot(2,1,1)
    plot(Time,TorqueDesQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorquePassiveOptimalQ1,'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueBicepsOptimalQ1,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    plot(Time,TorqueActive(1,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    L=legend('Required Torque','Mono Spring Torque','Biceps Spring Torque','Actuator Torque');
    set(L,'orientation','horizontal')
    subplot(2,1,2)
    plot(Time,TorqueDesQ2,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorquePassiveOptimalQ2,'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueBicepsOptimalQ1,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    plot(Time,TorqueActive(2,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    L=legend('Required Torque','Mono Spring Torque','Biceps Spring Torque','Actuator Torque');
    set(L,'orientation','horizontal')
    
figure('name','Substitude Toruq')
    subplot(2,2,1)
    plot(rad2deg(Q1),TorqueDesQ1-TorqueBicepsOptimalQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(rad2deg(Q1),TorquePassiveOptimalQ1,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    title('Mono Parallel Compliance','FontWeight','bold','FontSize',15,'FontName','mwa_cmb10');
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_r_1-u_p_b_1^* (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Spring Torque')
    subplot(2,2,2)
    plot(rad2deg(Qhat1),TorqueDesQ1-TorquePassiveOptimalQ1,'linewidth',2)
    hold on
    plot(rad2deg(Qhat1),TorqueBicepsOptimalQ1,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    hold off
    grid on
    title('Biceps Parallel Compliance','FontWeight','bold','FontSize',15,'FontName','mwa_cmb10');
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('q_1+q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_1-up_m_o_n_o_1^* (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Biceps compliance Torque')
    
    
    subplot(2,2,3)
    plot(rad2deg(Q2),TorqueDesQ2-TorqueBicepsOptimalQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(rad2deg(Q2),TorquePassiveOptimalQ2,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_r_2-u_p_b_1^*','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Spring Torque')
    subplot(2,2,4)
    plot(rad2deg(Qhat1),TorqueDesQ2+TorqueDesQ1-TorquePassiveOptimalQ2-TorquePassiveOptimalQ1,'linewidth',2)
    hold on
    plot(rad2deg(Qhat1),2*TorqueBicepsOptimalQ1,'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    hold off
    grid on
    title('Biceps Parallel Compliance','FontWeight','bold','FontSize',15,'FontName','mwa_cmb10');
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    xlabel('q_1+q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_r_1+u_r_2 - u_p_u_1^* - u_p_u_2^* (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    legend('Required Torque','Biceps compliance Torque')
 