function CF3b_Main()



%% Initialization

clear
close all
home
rand('twister', sum(100*clock));

%% create a active sample motion with 3 active joints as initial input of optimization
home
% clear
clc
% envirment
g=9.81*1;
% parameters of robot
m=1;
L=1;
% EF motion
A=.75;
% 
% % Circle motion
% f=0.5;
% phi=pi/2;
% Tres=0.005;
% time=0:Tres:11/f;
% Middle1=ceil(9*length(time)/11);
% Middle2=ceil(10*length(time)/11);
% xef=A*cos(2*pi*f*time)+0;
% yef=A*sin(2*pi*f*time)+2;
% Dxef=-(2*pi*f)*A*sin(2*pi*f*time);
% Dyef= (2*pi*f)*A*cos(2*pi*f*time);
% D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time);
% D2yef=-(2*pi*f)^2*A*sin(2*pi*f*time);
% q1=deg2rad(-60);
% q2=deg2rad(48.031);
% q3=deg2rad(-51.317);
  
% % Line motion: Horizontal
% f=1;
% phi=0;
% Tres=0.002;
% time=0:Tres:21/f;
% Middle1=ceil(19*length(time)/21);
% Middle2=ceil(20*length(time)/21);
% xef=A*cos(2*pi*f*time+phi);
% yef=2.5*ones(size(time));
% Dxef=-(2*pi*f)*A*sin(2*pi*f*time+phi);
% Dyef= 2*zeros(size(time));
% D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time+phi);
% D2yef= 2*zeros(size(time));
% q1=deg2rad( 22.7023); % 41.7023 deg for y=2.5m  ,  22.2 deg for y=2.0m
% q2=deg2rad(29.2663); % 18.2663 deg for y=2.5m  ,  29.468 deg for y=2.0m
% q3=deg2rad(80.3377); % 44.3377 deg for y=2.5m  ,  71.431 deg for y=2.0m

% Ellipose motion
f=0.5;
phi=pi/2;
Tres=0.01;
time=0:Tres:41/f;
Middle1=ceil(39*length(time)/41);
Middle2=ceil(40*length(time)/41);
Start=20;
xef=A*cos(2*pi*f*time)+0;
yef=A/2*sin(2*pi*f*time)+1.5;
Dxef=-(2*pi*f)*A*sin(2*pi*f*time);
Dyef= (2*pi*f)*A/2*cos(2*pi*f*time);
D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time);
D2yef=-(2*pi*f)^2*A/2*sin(2*pi*f*time);
q1=deg2rad(-60);
q2=deg2rad(48.031);
q3=deg2rad(-48.031);




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

rPosVal=L*[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
        sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];

    
% Required torque for desired path
q=[q1;q2;q3];
Dq=differential(q,time,Tres);
D2q=differential(Dq,time,Tres);
torqueDesire=TorqueCalculator3(D2q,Dq,q,g,m,m,m,L,L,L);
%%
% show    
% close all
figure('name','Path (3DoF)')
    plot(xef(Middle1:Middle2),yef(Middle1:Middle2),'linewidth',2.5,'linestyle','-','color','g')
    hold on 
    plot(rPosVal(1,Middle1:Middle2),rPosVal(2,Middle1:Middle2),'linewidth',2,'linestyle','-.','color','b')
    legend('Desired','Static Path Planing')
    hold off
    axis equal

figure('name','Joint Trajectories')
    subplot(3,1,1)
    plot(time(Middle1:Middle2),q1(Middle1:Middle2))
    grid on
    subplot(3,1,2)
    plot(time(Middle1:Middle2),q2(Middle1:Middle2))
    grid on
    subplot(3,1,3)
    plot(time(Middle1:Middle2),q3(Middle1:Middle2))
    grid on
    

figure('name','Desired Power vs Time')
%     subplot(3,1,1)
%     plot(time(Middle1:Middle2),torqueDesire(1,Middle1:Middle2).*Dq(1,Middle1:Middle2))
%     grid on
%     subplot(3,1,2)
%     plot(time(Middle1:Middle2),torqueDesire(2,Middle1:Middle2).*Dq(2,Middle1:Middle2))
%     grid on
%     subplot(3,1,3)
%     plot(time(Middle1:Middle2),torqueDesire(3,Middle1:Middle2).*Dq(3,Middle1:Middle2))
%     grid on
    plot(time(Middle1:Middle2),torqueDesire(:,Middle1:Middle2).*Dq(:,Middle1:Middle2),'linewidth',2)
    set(gca,'fontsize',12,'FontWeight','bold')
    legend('u_1','u_2','u_3')
    xlabel('time (s)','fontsize',14,'FontWeight','bold')
    ylabel('u','fontsize',14,'FontWeight','bold')
    grid on
    
    
figure('name','Desired Torque vs Time')
    plot(time(Middle1:Middle2),torqueDesire(:,Middle1:Middle2))
    legend('u_1','u_2','u_3')
    xlabel('time (s)')
    ylabel('u')
    grid on
    


figure('name','Desired Torque vs Angle')
    subplot(3,1,1)
    plot(rad2deg(q1(Middle1:Middle2)),torqueDesire(1,Middle1:Middle2),...
        'Color',[0.87058824300766 0.490196079015732 0],'linewidth',2)
    title('Initial Required Torque-Angle Profile','FontSize',16);
    xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
    subplot(3,1,2)
    plot(InRangeShifter(rad2deg(q2(Middle1:Middle2))),torqueDesire(2,Middle1:Middle2),...
            'Color',[0.87058824300766 0.490196079015732 0],'linewidth',2)
    xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
    subplot(3,1,3)
    plot(rad2deg(q3(Middle1:Middle2)),torqueDesire(3,Middle1:Middle2),...
        'Color',[0.87058824300766 0.490196079015732 0],'linewidth',2)
    xlabel('q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
    
Y=[q1(1:end)',q2(1:end)',q3(1:end)'];
% AnimBot3DOF(time(1:end),Y,L);

IntAbsUdqDesire_Opt=sum(sum(abs(torqueDesire(:,Middle1:Middle2).*Dq(:,Middle1:Middle2)),2))*Tres

%%  Generate Initial value for Optimization
% close all

tic
% DoF system
nn=3; % number of joints
% DoF of Optimization 
rQ=8; % Degree of joint trajectory
rU=3; % Degree of passive torque
rB=0;
% WeightMatrix
Weight=[3 2 1]';
Weight=Weight*1;
Sat=[1,1,1];
SampleRate=10;

% Landa for [DQ  D2q ]
SeletcStr={'DQ','D2Q'};
SelectLanda=[0 1];
Landa=[.05 1e-7 .05 1e-7];   % [landa_1* Beta'*Beta    Landa_2*(D2Q*Beta)'*(D2Q*Beta)
                           %  landa_3* Theta'*Theta  Landa_4*(D2Qhat*Theta)'*(D2Qhat*Theta) ]          
% Landa=[.000002 1e-7 .000002 1e-7];   % [landa_1* Beta'*Beta    Landa_2*(D2Q*Beta)'*(D2Q*Beta)
                           %  landa_3* Theta'*Theta  Landa_4*(D2Qhat*Theta)'*(D2Qhat*Theta) ]          


Time=time(Middle1:Middle2)-time(Middle1);
Q1=q1(Middle1:Middle2);
Q2=q2(Middle1:Middle2);
Q3=q3(Middle1:Middle2);
QJ=[Q1;Q2;Q3];

Qhat1=Q1+Q2;
Qhat2=Q2+Q3;
QhatJ=[Qhat1;Qhat2];

% Q_CF1_U2=load('temp2','Q_Opt');
% Q1=Q_Opt(1,:);
% Q2=Q_Opt(2,:);
% Q3=Q_Opt(3,:);
XEF=xef(Middle1:Middle2);
YEF=yef(Middle1:Middle2);
DXEF=Dxef(Middle1:Middle2);
DYEF=Dyef(Middle1:Middle2);
D2XEF=D2xef(Middle1:Middle2);
D2YEF=D2yef(Middle1:Middle2);
TorqueDesQ1=torqueDesire(1,Middle1:Middle2);
TorqueDesQ2=torqueDesire(2,Middle1:Middle2);
TorqueDesQ3=torqueDesire(3,Middle1:Middle2);
TorqueDesire=[TorqueDesQ1;TorqueDesQ2;TorqueDesQ3];


[Alpha_Q1,BezireCoef_q1]= BezierCoeffinet(Time,Q1,rQ);
[Alpha_Q2,BezireCoef_q2]= BezierCoeffinet(Time,Q2,rQ);
[Alpha_Q3,BezireCoef_q3]= BezierCoeffinet(Time,Q3,rQ);

Q1val=polyval(Alpha_Q1,Time);
Q2val=polyval(Alpha_Q2,Time);
Q3val=polyval(Alpha_Q3,Time);
QVal=[Q1val;Q2val;Q3val];


RPosVal=L*[cos(Q1val)+cos(Q1val+Q2val)+cos(Q1val+Q2val+Q3val);
        sin(Q1val)+sin(Q1val+Q2val)+sin(Q1val+Q2val+Q3val)];

    
[BetaOptimal,ThetaOptimal,CostActuation,CostD2Q,CostParaReg,TorquePassiveOptimal,TorqueBicepsOptimal]= ...
            OptimalParam(QJ,QhatJ,TorqueDesire,nn,rU,rB,Landa,Weight,SampleRate);

Cost=(CostActuation+ Landa(1:2:4)*CostParaReg+Landa(2:2:4)*CostD2Q)*Tres;
% Y=[Q1(1:end)',Q2(1:end)',Q3(1:end)'];
% AnimBot3DOF(Time,Y,L);


TorqueActive=TorqueDesire-TorquePassiveOptimal-[TorqueBicepsOptimal;zeros(size(Q1))]-[zeros(size(Q1));TorqueBicepsOptimal];
WorkActuator    =sum(sum(abs(TorqueActive.*Dq(:,Middle1:Middle2)),2))*Tres;
WorkDesire      =sum(sum(abs(TorqueDesire.*Dq(:,Middle1:Middle2)),2))*Tres;
Title=sprintf('%10s  %15s %15s'  ,   'Cost','Req. Work','Act. Work');
Result=sprintf('% 11.2f  % 11.2f % 15.2f'  ,   Cost,WorkDesire,WorkActuator);
disp('-------------------------------------------------------')
disp(Title)
disp(Result)


figure('name','Path (3DoF)')
    plot(XEF,YEF,'linewidth',2.5,'linestyle','-','color','g')
    hold on 
    plot(RPosVal(1,:),RPosVal(2,:),'linewidth',2,'linestyle','-.','color','b')
    legend('Desired','Static Path Planing')
    hold off
    axis equal

figure('name','Passive Torques')
    subplot(3,3,1)
    plot(Q1,TorqueDesQ1-TorqueBicepsOptimal(1,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q1,TorquePassiveOptimal(1,:),'linewidth',2,'linestyle','-.','color','r')
    hold off
    grid on
    xlabel('q_1')
    ylabel('u_r_1-u_b_1')
%     legend('Desired Torque','PassiveTorque')
    subplot(3,3,2)
    plot(Q1,TorqueActive(1,:),'linewidth',2)
    grid on
    xlabel('q_1')
    ylabel('u_a_1')
    title('Active Torque')
    
    subplot(3,3,4)
    plot(Q2,TorqueDesQ2-TorqueBicepsOptimal(1,:)-TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q2,TorquePassiveOptimal(2,:),'linewidth',2,'linestyle','-.','color','r')
    hold off
    grid on
    xlabel('q_2')
    ylabel('u_r_2-u_b_1-u_b_2')
%     legend('Desired Torque','PassiveTorque')
    subplot(3,3,5)
    plot(Q2,TorqueActive(2,:),'linewidth',2)
    grid on
    xlabel('q_2')
    ylabel('u_a_2')
    title('Active Torque')
    
    subplot(3,3,7)
    plot(Q3,TorqueDesQ3-TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q3,TorquePassiveOptimal(3,:),'linewidth',2,'linestyle','-.','color','r')
    hold off
    grid on
    xlabel('q_3')
    ylabel('u_r_3-u_b_2')
%     legend('Desired Torque','PassiveTorque')
    subplot(3,3,8)
    plot(Q3,TorqueActive(3,:),'linewidth',2)
    grid on
    xlabel('q_3')
    ylabel('u_a_3')
    title('Active Torque')
    
    subplot(3,3,3)
    plot(Qhat1,TorqueDesQ1+TorqueDesQ2-TorquePassiveOptimal(1,:)-TorquePassiveOptimal(2,:)-TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Qhat1,2*TorqueBicepsOptimal(1,:),'linewidth',2,'linestyle','-.','color','r')
    hold off
    grid on
    xlabel('q_1+q_2')
    ylabel('u_r_1+u_r_2-u_u_1-u_u_2-u_b_2')
%     legend('Desired Torque','PassiveTorque')

    subplot(3,3,6)
    plot(Qhat2,TorqueDesQ2+TorqueDesQ3-TorquePassiveOptimal(2,:)-TorquePassiveOptimal(3,:)-TorqueBicepsOptimal(1,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Qhat2,2*TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-.','color','r')
    hold off
    grid on
    xlabel('q_2+q_3')
    ylabel('u_r_2+u_r_3-u_u_2-u_u_3-u_b_1')

    
figure('name','Compare Torques')
    subplot(3,1,1)
    plot(Q1,TorqueDesQ1,'linewidth',2,'linestyle','-','color','b')
    xlabel('q_1')
    ylabel('u_1')
    hold on
    plot(Q1,TorqueActive(1,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    subplot(3,1,2)
    plot(Q2,TorqueDesQ2,'linewidth',2,'linestyle','-','color','b')
    xlabel('q_2')
    ylabel('u_2')
    hold on
    plot(Q2,TorqueActive(2,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    subplot(3,1,3)
    plot(Q3,TorqueDesQ3,'linewidth',2,'linestyle','-','color','b')
    xlabel('q_3')
    ylabel('u_3')
    hold on
    plot(Q3,TorqueActive(3,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    
    

    
figure('name','Time Torques')
    subplot(3,1,1)
    plot(Time,TorqueDesQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorquePassiveOptimal(1,:),'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueBicepsOptimal(1,:),'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    plot(Time,TorqueActive(1,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    Llg=legend('u_r_1','u_u_1','u_b_1','u_a_1');
    set(Llg,'orientation','horizontal')
    xlim([Time(1) Time(end)])
    subplot(3,1,2)
    plot(Time,TorqueDesQ2,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorquePassiveOptimal(2,:),'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueBicepsOptimal(1,:),'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    plot(Time,TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-.','Color',[0.75 .75 0.75])
    plot(Time,TorqueActive(2,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    Llg=legend('u_r_2','u_u_2','u_b_1','u_b_2','u_a_2');
    set(Llg,'orientation','horizontal')
    xlim([Time(1) Time(end)])
    subplot(3,1,3)
    plot(Time,TorqueDesQ3,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorquePassiveOptimal(3,:),'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-.','Color',[0.75 .75 0.75])
    plot(Time,TorqueActive(3,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    Llg=legend('u_r_3','u_u_3','u_b_2','u_a_3');
    set(Llg,'orientation','horizontal')
    xlim([Time(1) Time(end)])
    
%  toc
 
 

%% Optimization

Degree=[nn rQ rU rB];
% Initial=[Alpha_Q1 Alpha_Q2 Alpha_Q3];
% % 
%  Initial2=Initial;
%  Initial=x;

% WeightMatrix
% Weight=[ 10 1 1]';

% dynamic parameters
mL1=m;
mL2=m;
mL3=m;
LL1=L;
LL2=L;
LL3=L;

tic
MaxFunEvals_Data=5000*(rQ);
MaxIter_Data=1000;
TolFun_Data=1e-9;
TolX_Data=1e-9;
TolCon_Data=1e-7;
Algorithm='sqp';
Algorithm='interior-point';
Rand=5000*1e-6;
NewInit=Initial+Rand*(randn(1,3*(rQ+length(rQ))));

CostFun   = @(Alpha)CF3b_TorqueCost(Alpha,Time,Degree,Tres,Weight,Landa,SampleRate,g,mL1,mL2,mL3,LL1,LL2,LL3);
NonConstr = @(Alpha)CF3b_NonLinearConstraint(Alpha,Time,Tres,Degree,L,XEF,YEF,g,mL1,mL2,mL3,LL1,LL2,LL3);


[x,fval,exitflag,output,lambda,grad,hessian] = ...
    Op_FmisCon_SQP(CostFun,NonConstr,NewInit,MaxFunEvals_Data,MaxIter_Data,TolFun_Data,TolX_Data,TolCon_Data,Algorithm);




%% ShowTime
Degree=[nn rQ rU rB];

[TorqueDesire_X0,TorqueActive_X0,Q_X0,D1Q_X0,D2Q_X0,BetaOptimal_X0,ThetaOptimal_X0,IntU2_X0,IntUdq_X0,IntAbsUdq_X0,IntAbsUdqDesire_X0,CostSlopeD1Q_X0,CostSlopeD2Q_X0,CostParam_X0,RMSError_X0]=...    
                        ShowTime(Initial,Time,Tres,Degree,Weight,Landa,SampleRate,[],[],[],XEF,YEF,m,L,g,[],'DntShow','2Cycle','CostCc','Initial');
[TorqueDesire_Opt,TorqueActive_Opt,Q_Opt,D1Q_Opt,D2Q_Opt,BetaOptimal_Opt,ThetaOptimal_Opt,IntU2_Opt,IntUdq_Opt,IntAbsUdq_Opt,IntAbsUdqDesire_Opt,CostSlopeD1Q_Opt,CostSlopeD2Q_Opt,CostParam_Opt,RMSError_Opt]=...
                        ShowTime(x      ,Time,Tres,Degree,Weight,Landa,SampleRate,[],[],[],XEF,YEF,m,L,g,[],'DntShow','2Cycle','CostCc','Optimized');
                    

                    
TotalCost_X0  = (IntU2_X0    + Landa(1:2:3)*CostParam_X0 + Landa(2:2:4)*CostSlopeD2Q_X0 )*Tres;
TotalCost_Opt = (IntU2_Opt   + Landa(1:2:3)*CostParam_Opt+ Landa(2:2:4)*CostSlopeD2Q_Opt)*Tres;

IntAbsUdqDesire_X0=sum(sum(abs(TorqueDesire_X0.*D1Q_X0),2))*Tres;
IntAbsUdqDesire_Opt=sum(sum(abs(TorqueDesire_Opt.*D1Q_Opt),2))*Tres;
IntAbsUdqActive_X0=sum(sum(abs(TorqueActive_X0.*D1Q_X0),2))*Tres;
IntAbsUdqActive_Opt=sum(sum(abs(TorqueActive_Opt.*D1Q_Opt),2))*Tres;

DegreeStr=sprintf('\n   rQ:%3d\n   rU:%3d\n   Cost Type:  %s\n   Requlation Type:  %s\n ',rQ,rU, 'Cast 2b',SeletcStr{find(SelectLanda)});
Title=sprintf('%22s  % 11s %11s % 8s % 12s % 15s % 15s'  ,   'IntU2','C.S.(D2)','ParamReg','Total','RMS Err','Req. Work','Actuator Work');
Result_X0 =sprintf('%-11s %11.2e %11.2e %11.2e % 10.2e % 10.2e % 15.2e % 14.2e  ',   'Initial:',  IntU2_X0, sum(CostSlopeD2Q_X0) , sum(CostParam_X0)   , TotalCost_X0 , RMSError_X0, sum(IntAbsUdqDesire_X0), IntAbsUdqActive_X0);
Result_Opt=sprintf('%-11s %11.2e %11.2e %11.2e % 10.2e % 10.2e % 15.2e % 14.2e\n',   'Optimized:',IntU2_Opt,sum(CostSlopeD2Q_Opt), sum(CostParam_Opt)  , TotalCost_Opt, RMSError_Opt,sum(IntAbsUdqDesire_Opt),IntAbsUdqActive_Opt);
% display(output.message)
disp(DegreeStr)
disp(Title)
disp(Result_X0)
disp(Result_Opt)

TorqueActive=TorqueDesire-TorquePassiveOptimal-[TorqueBicepsOptimal;zeros(size(Q1))]-[zeros(size(Q1));TorqueBicepsOptimal];
IntAbsUdqActive=sum(sum(abs(TorqueActive.*Dq(:,Middle1:Middle2)),2))*Tres;
IntAbsUdqDesire=sum(sum(abs(TorqueDesire.*Dq(:,Middle1:Middle2)),2))*Tres;
IntAbsUdqActive=sum(sum(abs(TorqueActive.*Dq(:,Middle1:Middle2)),2))*Tres;
IntAbsUdqDesire=sum(sum(abs(TorqueDesire.*Dq(:,Middle1:Middle2)),2))*Tres;
% % RmsErrorX0=sqrt(sum(sum((rPosVal(:,Middle1:Middle2)-pos(:,Middle-1:end-1)).^2)*Tres/(time(end)-time(Middle))));
% % MaxErrorX0=max(sqrt(sum((rPosVal(:,Middle1:Middle2)-pos(:,Middle-1:end-1)).^2)));
% disp([RmsErrorX0 MaxErrorX0 ]*100)
disp([IntAbsUdqActive IntAbsUdqDesire ])

% toc


%% Scale and shift profile for schmit model

ThetaS=[];
ThetaStep=[];
tau=[];

for Joint=1:nn


    ThetaStep{Joint}=( (max(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)))  - min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)))) / (floor(size(Q_Opt,2)/2)));
    thetaS=min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2))) :ThetaStep{Joint}:max(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)));
    thetaScale=thetaS-(min(thetaS)+(max(thetaS)-min(thetaS))/2)+pi/2;
    ThetaS{Joint}=thetaScale;
    

    tau{Joint}=polyval(BetaOptimal_Opt((rU+1)*(Joint-1)+1: (Joint)*(rU+1)),thetaS);    
    

    figure('name',['Joint ', num2str( Joint)])
        subplot(2,1,1)
        plot(rad2deg( ThetaS{Joint}),tau{Joint})
        xlabel('ThetaS')
        ylabel('tau')
        grid on

        subplot(2,1,2)
        plot(ThetaS{Joint}(1:end-1),diff(tau{Joint})./diff(ThetaS{Joint}))
        xlabel('ThetaS')
        ylabel('Dtau')
        grid on
        
end

%%
% ParamA={[10000,0.02,.1,.14];
%         [10000,0.02,.1,.14];
%         [10000,0.02,.1,.14]};
ParamA={[5000,0.02,.1,.14];
        [5000,0.02,.1,.14];
        [5000,0.02,.1,.14]};
    
for Joint=1:nn
    K=ParamA{Joint}(1);
    R=ParamA{Joint}(2);
    l0=ParamA{Joint}(3);
    l0de=ParamA{Joint}(4);

    run('../../NonLinearSpring/RunMe');
end
%% Ga for find best Param
% 
% k0=2000;
% R0=1;
% q00=1;
% 
% nvars=3;
% lb=[50 5e-2 5e-2];
% PopInitRange=[lb; k0 R0 q00];
% PopulationSize=200;
% InitialPopulation=[k0*rand(PopulationSize,1) R0*rand(PopulationSize,1) q00*rand(PopulationSize,1)];
% 
% ParamA=[];
% for Joint=1:nn
%     disp(['Joint:',num2str(Joint)])
%     CostParam=@(Param)GA_CostParamNonLinearSpring(Param,ThetaStepscale{Joint},ThetaShiftScale{Joint},tauShiftScale{Joint});
% %     CostParam=@(Param)GA_CostParamNonLinearSpring(Param,ThetaStepMain{Joint},ThetaShiftMain{Joint},tauShiftScale{Joint});
% 
%     [ParamA{Joint},fval,exitflag,output,population,score] = ...
%         GA_FminCost(CostParam,nvars,lb,PopInitRange,PopulationSize,InitialPopulation);
%     disp(output.message)
% end
% %%
% % ParamA=[kg+600,Rg-.85,qg0+.25]
% for Joint=1:nn
%     k=ParamA{Joint}(1)*max(tauShiftMain{Joint});
%     R=ParamA{Joint}(2)*2;
%     q0=ParamA{Joint}(3);
% 
%     NonLinearSpring(ThetaStepscale{Joint},ThetaShiftScale{Joint},tauShiftScale{Joint}*max(tauShiftMain{Joint}),k,R,q0,.1,['Joint ',num2str(Joint)])
%     % NonLinearSpring(ThetaStep,ThetaS,tauShift,k,R,q0,.1)
% end
%% Control

InitState=([ Q_Opt(1:3,1)' , D1Q_Opt(1:3,1)', 0 0]);
Pos  =[[XEF   ,XEF(2:end)]   ;[YEF   ,YEF(2:end)]   ];
% DPos =[[DXEF  ,DXEF(2:end)]  ;[DYEF  ,DYEF(2:end)]  ];
% D2Pos=[[D2XEF ,D2XEF(2:end)] ;[D2YEF ,D2YEF(2:end)] ];
Time2=[Time , Time(2:end)+Time(end)];
% TorqueFF=[TorqueDesire_Opt,TorqueDesire_Opt(:,2:end)];
TorqueFF=[TorqueActive_Opt,TorqueActive_Opt(:,2:end)];
Qref=[Q_Opt, Q_Opt(:,2:end)];
DQref=[D1Q_Opt, D1Q_Opt(:,2:end)];
Mid=ceil(length(Time2)/1);

Kp=1*4000;
Kd=1*400;
% Kp=1*20000;
% Kd=1*400;


OdeOpt= odeset('RelTol',1e-8,'AbsTol',1e-8,'maxstep',1e-3,...
               'OutputFcn',@(t,Y,flag)TorqueCalculatorControl(t,Y,flag,g,L,m,Qref,DQref,Time2,Kp,Kd,BetaOptimal_Opt*1,rU));
[T,Q_ode] =...
    ode15s(@(t,Y)SirDyn(t,Y,g,L,m,Qref,DQref,Time2,Kp,Kd,BetaOptimal_Opt*1,rU),...
            Time2,InitState,OdeOpt);

QO1=Q_ode(:,1)';
QO2=Q_ode(:,2)';
QO3=Q_ode(:,3)';
RPosODE=L*[cos(QO1)+cos(QO1+QO2)+cos(QO1+QO2+QO3);
        sin(QO1)+sin(QO1+QO2)+sin(QO1+QO2+QO3)];

% EnergyABS_Ode=Q_ode(end,end-1)
% EnergyACT_Ode=Q_ode(end,end)

Torque_Ode=TorqueCalculatorControl([],[],'done');
TorqueActive_Ode=Torque_Ode(1:3,:);
TorquePassive_Ode=Torque_Ode(4:6,:);

EnergyABS=sum(sum(abs((TorqueActive_Ode+TorquePassive_Ode)'.*Q_ode(:,4:6)),2))*Tres;
EnergyACT=sum(sum(abs(TorqueActive_Ode'.*Q_ode(:,4:6)),2))*Tres;
disp(['Perfect: ',num2str(sum(IntAbsUdqDesire_Opt)*2),'  ,  ',num2str( IntAbsUdqActive_Opt*2)]);
disp(['PID:     ',num2str(EnergyABS),'  ,  ',num2str( EnergyACT)]);
disp(['Initial: ',num2str(sum(IntAbsUdqDesire_X0)*2),'  ,  ',num2str( IntAbsUdqActive_X0*2)]);

%%
figure
    plot(Time2,rad2deg(Qref(1,:)-Q_ode(:,1)'),'linewidth',2,'color','r')
    hold on
    plot(Time2,rad2deg(Qref(2,:)-Q_ode(:,2)'),'linewidth',2,'color','b')
    plot(Time2,rad2deg(Qref(3,:)-Q_ode(:,3)'),'linewidth',2,'color','g')
    set(gca,'fontsize',14,'YMinorGrid','on')    
    grid on
    xlabel('Time (s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10')
    ylabel('Tracking Error (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10')
    Le1=legend('e_1','e_2','e_3');
    set(Le1,'FontSize',12,'Orientation','horizontal')
    



figure
subplot(3,1,1)
    plot(Time2,Qref(1,:),'linewidth',2,'linestyle','--','color','r')
    hold on
    plot(Time2,Q_ode(:,1),'linewidth',2,'linestyle','--','color','b')
    legend('{q_1}^*','q_r')
   
subplot(3,1,2)
    plot(Time2,Qref(2,:))
    hold all
    plot(Time2,Q_ode(:,2))

subplot(3,1,3)
    plot(Time2,Qref(3,:))
    hold all
    plot(Time2,Q_ode(:,3))

figure
plot(Pos(1,:),Pos(2,:),'linewidth',4,'linestyle','-','color','b')
hold on
plot(RPosODE(1,1:Mid),RPosODE(2,1:Mid),'linewidth',3,'linestyle','--','color','r')
plot(RPosODE(1,Mid:end),RPosODE(2,Mid:end),'linewidth',2,'linestyle','--','color','g')
xlabel('x')
ylabel('y')
hold off
% axis equal
legend('Desired','dynamic path planing')

figure
subplot(3,1,1)
    plot(Time2,TorqueFF(1,:),'linewidth',2,'linestyle','--','color','r')
    hold on
    plot(Time2,TorqueActive_Ode(1,:),'linewidth',2,'linestyle','-','color','b')
    grid on
    set(gca,'fontsize',14,'YMinorGrid','on')
    L1=legend('Ideal Controller','PD Controller');
    set(L1,'FontSize',12,'Orientation','horizontal')
    ylabel('{ua}_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10')
    title('Actuator Torque')
    
    
subplot(3,1,2)
    plot(Time2,TorqueFF(2,:),'linewidth',2,'linestyle','--','color','r')
    hold on
    plot(Time2,TorqueActive_Ode(2,:),'linewidth',2,'linestyle','-','color','b')
    grid on
    set(gca,'fontsize',14,'YMinorGrid','on')
    L2=legend('Ideal Controller','PD Controller');
    ylabel('{ua}_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10')
    set(L2,'FontSize',12,'Orientation','horizontal')
%    title('Actuator Torque')
    
subplot(3,1,3)
    plot(Time2,TorqueFF(3,:),'linewidth',2,'linestyle','--','color','r')
    hold on
    plot(Time2,TorqueActive_Ode(3,:),'linewidth',2,'linestyle','-','color','b')
    grid on
    set(gca,'fontsize',14,'YMinorGrid','on')
    L3=legend('Ideal Controller','PD Controller');
    set(L3,'FontSize',12,'Orientation','horizontal')
    ylabel('{ua}_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10')
%     title('Actuator Torque')
    xlabel('Time (S)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10')
    
