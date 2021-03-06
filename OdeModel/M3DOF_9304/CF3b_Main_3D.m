function CF3b_Main_3D()



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
LL1=.20;
LL2=.56;
LL3=.51;
LL4=.20;
m_total=51;
mL1=1.1*LL1/(1.1*.4+1*LL2+.9*LL3+.8*LL4)*m_total;
mL2=1.0*LL2/(1.1*.4+1*LL2+.9*LL3+.8*LL4)*m_total;
mL3=0.9*LL3/(1.1*.4+1*LL2+.9*LL3+.8*LL4)*m_total;
mL4=0.8*LL4/(1.1*.4+1*LL2+.9*LL3+.8*LL4)*m_total;

% EF motion
A=.3;
 
% 3D Ellipose motion
f=0.5;
Tres=0.005;
Repeat=10;
time=0:Tres:Repeat/f;
Middle1=ceil((Repeat-2) *length(time)/Repeat);
Middle2=ceil((Repeat-1)*length(time)/Repeat);
xef= A/2*sin(2*pi*f*time)+.7;
yef= A*cos(2*pi*f*time)+.1;
zef=-A/4*cos(4*pi*f*time)+.5;

Dxef= (2*pi*f)*A/2*cos(2*pi*f*time);
Dyef=-(2*pi*f)*A*sin(2*pi*f*time);
Dzef= (4*pi*f)*A/4*sin(4*pi*f*time);

D2xef=-(2*pi*f)^2*A/2*sin(2*pi*f*time);
D2yef=-(2*pi*f)^2*A*cos(2*pi*f*time);
D2zef= (4*pi*f)^2*A/4*cos(2*pi*f*time);

q1=deg2rad(0);
q2=deg2rad(0);
q3=deg2rad(90);
q4=deg2rad(90);


plot3(xef,yef,zef)

pos=[xef;yef;zef];
Dpos=[Dxef;Dyef;Dzef];
D2pos=[D2xef;D2yef;D2zef];


% Joint Trajectories

for tt=1:length(time)-1
%     tt
    lPos=FK_RzRyRyRy_3D(q1(tt),q2(tt),q3(tt),q4(tt),LL1,LL2,LL3,LL4);

    cPos=pos(:,tt);
    dx=cPos-lPos;

    JJ=JACOBIAN_RzRyRyRy_3D(q1(tt),q2(tt),q3(tt),q4(tt),LL1,LL2,LL3,LL4);

    dq=JJ'*(JJ*JJ')^-1*dx;

    q1(tt+1)=q1(tt)+dq(1);
    q2(tt+1)=q2(tt)+dq(2);
    q3(tt+1)=q3(tt)+dq(3);
    q4(tt+1)=q4(tt)+dq(4);
    
end

for tt=1:length(time)-1
    rPosVal(:,tt)=FK_RzRyRyRy_3D(q1(tt),q2(tt),q3(tt),q4(tt),LL1,LL2,LL3,LL4);

end
    
% Required torque for desired path
q=[q1;q2;q3;q4];
Dq=differential(q,time,Tres);
D2q=differential(Dq,time,Tres);
torqueDesire=TorqueCalculator_4R_3D(D2q,Dq,q,g,mL1,mL2,mL3,mL4,LL1,LL2,LL3,LL4);
%%
% show    
% close all
figure('name','Path (3DoF)')
    plot3(xef(Middle1:Middle2),yef(Middle1:Middle2),zef(Middle1:Middle2),'linewidth',2.5,'linestyle','-','color','g')
    hold on 
    plot3(rPosVal(1,Middle1:Middle2),rPosVal(2,Middle1:Middle2),rPosVal(3,Middle1:Middle2),'linewidth',2,'linestyle','-.','color','b')
    legend('Desired','Static Path Planing')
    hold off
    axis equal

figure('name','Joint Trajectories')
    subplot(4,1,1)
    plot(time(Middle1:Middle2),q1(Middle1:Middle2))
    grid on
    subplot(4,1,2)
    plot(time(Middle1:Middle2),q2(Middle1:Middle2))
    grid on
    subplot(4,1,3)
    plot(time(Middle1:Middle2),q3(Middle1:Middle2))
    grid on
    subplot(4,1,4)
    plot(time(Middle1:Middle2),q4(Middle1:Middle2))
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
    legend('p_1','p_2','p_3','p_4')
    xlabel('time (s)','fontsize',14,'FontWeight','bold')
    ylabel('power','fontsize',14,'FontWeight','bold')
    grid on
    
    
figure('name','Desired Torque vs Time')
    plot(time(Middle1:Middle2),torqueDesire(:,Middle1:Middle2),'linewidth',2)
    set(gca,'fontsize',12,'FontWeight','bold')
    legend('u_1','u_2','u_3','u_4')
    xlabel('time (s)','fontsize',14,'FontWeight','bold')
    ylabel('u','fontsize',14,'FontWeight','bold')
    grid on
    


figure('name','Desired Torque vs Angle')
    subplot(1,4,1)
    plot(rad2deg(q1(Middle1:Middle2)),torqueDesire(1,Middle1:Middle2),...
        'Color',[0.87058824300766 0.490196079015732 0],'linewidth',2)
    title('Initial Required Torque-Angle Profile','FontSize',16);
    xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
    subplot(1,4,2)
    plot(rad2deg(InRangeShifter(q2(Middle1:Middle2))),torqueDesire(2,Middle1:Middle2),...
            'Color',[0.87058824300766 0.490196079015732 0],'linewidth',2)
    xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
    subplot(1,4,3)
    plot(( rad2deg(InRangeShifter(q3(Middle1:Middle2)))),torqueDesire(3,Middle1:Middle2),...
        'Color',[0.87058824300766 0.490196079015732 0],'linewidth',2)
    xlabel('q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
    subplot(1,4,4)
    plot(( rad2deg(InRangeShifter(q4(Middle1:Middle2)))),torqueDesire(4,Middle1:Middle2),...
        'Color',[0.87058824300766 0.490196079015732 0],'linewidth',2)
    xlabel('q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
Y=[q1(Middle1:end)',q2(Middle1:end)',q3(Middle1:end)' , q4(Middle1:end)' ];
AnimBot4DOF_4R_3D(time(Middle1:end),Y,LL1,LL2,LL3,LL4);

IntAbsUdqDesire=sum(sum(abs(torqueDesire(:,Middle1:Middle2).*Dq(:,Middle1:Middle2)),2))*Tres
IntAbsUdqDesire=sum(sum((torqueDesire(:,Middle1:Middle2).*torqueDesire(:,Middle1:Middle2)),2))*Tres

%%  Generate Initial value for Optimization
% close all

tic
% DoF system
nn=4; % number of joints
% DoF of Optimization 
rQ=10; % Degree of joint trajectory
rM=4; % Degree of passive torque
rB=2;
% WeightMatrix
Weight=[1 3 2 1]';
Weight=Weight*1;
SampleRate=25;

% Landa for [DQ  D2q ]
Landa=[1e-1 1e-3 1e-1 1e-3]*1;   % [landa_1* Beta'*Beta    Landa_2*(D2Q*Beta)'*(D2Q*Beta)
                           %  landa_3* Theta'*Theta  Landa_4*(D2Qhat*Theta)'*(D2Qhat*Theta) ]          
%  Landa=[1e-6 1e-7 1e-6 1e-7]*1; % for linear

%%
Time=time(Middle1:Middle2)-time(Middle1);
Q1=q1(Middle1:Middle2);
Q2=q2(Middle1:Middle2);
Q3=q3(Middle1:Middle2);
Q4=q4(Middle1:Middle2);
 
QJ=[Q1;Q2;Q3;Q4];

Qhat1=zeros(size(Q1));
Qhat2=Q2+Q3;
Qhat3=Q3+Q4;
QhatJ=[Qhat1;Qhat2;Qhat3];

XEF=xef(Middle1:Middle2);
YEF=yef(Middle1:Middle2);
ZEF=zef(Middle1:Middle2);
DXEF=Dxef(Middle1:Middle2);
DYEF=Dyef(Middle1:Middle2);
DZEF=Dzef(Middle1:Middle2);
D2XEF=D2xef(Middle1:Middle2);
D2YEF=D2yef(Middle1:Middle2);
D2ZEF=D2zef(Middle1:Middle2);
TorqueDesQ1=torqueDesire(1,Middle1:Middle2);
TorqueDesQ2=torqueDesire(2,Middle1:Middle2);
TorqueDesQ3=torqueDesire(3,Middle1:Middle2);
TorqueDesQ4=torqueDesire(4,Middle1:Middle2);
TorqueDesire=[TorqueDesQ1;TorqueDesQ2;TorqueDesQ3;TorqueDesQ4];


% [Alpha_Q1,BezireCoef_q1]= BezierCoeffinet([Time Time(2:end)+Time(1) ],[Q1 Q1(2:end)],rQ);
% [Alpha_Q2,BezireCoef_q2]= BezierCoeffinet(Time,Q2,rQ);
% [Alpha_Q3,BezireCoef_q3]= BezierCoeffinet(Time,Q3,rQ);
% [Alpha_Q4,BezireCoef_q4]= BezierCoeffinet(Time,Q4,rQ);
[Alpha_Q1,BezireCoef_q1]= BezierCoeffinet([Time Time(2:end)+Time(1) ],[Q1 Q1(2:end)],rQ);
[Alpha_Q2,BezireCoef_q2]= BezierCoeffinet([Time Time(2:end)+Time(1) ],[Q2 Q2(2:end)],rQ);
[Alpha_Q3,BezireCoef_q3]= BezierCoeffinet([Time Time(2:end)+Time(1) ],[Q3 Q3(2:end)],rQ);
[Alpha_Q4,BezireCoef_q4]= BezierCoeffinet([Time Time(2:end)+Time(1) ],[Q4 Q4(2:end)],rQ);

Q1val=polyval(Alpha_Q1,Time);
Q2val=polyval(Alpha_Q2,Time);
Q3val=polyval(Alpha_Q3,Time);
Q4val=polyval(Alpha_Q4,Time);
QVal=[Q1val;Q2val;Q3val;Q4val];


RPosVal=FK_RzRyRyRy_3D(Q1val,Q2val,Q3val,Q4val,LL1,LL2,LL3,LL4);

    
[BetaOptimal,ThetaOptimal,CostActuation,CostD2Q,CostParaReg,TorqueMonoOptimal,TorqueBicepsOptimal]= ...
            OptimalParam_4R_3D(QJ,QhatJ,TorqueDesire,nn,rM,rB,Landa,Weight,SampleRate);

CostRequired=0;
 for i=1:nn
    CostRequired=CostRequired + ...
          Weight(i)*( 1/2*(TorqueDesire(i,:)')'*TorqueDesire(i,:)');
 end
CostRequired*Tres;
           
        
Cost=(CostActuation+ Landa(1:2:4)*CostParaReg+Landa(2:2:4)*CostD2Q)*Tres;
ActCost=(CostActuation)*Tres;


TorqueActive=TorqueDesire-TorqueMonoOptimal-[TorqueBicepsOptimal;zeros(size(Q1))]-[zeros(size(Q1));TorqueBicepsOptimal];
WorkActuator    =sum(sum(abs(TorqueActive.*Dq(:,Middle1:Middle2)),2))*Tres;
WorkDesire      =sum(sum(abs(TorqueDesire.*Dq(:,Middle1:Middle2)),2))*Tres;
Title=sprintf('%11s  %14s  %15s %15s'  ,   'Cost','Act. Cost','Req. Work','Act. Work');
Result=sprintf('% 11.2f  % 11.2f  % 15.2f % 15.2f'  ,   Cost,ActCost,WorkDesire,WorkActuator);
disp('-------------------------------------------------------')
disp(Title)
disp(Result)


figure('name','Path (3DoF)')
    plot3(XEF,YEF,ZEF,'linewidth',2.5,'linestyle','-','color','g')
    hold on 
    plot3(RPosVal(1,:),RPosVal(2,:),RPosVal(3,:),'linewidth',2,'linestyle','-.','color','b')
    legend('Desired','Static Path Planing')
    hold off
    axis equal

figure('name','Passive Torques')
    subplot(4,3,1)
    plot(Q1,TorqueDesQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q1,TorqueMonoOptimal(1,:),'linewidth',2,'linestyle','-.','color','r')
    hold off
    grid on
    xlabel('q_1')
    ylabel('u_r_1')
%     legend('Desired Torque','PassiveTorque')
    subplot(4,3,2)
    plot(Q1,TorqueActive(1,:),'linewidth',2)
    grid on
    xlabel('q_1')
    ylabel('u_a_1')
    title('Active Torque')
    
    subplot(4,3,4)
    plot(Q2,TorqueDesQ2-TorqueBicepsOptimal(1,:)-TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q2,TorqueMonoOptimal(2,:),'linewidth',2,'linestyle','-.','color','r')
    hold off
    grid on
    xlabel('q_2')
    ylabel('u_r_2-u_b_2')
%     legend('Desired Torque','PassiveTorque')
    subplot(4,3,5)
    plot(Q2,TorqueActive(2,:),'linewidth',2)
    grid on
    xlabel('q_2')
    ylabel('u_a_2')
    title('Active Torque')
    
    subplot(4,3,7)
    plot(Q3,TorqueDesQ3-TorqueBicepsOptimal(2,:)-TorqueBicepsOptimal(3,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q3,TorqueMonoOptimal(3,:),'linewidth',2,'linestyle','-.','color','r')
    hold off
    grid on
    xlabel('q_3')
    ylabel('u_r_3-u_b_2-u_b_3')
%     legend('Desired Torque','PassiveTorque')
    subplot(4,3,8)
    plot(Q3,TorqueActive(3,:),'linewidth',2)
    grid on
    xlabel('q_3')
    ylabel('u_a_3')
    title('Active Torque')
    
    subplot(4,3,10)
    plot(Q4,TorqueDesQ4-TorqueBicepsOptimal(3,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q4,TorqueMonoOptimal(4,:),'linewidth',2,'linestyle','-.','color','r')
    hold off
    grid on
    xlabel('q_4')
    ylabel('u_r_4-u_b_3')
%     legend('Desired Torque','PassiveTorque')
    subplot(4,3,11)
    plot(Q4,TorqueActive(4,:),'linewidth',2)
    grid on
    xlabel('q_4')
    ylabel('u_a_4')
    title('Active Torque')
    
    subplot(4,3,6)
    plot(Qhat2,TorqueDesQ2+TorqueDesQ3-TorqueMonoOptimal(2,:)-TorqueMonoOptimal(3,:)-TorqueBicepsOptimal(3,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Qhat2,2*TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-.','color','r')
    hold off
    grid on
    xlabel('q_2+q_3')
    ylabel('u_r_2+u_r_3-u_m_2-u_m_3-u_b_3')
%     legend('Desired Torque','PassiveTorque')

    subplot(4,3,12)
    plot(Qhat3,TorqueDesQ3+TorqueDesQ4-TorqueMonoOptimal(3,:)-TorqueMonoOptimal(4,:)-TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Qhat3,2*TorqueBicepsOptimal(3,:),'linewidth',2,'linestyle','-.','color','r')
    hold off
    grid on
    xlabel('q_3+q_4')
    ylabel('u_r_3+u_r_4-u_m_3-u_m_4-u_b_2')

    
figure('name','Compare Torques')
    subplot(4,1,1)
    plot(Q1,TorqueDesQ1,'linewidth',2,'linestyle','-','color','b')
    xlabel('q_1')
    ylabel('u_1')
    hold on
    plot(Q1,TorqueActive(1,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    subplot(4,1,2)
    plot(Q2,TorqueDesQ2,'linewidth',2,'linestyle','-','color','b')
    xlabel('q_2')
    ylabel('u_2')
    hold on
    plot(Q2,TorqueActive(2,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    subplot(4,1,3)
    plot(Q3,TorqueDesQ3,'linewidth',2,'linestyle','-','color','b')
    xlabel('q_3')
    ylabel('u_3')
    hold on
    plot(Q3,TorqueActive(3,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    subplot(4,1,4)
    plot(Q4,TorqueDesQ4,'linewidth',2,'linestyle','-','color','b')
    xlabel('q_4')
    ylabel('u_4')
    hold on
    plot(Q4,TorqueActive(4,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    
    

    
figure('name','Time Torques')
    subplot(4,1,1)
    plot(Time,TorqueDesQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorqueMonoOptimal(1,:),'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueActive(1,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    Llg=legend('u_r_1','u_m_1','u_a_1');
    set(Llg,'orientation','horizontal')
    xlim([Time(1) Time(end)])
    subplot(4,1,2)
    plot(Time,TorqueDesQ2,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorqueMonoOptimal(2,:),'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    plot(Time,TorqueActive(2,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    Llg=legend('u_r_2','u_m_2','u_b_2','u_a_2');
    set(Llg,'orientation','horizontal')
    xlim([Time(1) Time(end)])
    subplot(4,1,3)
    plot(Time,TorqueDesQ3,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorqueMonoOptimal(3,:),'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
    plot(Time,TorqueBicepsOptimal(3,:),'linewidth',2,'linestyle','-.','Color',[0.75 .75 0.75])
    plot(Time,TorqueActive(3,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    Llg=legend('u_r_3','u_m_3','u_b_2','u_b_3','u_a_3');
    set(Llg,'orientation','horizontal')
    xlim([Time(1) Time(end)])
    subplot(4,1,4)
    plot(Time,TorqueDesQ4,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Time,TorqueMonoOptimal(4,:),'linewidth',2,'linestyle','-.','color','g')
    plot(Time,TorqueBicepsOptimal(3,:),'linewidth',2,'linestyle','-.','Color',[0.75 .75 0.75])
    plot(Time,TorqueActive(4,:),'linewidth',2,'linestyle','-','color','r')
    hold off
    grid on
    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('u_4 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    Llg=legend('u_r_4','u_m_4','u_b_3','u_a_4');
    set(Llg,'orientation','horizontal')
    xlim([Time(1) Time(end)])
    
%  toc
 
 

%% Optimization

Degree=[nn rQ rM rB];
% Initial=[Alpha_Q1 Alpha_Q2 Alpha_Q3 Alpha_Q4];
% Initial=[Alpha_Q2 Alpha_Q3 Alpha_Q4];

QLimit=deg2rad([-170,170;-96,135;-5,156;-20,110]);
DQLimit=deg2rad([-300,300;-225,225;-225,225;-300,300]);
% % 
%  Initial2=Initial;
%  Initial=x;

% WeightMatrix
% Weight=[ 10 1 1]';

% dynamic parameters
tic
MaxFunEvals_Data=5000*(rQ);
MaxIter_Data=1000;
TolFun_Data=1e-9;
TolX_Data=1e-9;
TolCon_Data=1e-8;
Algorithm='sqp';
Algorithm='interior-point';
Rand=5000*1e-10;
NewInit=Initial+Rand*(randn(1,4*(rQ+length(rQ))));

CostFun   = @(Alpha)CF3b_TorqueCost_4R_3D(Alpha,Alpha_Q1,Time,Degree,Tres,Weight,Landa,SampleRate,g,mL1,mL2,mL3,mL4,LL1,LL2,LL3,LL4);
% NonConstr = @(Alpha)CF3b_NonLinearConstraint(Alpha,Time,Tres,Degree,L,XEF,YEF,g,mL1,mL2,mL3,LL1,LL2,LL3);
NonConstr = @(Alpha)CF3b_NonLinearConstraint_WithJointLimit_4R_3D(Alpha,Alpha_Q1,Time,Tres,Degree,XEF,YEF,ZEF,g,mL1,mL2,mL3,mL4,LL1,LL2,LL3,LL4,QLimit,DQLimit,SampleRate);


[x,fval,exitflag,output,lambda,grad,hessian] = ...
    Op_FmisCon_SQP(CostFun,NonConstr,NewInit,MaxFunEvals_Data,MaxIter_Data,TolFun_Data,TolX_Data,TolCon_Data,Algorithm);




%% ShowTime
Degree=[nn rQ rM rB];
Initial=[Alpha_Q1 Alpha_Q2 Alpha_Q3 Alpha_Q4];

[TorqueDesire_X0,TorqueActive_X0,TorqueMonoOptimal_X0,TorqueBicepsOptimal_X0,Q_X0,D1Q_X0,D2Q_X0,BetaOptimal_X0,ThetaOptimal_X0,IntU2_X0,IntUdq_X0,IntAbsUdq_X0,IntAbsUdqDesire_X0,CostSlopeD1Q_X0,CostSlopeD2Q_X0,CostParam_X0,RMSError_X0]=...    
                        ShowTime_4R_3D([ Initial],Time,Tres,Degree,Weight,Landa,SampleRate,[],[],[],XEF,YEF,ZEF,mL1,mL2,mL3,mL4,LL1,LL2,LL3,LL4,g,[],'DntShow','2Cycle','CostCc','Initial');
[TorqueDesire_Opt,TorqueActive_Opt,TorqueMonoOptimal_Opt,TorqueBicepsOptimal_Opt,Q_Opt,D1Q_Opt,D2Q_Opt,BetaOptimal_Opt,ThetaOptimal_Opt,IntU2_Opt,IntUdq_Opt,IntAbsUdq_Opt,IntAbsUdqDesire_Opt,CostSlopeD1Q_Opt,CostSlopeD2Q_Opt,CostParam_Opt,RMSError_Opt]=...
                       ShowTime_4R_3D([ x]      ,Time,Tres,Degree,Weight,Landa,SampleRate,[],[],[],XEF,YEF,ZEF,mL1,mL2,mL3,mL4,LL1,LL2,LL3,LL4,g,[],'Show','2Cycle','CostCc','Optimized');
                    

                    
%%
TotalCost_X0  = (IntU2_X0    + Landa(1:2:3)*CostParam_X0 + Landa(2:2:4)*CostSlopeD2Q_X0 )*Tres;
TotalCost_Opt = (IntU2_Opt   + Landa(1:2:3)*CostParam_Opt+ Landa(2:2:4)*CostSlopeD2Q_Opt)*Tres;

IntAbsUdqDesire_X0=sum(sum(abs(TorqueDesire_X0.*D1Q_X0),2))*Tres;
IntAbsUdqDesire_Opt=sum(sum(abs(TorqueDesire_Opt.*D1Q_Opt),2))*Tres;
IntAbsUdqActive_X0=sum(sum(abs(TorqueActive_X0.*D1Q_X0),2))*Tres;
IntAbsUdqActive_Opt=sum(sum(abs(TorqueActive_Opt.*D1Q_Opt),2))*Tres;

DegreeStr=sprintf('\n   rQ:%3d\n   rM:%3d\n   rB:%3d\n   Cost Type:  %s\n   Requlation Type:  %s\n ',rQ,rM,rB, 'Cast Cc','D2Q');
Title=sprintf('%22s  % 11s %11s % 8s % 12s % 15s % 15s'  ,   'IntU2','C.S.(D2)','ParamReg','Total','RMS Err','Req. Work','Actuator Work');
Result_X0 =sprintf('%-11s %11.2e %11.2e %11.2e % 10.2e % 10.2e % 15.2e % 14.2e  ',   'Initial:',  IntU2_X0*Tres, sum(CostSlopeD2Q_X0)*Tres , sum(CostParam_X0)*Tres   , TotalCost_X0 , RMSError_X0, sum(IntAbsUdqDesire_X0), IntAbsUdqActive_X0);
Result_Opt=sprintf('%-11s %11.2e %11.2e %11.2e % 10.2e % 10.2e % 15.2e % 14.2e\n',   'Optimized:',IntU2_Opt*Tres,sum(CostSlopeD2Q_Opt)*Tres, sum(CostParam_Opt)*Tres  , TotalCost_Opt, RMSError_Opt,sum(IntAbsUdqDesire_Opt),IntAbsUdqActive_Opt);
% display(output.message)
disp(DegreeStr)
disp(Title)
disp(Result_X0)
disp(Result_Opt)

TorqueActive=TorqueDesire-TorqueMonoOptimal-[TorqueBicepsOptimal;zeros(size(Q1))]-[zeros(size(Q1));TorqueBicepsOptimal];
IntAbsUdqActive=sum(sum(abs(TorqueActive.*Dq(:,Middle1:Middle2)),2))*Tres;
IntAbsUdqDesire=sum(sum(abs(TorqueDesire.*Dq(:,Middle1:Middle2)),2))*Tres;
IntAbsUdqActive=sum(sum(abs(TorqueActive.*Dq(:,Middle1:Middle2)),2))*Tres;
IntAbsUdqDesire=sum(sum(abs(TorqueDesire.*Dq(:,Middle1:Middle2)),2))*Tres;
% % RmsErrorX0=sqrt(sum(sum((rPosVal(:,Middle1:Middle2)-pos(:,Middle-1:end-1)).^2)*Tres/(time(end)-time(Middle))));
% % MaxErrorX0=max(sqrt(sum((rPosVal(:,Middle1:Middle2)-pos(:,Middle-1:end-1)).^2)));
% disp([RmsErrorX0 MaxErrorX0 ]*100)
disp([IntAbsUdqActive IntAbsUdqDesire ])

% toc
disp(sum(   max(D1Q_Opt')>DQLimit(:,2)' +min(D1Q_Opt')<DQLimit(:,1)'+max(Q_Opt')>QLimit(:,2)'+min(Q_Opt')<QLimit(:,1)' ) )

%% Scale profile for Nonlinear Spring

ThetaS=[];
ThetaStep=[];
tau=[];

for Joint=1:3

    ThetaStep{Joint}=( (max(Q_Opt(Joint, 1: floor(size(Q_Opt,2))))  - min(Q_Opt(Joint, 1: floor(size(Q_Opt,2))))) / (floor(size(Q_Opt,2)/2)));
    thetaS=min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)))) :ThetaStep{Joint}:max(Q_Opt(Joint, 1: floor(size(Q_Opt,2))));
    thetaScale=thetaS-(min(thetaS)+(max(thetaS)-min(thetaS))/2)+pi/2;
    ThetaS{Joint}=thetaScale;
    

    tau{Joint}=polyval(BetaOptimal_Opt((rM+1)*(Joint-1)+1: (Joint)*(rM+1)),thetaS);    
    

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

%% Desing Cam
ParamA={[10000,0.02,.1,.14];
        [10000,0.02,.1,.14];
        [10000,0.02,.1,.14]};
    
for Joint=1:3
    K=ParamA{Joint}(1);
    R=ParamA{Joint}(2);
    l0=ParamA{Joint}(3);
    l0de=ParamA{Joint}(4);

    run('../../NonLinearSpring/RunMe');
end
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
Kp=1*20000000;
Kd=1*400000;


OdeOpt= odeset('RelTol',1e-4,'AbsTol',1e-4,'maxstep',1e-2,...
               'OutputFcn',@(t,Y,flag)TorqueCalculatorControlBi(t,Y,flag,g,L,m,Qref,DQref,Time2,Kp,Kd,x,rM,rB));
[T,Q_ode] =...
    ode15s(@(t,Y)SirDynBi(t,Y,g,L,m,Qref,DQref,Time2,Kp,Kd,x,rM,0),...
            Time2,InitState,OdeOpt);

QO1=Q_ode(:,1)';
QO2=Q_ode(:,2)';
QO3=Q_ode(:,3)';
RPosODE=L*[cos(QO1)+cos(QO1+QO2)+cos(QO1+QO2+QO3);
        sin(QO1)+sin(QO1+QO2)+sin(QO1+QO2+QO3)];

% EnergyABS_Ode=Q_ode(end,end-1)
% EnergyACT_Ode=Q_ode(end,end)

Torque_Ode=TorqueCalculatorControlBi([],[],'done');
TorqueActive_Ode=Torque_Ode(1:3,:);
TorquePassiveUni_Ode=Torque_Ode(4:6,:);
TorquePassiveBi_Ode=Torque_Ode(7:9,:);

EnergyABS=sum(sum(abs((TorqueActive_Ode+TorquePassiveUni_Ode+TorquePa)'.*Q_ode(:,4:6)),2))*Tres;
EnergyACT=sum(sum(abs(TorqueActive_Ode'.*Q_ode(:,4:6)),2))*Tres;
disp(['Perfect: ',num2str(sum(IntAbsUdqDesire_Opt)*2),'  ,  ',num2str( IntAbsUdqActive_Opt*2)]);
disp(['PID:     ',num2str(EnergyABS),'  ,  ',num2str( EnergyACT)]);
disp(['Initial: ',num2str(sum(IntAbsUdqDesire_X0)*2),'  ,  ',num2str( IntAbsUdqActive_X0*2)]);

%% ShowTime for Control
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
    LL1=legend('Ideal Controller','PD Controller');
    set(LL1,'FontSize',12,'Orientation','horizontal')
    ylabel('{ua}_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10')
    title('Actuator Torque')
    
    
subplot(3,1,2)
    plot(Time2,TorqueFF(2,:),'linewidth',2,'linestyle','--','color','r')
    hold on
    plot(Time2,TorqueActive_Ode(2,:),'linewidth',2,'linestyle','-','color','b')
    grid on
    set(gca,'fontsize',14,'YMinorGrid','on')
    LL2=legend('Ideal Controller','PD Controller');
    ylabel('{ua}_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10')
    set(LL2,'FontSize',12,'Orientation','horizontal')
%    title('Actuator Torque')
    
subplot(3,1,3)
    plot(Time2,TorqueFF(3,:),'linewidth',2,'linestyle','--','color','r')
    hold on
    plot(Time2,TorqueActive_Ode(3,:),'linewidth',2,'linestyle','-','color','b')
    grid on
    set(gca,'fontsize',14,'YMinorGrid','on')
    LL3=legend('Ideal Controller','PD Controller');
    set(LL3,'FontSize',12,'Orientation','horizontal')
    ylabel('{ua}_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10')
%     title('Actuator Torque')
    xlabel('Time (S)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10')
    
