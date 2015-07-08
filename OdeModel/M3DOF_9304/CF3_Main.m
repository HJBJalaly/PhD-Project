function CF3_Main()



%% Initialization

clear
close all
home
rand('twister', sum(100*clock));

%% create a sample motion with 2 active and 1 passive joints


% OdeOpt= odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,7));
% InitState=deg2rad([45 0  0 ,0 0 0, 0]);
% 
% [T,Y] = ode15s(@(t,Y)SirDynSample(t,Y,g,L,m,KKmat), time,InitState,OdeOpt);
% 
% AnimBot3DOF(T,Y,L);
% 
% EnergyABS=Y(end,end)
% 
% q1=Y(:,1)';
% q2=Y(:,2)';
% q3=Y(:,3)';
% RPos=L*[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
%         sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];
% plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','--','color','r')
% xlabel('x')
% ylabel('y')
% hold off
% axis equal
% legend('Desired','Static path planing','dynamic path planing')
% 
% figure
% subplot(3,1,1)
% plot(T,q1)
% subplot(3,1,2)
% plot(T,q2)
% subplot(3,1,3)
% plot(T,q3)

%% create a sample passive motion with 3 passive joints
home
clear
close all

% envirment
g=9.81*0;

% parameters of robot
m=1;
L=1;

KKmat=[50;500;100];

% Desired motion
f=1;
A=1;
Tres=0.005;
time=0:Tres:6/f;
q0=deg2rad([0 -0 01]);

OdeOpt= odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,7));
InitState=[q0 ,0 0 0, 0];

[T,Y] = ode15s(@(t,Y)SirDynSample(t,Y,g,L,m,KKmat), time,InitState,OdeOpt);

EnergyABS=Y(end,end)

q1p=Y(:,1)';
q2p=Y(:,2)';
q3p=Y(:,3)';
xef=L*(cos(q1p)+cos(q1p+q2p)+cos(q1p+q2p+q3p));
yef=L*(sin(q1p)+sin(q1p+q2p)+sin(q1p+q2p+q3p));


pos=[xef;yef];


% Joint Trajectories
q1=q0(1);
q2=q0(2);
q3=q0(3);

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

RPosVal=L*[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
        sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];

    
% Required torque for desired path
q=[q1;q2;q3];
Dq=differential(q,time,Tres);
D2q=differential(Dq,time,Tres);
torqueDesire=TorqueCalculator(D2q,Dq,q,g,m,m,m,L,L,L);

Middle=ceil(length(time)/2);

% show    
close all
figure('name','Path (3DoF)')
    plot(xef,yef,'linewidth',2.5,'linestyle','-','color','g')
    hold on 
    plot(RPosVal(1,Middle:end),RPosVal(2,Middle:end),'linewidth',2,'linestyle','-.','color','b')
    legend('Desired','Static Path Planing')
    hold off
    axis equal

figure('name','Joint Trajectories')
    subplot(3,1,1)
    plot(time(Middle:end),q1(Middle:end))
    grid on
    hold all
    plot(T(Middle:end),q1p(Middle:end))
    hold off
    subplot(3,1,2)
    plot(time(Middle:end),q2(Middle:end))
    grid on
    hold all
    plot(T(Middle:end),q2p(Middle:end))
    hold off
    subplot(3,1,3)
    plot(time(Middle:end),q3(Middle:end))
    grid on
    hold all
    plot(T(Middle:end),q3p(Middle:end))
    hold off
    
figure('name','Torque vs Time')
    plot(time(Middle:end),torqueDesire(:,Middle:end))
    legend('\tau_1','\tau_2','\tau_3')
    xlabel('time (s)')
    ylabel('\tau')
    grid on
    


figure('name','Torque vs Angle')
    subplot(3,1,1)
    plot(q1(Middle:end),torqueDesire(1,Middle:end))
    grid on
    subplot(3,1,2)
    plot(q2(Middle:end),torqueDesire(2,Middle:end))
    grid on
    subplot(3,1,3)
    plot(q3(Middle:end),torqueDesire(3,Middle:end))
    grid on
    
Y=[q1(1:end)',q2(1:end)',q3(1:end)'];
AnimBot3DOF(time(1:end),Y,L);



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
f=1;
A=.75;
phi=pi/2;
Tres=0.002;
time=0:Tres:10/f;

% % Circle motion
% f=0.5;
% time=0:Tres:10/f;
% xef=A*cos(2*pi*f*time)+0;
% yef=A*sin(2*pi*f*time)+2;
% Dxef=-(2*pi*f)*A*sin(2*pi*f*time);
% Dyef= (2*pi*f)*A*cos(2*pi*f*time);
% D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time);
% D2yef=-(2*pi*f)^2*A*sin(2*pi*f*time);
% q1=deg2rad(-60);
% q2=deg2rad(-88.031);
% q3=deg2rad(11.317);

% Line motion: Horizontal
xef=A*cos(2*pi*f*time+phi);
yef=2.5*ones(size(time));
Dxef=-(2*pi*f)*A*sin(2*pi*f*time+phi);
Dyef= 2*zeros(size(time));
D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time+phi);
D2yef= 2*zeros(size(time));
q1=deg2rad( 42.7023); % 41.7023 deg for y=2.5m  ,  22.2 deg for y=2.0m
q2=deg2rad(19.2663); % 18.2663 deg for y=2.5m  ,  29.468 deg for y=2.0m
q3=deg2rad(80.3377); % 44.3377 deg for y=2.5m  ,  71.431 deg for y=2.0m


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

RPosVal=L*[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
        sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];

    
% Required torque for desired path
q=[q1;q2;q3];
Dq=differential(q,time,Tres);
D2q=differential(Dq,time,Tres);
torqueDesire=TorqueCalculator(D2q,Dq,q,g,m,m,m,L,L,L);

Middle=1;

% show    
% close all
Middle=ceil(9*length(time)/10);
figure('name','Path (3DoF)')
    plot(xef(Middle:end),yef(Middle:end),'linewidth',2.5,'linestyle','-','color','g')
    hold on 
    plot(RPosVal(1,Middle:end),RPosVal(2,Middle:end),'linewidth',2,'linestyle','-.','color','b')
    legend('Desired','Static Path Planing')
    hold off
    axis equal

figure('name','Joint Trajectories')
    subplot(3,1,1)
    plot(time(Middle:end),q1(Middle:end))
    grid on
    subplot(3,1,2)
    plot(time(Middle:end),q2(Middle:end))
    grid on
    subplot(3,1,3)
    plot(time(Middle:end),q3(Middle:end))
    grid on
    

figure('name','Desired Power vs Time')
    subplot(3,1,1)
    plot(time(Middle:end),torqueDesire(1,Middle:end).*Dq(1,Middle:end))
    grid on
    subplot(3,1,2)
    plot(time(Middle:end),torqueDesire(2,Middle:end).*Dq(2,Middle:end))
    grid on
    subplot(3,1,3)
    plot(time(Middle:end),torqueDesire(3,Middle:end).*Dq(3,Middle:end))
    grid on
    
    
figure('name','Desired Torque vs Time')
    plot(time(Middle:end),torqueDesire(:,Middle:end))
    legend('\tau_1','\tau_2','\tau_3')
    xlabel('time (s)')
    ylabel('\tau')
    grid on
    


figure('name','Desired Torque vs Angle')
    subplot(3,1,1)
    plot(rad2deg(q1(Middle:end)),torqueDesire(1,Middle:end),'linewidth',2)
    title('Initial Desired Torque-Angle Profile','FontWeight','normal','FontSize',16,'FontName','Times');
    xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('\tau_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
    subplot(3,1,2)
    plot(InRangeShifter(rad2deg(q2(Middle:end))),torqueDesire(2,Middle:end),'linewidth',2)
    xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('\tau_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
    subplot(3,1,3)
    plot(rad2deg(q3(Middle:end)),torqueDesire(3,Middle:end),'linewidth',2)
    xlabel('q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('\tau_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
    
Y=[q1(1:end)',q2(1:end)',q3(1:end)'];
% AnimBot3DOF(time(1:end),Y,L);

IntAbsUdqDesire_Opt=sum(sum(abs(torqueDesire(:,Middle:end).*Dq(:,Middle:end)),2))*Tres

%%  Generate Initial value for Optimization
% close all

tic
% DoF system
nn=3; % number of joints
% DoF of Optimization 
rQ=10; % Degree of joint trajectory
rU=4; % Degree of passive torque
% B matrix
B=eye(nn);
% WeightMatrix
Weight=[ 2 1 1]';
Sat=[1,1,1];

% Landa for [DQ  D2q ]
SeletcStr={'DQ','D2Q'};
SelectLanda=[0 1];
Landa=[1e-7 1e-7];


Time=time(Middle:end)-time(Middle);
Q1=q1(Middle:end);
Q2=q2(Middle:end);
Q3=q3(Middle:end);
% Q_CF1_U2=load('temp2','Q_Opt');
% Q1=Q_Opt(1,:);
% Q2=Q_Opt(2,:);
% Q3=Q_Opt(3,:);
XEF=xef(Middle:end);
YEF=yef(Middle:end);
DXEF=Dxef(Middle:end);
DYEF=Dyef(Middle:end);
D2XEF=D2xef(Middle:end);
D2YEF=D2yef(Middle:end);
TorqueDesQ1=torqueDesire(1,Middle:end);
TorqueDesQ2=torqueDesire(2,Middle:end);
TorqueDesQ3=torqueDesire(3,Middle:end);
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

    
CoefBLS_UPassive1 = LSParamPoly(Q1',TorqueDesQ1',rU,(Landa.*SelectLanda),Sat(1));    
CoefBLS_UPassive2 = LSParamPoly(Q2',TorqueDesQ2',rU,(Landa.*SelectLanda),Sat(2));    
CoefBLS_UPassive3 = LSParamPoly(Q3',TorqueDesQ3',rU,(Landa.*SelectLanda),Sat(3));    
CoefBLS_UPassive=[CoefBLS_UPassive1;CoefBLS_UPassive2;CoefBLS_UPassive3];


TorquePassiveQ1val=polyval(CoefBLS_UPassive1,Q1);
TorquePassiveQ2val=polyval(CoefBLS_UPassive2,Q2);
TorquePassiveQ3val=polyval(CoefBLS_UPassive3,Q3);
TorquePassiveVal=[TorquePassiveQ1val; TorquePassiveQ2val; TorquePassiveQ3val];
    

% Integral Matrix
Cost=0;
CostSub=0;
for i=1:nn
    QQ=[];
    DQ=[];
    for j=1:length(QVal(i,:))
         QQ(j,:) = QVal(i,j).^(rU:-1:0)';
         DQ(j,:) = ([QVal(i,j).^(rU-1:-1:0) 0].*(rU:-1:0))';
    end
    
    Cost=Cost + ...
          Weight(i)*( 1/2*(TorqueDesire(i,:)' - QQ*CoefBLS_UPassive((i-1)*(rU+1)+1:(i)*(rU+1)) )'* (TorqueDesire(i,:)' - QQ*CoefBLS_UPassive((i-1)*(rU+1)+1:(i)*(rU+1)) )+ ...
                      sum(Landa.*SelectLanda)* CoefBLS_UPassive((i-1)*(rU+1)+1:(i)*(rU+1))'*(DQ'*DQ)*CoefBLS_UPassive((i-1)*(rU+1)+1:(i)*(rU+1)) );
%     Cost=Cost + ...
%           Weight(i)*1/2*( (TorqueDesire(i,:)' - QQ*CoefBLSI )'*(TorqueDesire(i,:)' - QQ*CoefBLSI)/Sat(i)^2 + ...
%                 Landa(1)* CoefBLSI'*(DQ'*DQ)*CoefBLSI )*Tres;%+...
%           
    CostSub=CostSub + ...
          Weight(i)*( 1/2*(TorqueDesire(i,:)')'*(TorqueDesire(i,:)')+...
                     -1/2*(TorqueDesire(i,:)')'*QQ*((QQ'*QQ)+sum(Landa.*SelectLanda)*(DQ'*DQ))^-1 * QQ' * (TorqueDesire(i,:)'));
end


BetaOptimal=[CoefBLS_UPassive1;CoefBLS_UPassive2;CoefBLS_UPassive3];
TorquePassiveQ1valOptimal=polyval(BetaOptimal(0*(rU+1)+1:1*(rU+1)),Q1);
TorquePassiveQ2valOptimal=polyval(BetaOptimal(1*(rU+1)+1:2*(rU+1)),Q2val);
TorquePassiveQ3valOptimal=polyval(BetaOptimal(2*(rU+1)+1:3*(rU+1)),Q3val);
TorquePassiveValOptimal=[TorquePassiveQ1valOptimal; TorquePassiveQ2valOptimal; TorquePassiveQ3valOptimal];


% show time
% close all


Y=[Q1(1:end)',Q2(1:end)',Q3(1:end)'];
AnimBot3DOF(Time,Y,L);


figure('name','compare trajectory')
    subplot(3,1,1)
    plot(Time,Q1)
    hold all
    plot(Time,Q1val,'-.')
    hold off
    grid on
    subplot(3,1,2)
    plot(Time,Q2)
    hold all
    plot(Time,Q2val,'-.')
    hold off
    grid on
    subplot(3,1,3)
    plot(Time,Q3)
    hold all
    plot(Time,Q3val,'-.')
    hold off
    grid on

figure('name','Path (3DoF)')
    plot(xef,yef,'linewidth',2.5,'linestyle','-','color','g')
    hold on 
    plot(RPosVal(1,:),RPosVal(2,:),'linewidth',2,'linestyle','-.','color','b')
    legend('Desired','Static Path Planing')
    hold off
    axis equal

figure('name','Passive Torques')
    subplot(3,2,1)
    plot(Q1,TorqueDesQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q1,TorquePassiveQ1val,'linewidth',2,'linestyle','-.','color','r')
    plot(Q1val,TorquePassiveQ1valOptimal,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    xlabel('q_1')
    ylabel('\tau_1')
    legend('Desired Torque','PassiveTorque')
    subplot(3,2,2)
    plot(Q1,TorqueDesQ1-TorquePassiveQ1val,'linewidth',2)
    grid on
    xlabel('q_1')
    ylabel('\tau_1')
    title('Active Torque')
    
    subplot(3,2,3)
    plot(Q2,TorqueDesQ2,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q2,TorquePassiveQ2val,'linewidth',2,'linestyle','-.','color','r')
    plot(Q2val,TorquePassiveQ2valOptimal,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    xlabel('q_2')
    ylabel('\tau_2')
    legend('Desired Torque','PassiveTorque')
    subplot(3,2,4)
    plot(Q2,TorqueDesQ2-TorquePassiveQ2val,'linewidth',2)
    grid on
    xlabel('q_2')
    ylabel('\tau_2')
    title('Active Torque')
    
    subplot(3,2,5)
    plot(Q3,TorqueDesQ3,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q3,TorquePassiveQ3val,'linewidth',2,'linestyle','-.','color','r')
    plot(Q3val,TorquePassiveQ3valOptimal,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    xlabel('q_3')
    ylabel('\tau_3')
    legend('Desired Torque','PassiveTorque')
    subplot(3,2,6)
    plot(Q3,TorqueDesQ3-TorquePassiveQ3val,'linewidth',2)
    grid on
    xlabel('q_3')
    ylabel('\tau_3')
    title('Active Torque')
    
 toc
 
 
TorqueActive=TorqueDesire-TorquePassiveVal;
IntAbsUdqActive=sum(sum(abs(TorqueActive.*Dq(:,Middle:end)),2))*Tres
IntAbsUdqDesire=sum(sum(abs(TorqueDesire.*Dq(:,Middle:end)),2))*Tres

%% Optimization

Degree=[nn rQ rU];
Initial=[Alpha_Q1 Alpha_Q2 Alpha_Q3];
% % 
% Initail2=Initial;
% Initial=x;

% WeightMatrix
% Weight=[ 10 1 1]';
% Landa=.85;
% B matrix
% B=eye(nn);

% dynamic parameters
mL1=m;
mL2=m;
mL3=m;
LL1=L;
LL2=L;
LL3=L;

tic
MaxFunEvals_Data=3000*(rQ);
MaxIter_Data=1000;
TolFun_Data=1e-8;
TolX_Data=1e-8;
TolCon_Data=1e-6;
Algorithm='sqp';
Algorithm='interior-point';
Rand=5000*1e-6;
NewInit=Initial+Rand*(randn(1,3*(rQ+length(rQ))));

CostFun   = @(Alpha)CF3_TorqueCost(Alpha,Time,Degree,Tres,Weight,(Landa.*SelectLanda),Sat,g,mL1,mL2,mL3,LL1,LL2,LL3);
NonConstr = @(Alpha)CF3_NonLinearConstraint(Alpha,Time,Tres,Degree,L,XEF,YEF);



[x,fval,exitflag,output,lambda,grad,hessian] = ...
    Op_FmisCon_SQP(CostFun,NonConstr,NewInit,MaxFunEvals_Data,MaxIter_Data,TolFun_Data,TolX_Data,TolCon_Data,Algorithm);




%% ShowTime
[TorqueDesire_X0,TorqueActive_X0,Q_X0,D1Q_X0,D2Q_X0,BetaOptimal_X0,IntU2_X0,IntUdq_X0,IntAbsUdq_X0,IntAbsUdqDesire_X0,CostSlopeD1Q_X0,CostSlopeD2Q_X0,RMSError_X0]=...    
                        ShowTime(Initial,Time,Tres,Degree,Weight,(Landa.*SelectLanda),Sat,[],[] ,XEF,YEF,m,L,g,[],'Show','2Cycle','CostC','Initial');
[TorqueDesire_Opt,TorqueActive_Opt,Q_Opt,D1Q_Opt,D2Q_Opt,BetaOptimal_Opt,IntU2_Opt,IntUdq_Opt,IntAbsUdq_Opt,IntAbsUdqDesire_Opt,CostSlopeD1Q_Opt,CostSlopeD2Q_Opt,RMSError_Opt]=...
                        ShowTime(x,Time,Tres,Degree,Weight,(Landa.*SelectLanda),Sat,[],[],XEF,YEF,m,L,g,[],'Show','2Cycle','CostC','Optimized');
                    
TotalCost_X0  = IntU2_X0    + sum(Landa.*SelectLanda.*[CostSlopeD1Q_X0 CostSlopeD2Q_X0]);
TotalCost_Opt = IntU2_Opt   + sum(Landa.*SelectLanda.*[CostSlopeD1Q_Opt CostSlopeD2Q_Opt]);
IntAbsUdqDesire_X0=sum(sum(abs(TorqueDesire_X0.*D1Q_X0),2))*Tres;
IntAbsUdqDesire_Opt=sum(sum(abs(TorqueDesire_Opt.*D1Q_Opt),2))*Tres;
IntAbsUdqActive_X0=sum(sum(abs(TorqueActive_X0.*D1Q_X0),2))*Tres;
IntAbsUdqActive_Opt=sum(sum(abs(TorqueActive_Opt.*D1Q_Opt),2))*Tres;

DegreeStr=sprintf('\n   rQ:%3d\n   rU:%3d\n   Cost Type:  %s\n   Requlation Type:  %s\n ',rQ,rU, 'Cast 2b',SeletcStr{find(SelectLanda)});
Title=sprintf('%22s  % 11s %11s % 8s % 12s % 15s % 15s'  ,   'IntU2','C.S.(D1)','C.S.(D2)','Total','RMS Err','Req. Work','Actuator Work');
Result_X0 =sprintf('%-11s %11.2e %11.2e %11.2e % 10.2e % 10.2e % 15.2e % 14.2e  ',   'Initial:',  IntU2_X0, CostSlopeD1Q_X0, CostSlopeD2Q_X0, TotalCost_X0, RMSError_X0, sum(IntAbsUdqDesire_X0), IntAbsUdqActive_X0);
Result_Opt=sprintf('%-11s %11.2e %11.2e %11.2e % 10.2e % 10.2e % 15.2e % 14.2e\n',   'Optimized:',IntU2_Opt,CostSlopeD1Q_Opt,CostSlopeD2Q_Opt,TotalCost_Opt,RMSError_Opt,sum(IntAbsUdqDesire_Opt),IntAbsUdqActive_Opt);
display(output.message)
disp(DegreeStr)
disp(Title)
disp(Result_X0)
disp(Result_Opt)

TorqueActive=TorqueDesire-TorquePassiveVal;
IntAbsUdqActive=sum(sum(abs(TorqueActive.*Dq(:,Middle:end)),2))*Tres;
IntAbsUdqDesire=sum(sum(abs(TorqueDesire.*Dq(:,Middle:end)),2))*Tres;
IntAbsUdqActive=sum(sum(abs(TorqueActive.*Dq(:,Middle:end)),2))*Tres;
IntAbsUdqDesire=sum(sum(abs(TorqueDesire.*Dq(:,Middle:end)),2))*Tres;
disp([IntAbsUdqActive IntAbsUdqDesire ])

% toc


%% Scale and shift profile

ThetaShiftScale=[];
ThetaStepscale=[];
tauShiftScale=[];
tauShiftMain=[];

for Joint=1:nn


    ThetaStep=( (max(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)))  - min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)))) / (floor(size(Q_Opt,2)/2)));
    ThetaS=min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2))) :ThetaStep:max(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)));
    ThetaShift=ThetaS-min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)));
    ThetaShiftScale{Joint} = ThetaShift* floor(deg2rad(270) /  max(ThetaShift));
    ThetaStepscale{Joint}  = ThetaStep* floor(deg2rad(270) /  max(ThetaShift));
    ThetaShiftMain{Joint} = ThetaShift;
    ThetaStepMain{Joint}  = ThetaStep;
%     ThetaShiftScale{Joint} = ThetaShiftMain{Joint};
%     



    tau=polyval(BetaOptimal_Opt((rU+1)*(Joint-1)+1: (Joint)*(rU+1)),ThetaS);    
    TauMin=abs(min(tau));
    tauShift=tau+TauMin*3;
    tauShiftMain{Joint}=(tauShift);
    tauShiftScale{Joint}=(tauShift)/max(tauShift);
    


    figure('name',['Joint ', num2str( Joint)])
        subplot(2,1,1)
        plot(ThetaShiftScale{Joint},tauShiftScale{Joint})
        xlabel('ThetaShiftScale')
        ylabel('tauShiftScale')
        grid on

        subplot(2,1,2)
        plot(ThetaShiftScale{Joint}(1:end-1),diff(tauShiftScale{Joint})./diff(ThetaShiftScale{Joint}))
        xlabel('ThetaShiftScale')
        ylabel('DtauShiftScale')
        grid on
        hold on
        plot(ThetaShiftScale{Joint},+ones(size(ThetaShiftScale{Joint})),'r-.')
        plot(ThetaShiftScale{Joint},-ones(size(ThetaShiftScale{Joint})),'r-.')
        hold off
end
%% Ga for find best Param

k0=2000;
R0=1;
q00=1;

nvars=3;
lb=[50 5e-2 5e-2];
PopInitRange=[lb; k0 R0 q00];
PopulationSize=200;
InitialPopulation=[k0*rand(PopulationSize,1) R0*rand(PopulationSize,1) q00*rand(PopulationSize,1)];

ParamA=[];
for Joint=1:nn
    disp(['Joint:',num2str(Joint)])
    CostParam=@(Param)GA_CostParamNonLinearSpring(Param,ThetaStepscale{Joint},ThetaShiftScale{Joint},tauShiftScale{Joint});
%     CostParam=@(Param)GA_CostParamNonLinearSpring(Param,ThetaStepMain{Joint},ThetaShiftMain{Joint},tauShiftScale{Joint});

    [ParamA{Joint},fval,exitflag,output,population,score] = ...
        GA_FminCost(CostParam,nvars,lb,PopInitRange,PopulationSize,InitialPopulation);
    disp(output.message)
end
%%
% ParamA=[kg+600,Rg-.85,qg0+.25]
for Joint=1:nn
    k=ParamA{Joint}(1)*max(tauShiftMain{Joint});
    R=ParamA{Joint}(2)*2;
    q0=ParamA{Joint}(3);

    NonLinearSpring(ThetaStepscale{Joint},ThetaShiftScale{Joint},tauShiftScale{Joint}*max(tauShiftMain{Joint}),k,R,q0,.1,['Joint ',num2str(Joint)])
    % NonLinearSpring(ThetaStep,ThetaS,tauShift,k,R,q0,.1)
end
%% Control

OdeOpt= odeset('RelTol',1e-15,'AbsTol',1e-15*ones(1,7),'maxstep',1e-4);
InitState=([ Q_Opt(1:3,1)' , D1Q_Opt(1:3,1)', 0]);
Pos  =[[XEF   ,XEF(2:end)]   ;[YEF   ,YEF(2:end)]   ];
DPos =[[DXEF  ,DXEF(2:end)]  ;[DYEF  ,DYEF(2:end)]  ];
D2Pos=[[D2XEF ,D2XEF(2:end)] ;[D2YEF ,D2YEF(2:end)] ];
Time2=[Time , Time(2:end)+Time(end)];
% TorqueFF=[TorqueDesire_Opt,TorqueDesire_Opt(:,2:end)];
TorqueFF=[TorqueActive_Opt,TorqueActive_Opt(:,2:end)];
Qref=[Q_Opt, Q_Opt(:,2:end)];
DQref=[D1Q_Opt, D1Q_Opt(:,2:end)];
Mid=ceil(length(Time2)/2);

Kp=0*10000;
Kd=0*200;

[T,Q_ode] =...
    ode15s(@(t,Y)SirDyn(t,Y,g,L,m,Pos,DPos,D2Pos,Time2,Kp,Kd,TorqueFF,BetaOptimal_Opt*1,rU),...
            Time,InitState,OdeOpt);

QO1=Q_ode(:,1)';
QO2=Q_ode(:,2)';
QO3=Q_ode(:,3)';
RPosODE=L*[cos(QO1)+cos(QO1+QO2)+cos(QO1+QO2+QO3);
        sin(QO1)+sin(QO1+QO2)+sin(QO1+QO2+QO3)];
%%
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
%%
plot(Q_ode)

[TorqueActive,TorquePassive,TorqueControl]=...
    TorqueCalculatorControl(T,Q_ode,g,m,m,m,L,L,L,Pos,DPos,D2Pos,Time2,Kp,Kd,TorqueFF,BetaOptimal_Opt,rU);
figure('name','Torque')
subplot(3,1,1)
plot(Time2(Mid:end),TorqueActive(:,Mid:end))
subplot(3,1,2)
plot(Time2(Mid:end),TorquePassive(:,Mid:end))
subplot(3,1,3)
plot(Time2(Mid:end),TorqueControl(:,Mid:end))
legend('\tau_1','\tau_2','\tau_3')
xlabel('time (s)')
ylabel('\tau')

EnergyABS=Q_ode(end,end)
EnergyABS=sum(sum(abs(Torque(Mid:end,:).*Q_ode(Mid:end,4:6)),2))*Tres
