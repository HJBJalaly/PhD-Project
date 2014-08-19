
%% Initialization

clear
close all
home

% envirment
g=9.81*0;

% parameters of robot
m=1;
L=1;

KKmat=[100;1000;500];

% Desired motion
f=1;
A=1;
Tres=0.01;
time=0:Tres:4/f;

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

%% create a active motion with 3 active joints as initial input of optimization
home

% EF motion

f=1;
A=1;
Tres=0.005;
time=0:Tres:2*f;

% % Circle motion
% Xef=A*cos(2*pi/f*time)+1.5;
% Yef=A*sin(2*pi/f*time)+1;
% DXef=-(2*pi/f)*A*sin(2*pi/f*time);
% DYef= (2*pi/f)*A*cos(2*pi/f*time);
% D2Xef=-(2*pi/f)^2*A*cos(2*pi/f*time);
% D2Yef=-(2*pi/f)^2*A*sin(2*pi/f*time);
% q1=deg2rad(0);
% q2=deg2rad(8.031);
% q3=deg2rad(51.317);

% Line motion
Xef=A*cos(2*pi/(1*f)*time);
Yef=2.5*ones(size(time));
DXef=-(2*pi/f)*A*sin(2*pi/(1*f)*time);
DYef= 2*zeros(size(time));
D2Xef=-(2*pi/f)^2*A*cos(2*pi/(1*f)*time);
D2Yef= 2*zeros(size(time));
q1=deg2rad( 41.7023); % 41.7023 deg for y=2.5m  ,  22.2 deg for y=2.0m
q2=deg2rad(18.2663); % 18.2663 deg for y=2.5m  ,  29.468 deg for y=2.0m
q3=deg2rad(44.3377); % 44.3377 deg for y=2.5m  ,  71.431 deg for y=2.0m


Pos=[Xef;Yef];
DPos=[DXef;DYef];
D2Pos=[D2Xef;D2Yef];


% Joint Trajectories

for i=1:length(time)-1
    lPos=L*[cos(q1(i))+cos(q1(i)+q2(i))+cos(q1(i)+q2(i)+q3(i));
            sin(q1(i))+sin(q1(i)+q2(i))+sin(q1(i)+q2(i)+q3(i))];

    cPos=Pos(:,i);
    dx=cPos-lPos;

    JJ=L*[-sin(q1(i))-sin(q1(i)+q2(i))-sin(q1(i)+q2(i)+q3(i)), -sin(q1(i)+q2(i))-sin(q1(i)+q2(i)+q3(i)), -sin(q1(i)+q2(i)+q3(i));
           cos(q1(i))+cos(q1(i)+q2(i))+cos(q1(i)+q2(i)+q3(i)),  cos(q1(i)+q2(i))+cos(q1(i)+q2(i)+q3(i)),  cos(q1(i)+q2(i)+q3(i))]; 

    dq=JJ'*(JJ*JJ')^-1*dx;

    q1(i+1)=q1(i)+dq(1);
    q2(i+1)=q2(i)+dq(2);
    q3(i+1)=q3(i)+dq(3);
    
end

RPos=L*[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
        sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];

    
% Required torque for desired path
q=[q1;q2;q3];
Dq=differential(q,time,Tres);
D2q=differential(Dq,time,Tres);
Torque=TorqueCalculator(D2q,Dq,q,g,m,m,m,L,L,L);

Middle=ceil(length(time)/2);

% show    
close all
figure('name','Path (3DoF)')
    plot(Xef,Yef,'linewidth',2.5,'linestyle','-','color','g')
    hold on 
    plot(RPos(1,Middle:end),RPos(2,Middle:end),'linewidth',2,'linestyle','-.','color','b')
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
    
figure('name','Torque vs Time')
    plot(time(Middle:end),Torque(:,Middle:end))
    legend('\tau_1','\tau_2','\tau_3')
    xlabel('time (s)')
    ylabel('\tau')
    grid on
    


figure('name','Torque vs Angle')
    subplot(3,1,1)
    plot(q1(Middle:end),Torque(1,Middle:end))
    grid on
    subplot(3,1,2)
    plot(q2(Middle:end),Torque(2,Middle:end))
    grid on
    subplot(3,1,3)
    plot(q3(Middle:end),Torque(3,Middle:end))
    grid on
    
Y=[q1(1:end)',q2(1:end)',q3(1:end)'];
AnimBot3DOF(time(1:end),Y,L);

%%  Generate Initial value for Optimization

Degree=6;

Time=time(Middle:end)-1;
Q1=q1(Middle:end);
Q2=q2(Middle:end);
Q3=q3(Middle:end);
XEF=Xef(Middle:end);
YEF=Yef(Middle:end);

[CoefP_q1,CoefB_q1]= BezierCoeffinet(Time,Q1,Degree);
[CoefP_q2,CoefB_q2]= BezierCoeffinet(Time,Q2,Degree);
[CoefP_q3,CoefB_q3]= BezierCoeffinet(Time,Q3,Degree);

Q1val=polyval(CoefP_q1,Time);
Q2val=polyval(CoefP_q2,Time);
Q3val=polyval(CoefP_q3,Time);

RPos=L*[cos(Q1val)+cos(Q1val+Q2val)+cos(Q1val+Q2val+Q3val);
        sin(Q1val)+sin(Q1val+Q2val)+sin(Q1val+Q2val+Q3val)];

% show    
close all
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

plot(Xef,Yef,'linewidth',2.5,'linestyle','-','color','g')
hold on 
plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','-.','color','b')
legend('Desired','Static Path Planing')
hold off
axis equal

Initial=[CoefP_q1 CoefP_q2 CoefP_q3];
mL1=m;
mL2=m;
mL3=m;
LL1=L;
LL2=L;
LL3=L;

%% Optimization

tic
MaxFunEvals_Data=1000*Degree;
MaxIter_Data=1000;
TolFun_Data=1e-8;
TolX_Data=1e-8;
TolCon_Data=1e-8;
Algorithm='sqp';
% Algorithm='interior-point';

Select=[ 0 0 1]; SeletcStr={'IntU2','IntAbsUdq','IntUdq'};
Weight=[ 1 1 1]';

CostFun=@(Coef)TorqueCost(Coef,Time,Degree,Tres,Select,Weight,g,mL1,mL2,mL3,LL1,LL2,LL3);
NonCons=@(Coef)NonLinearConstraint(Coef,Time,Tres,Degree,L,XEF,YEF);

[x,fval,exitflag,output,lambda,grad,hessian] = ...
    FmisCon_SQP(CostFun,NonCons,Initial,MaxFunEvals_Data,MaxIter_Data,TolFun_Data,TolX_Data,TolCon_Data,Algorithm);


[Torque_X0,Q_X0,D1Q_X0,D2Q_X0,IntU2_X0,IntUdq_X0,IntAbsUdq_X0,RMSError_X0]=...
                        ShowTime(Initial,Time,Tres,Degree,Weight,XEF,YEF,m,L,g,0);
[Torque_Opt,Q_Opt,D1Q_Opt,D2Q_Opt,IntU2_Opt,IntUdq_Opt,IntAbsUdq_Opt,RMSError_Opt]=...
                        ShowTime(x,Time,Tres,Degree,Weight,XEF,YEF,m,L,g,1);

DegreeStr=sprintf('\n   Degree:     %d\n   Cost Type:  %s\n',Degree,SeletcStr{find(Select)});
Title=sprintf('% 24s % 12s % 12s % 12s','IntU2','IntAbsUdq','IntUdq','RMS Err');
Result_X0 =sprintf('Initial:    % 12.2e % 12.2e % 12.2e % 12.2e',IntU2_X0,IntAbsUdq_X0,IntUdq_X0,RMSError_X0);
Result_Opt=sprintf('Optimized:  % 12.2e % 12.2e % 12.2e % 12.2e\n',IntU2_Opt,IntAbsUdq_Opt,IntUdq_Opt,RMSError_Opt);
display(output.message)
disp(DegreeStr)
disp(Title)
disp(Result_X0)
disp(Result_Opt)
toc


%% Make Nonlinear Spring

Joint=2;

ThetaStep=( (max(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)))  - min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2))))/200);
ThetaS=min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2))) :ThetaStep:max(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)));
ThetaShift=ThetaS-min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)));
ThetaShiftScale = ThetaShift* floor(deg2rad(270) /  max(ThetaShift));
tau=interp1(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)),Torque_Opt(Joint, 1: floor(size(Q_Opt,2)/2)),ThetaS);
TauMin=abs(min(tau));
tauShift=tau+TauMin*1.5;
tauShiftScale=(tauShift)/10;




%%

k0=1000;
R0=1;
q00=1;

nvars=3;
lb=[50 5e-2 5e-2];
PopInitRange=[lb; k0 R0 q00];
PopulationSize=1000;
InitialPopulation=[k0*rand(PopulationSize,1) R0*rand(PopulationSize,1) q00*rand(PopulationSize,1)];
CostParam=@(Param)FindBestParamCost(Param,ThetaStep,ThetaShiftScale,tauShiftScale);

[ParamA,fval,exitflag,output,population,score] = ...
    Ga_FindParamOfNonLinearSpring(CostParam,nvars,lb,PopInitRange,PopulationSize,InitialPopulation);
disp(output.message)

%%
k=ParamA(1);
R=ParamA(2);
q0=ParamA(3);

NonLinearSpring(ThetaStep,ThetaShiftScale,tauShiftScale,k,R,q0,.5)


