function CF2_Main()



%% Initialization

clear
close all
home


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
clear
clc
% envirment
g=9.81*0;
% parameters of robot
m=1;
L=1;
% EF motion
f=1;
A=1;
phi=pi;
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
xef=A*cos(2*pi/(1*f)*time+phi);
yef=2.5*ones(size(time));
Dxef=-(2*pi/f)*A*sin(2*pi/(1*f)*time+phi);
Dyef= 2*zeros(size(time));
D2xef=-(2*pi/f)^2*A*cos(2*pi/(1*f)*time+phi);
D2yef= 2*zeros(size(time));
q1=deg2rad( 41.7023); % 41.7023 deg for y=2.5m  ,  22.2 deg for y=2.0m
q2=deg2rad(18.2663); % 18.2663 deg for y=2.5m  ,  29.468 deg for y=2.0m
q3=deg2rad(44.3377); % 44.3377 deg for y=2.5m  ,  71.431 deg for y=2.0m


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
    subplot(3,1,2)
    plot(time(Middle:end),q2(Middle:end))
    grid on
    subplot(3,1,3)
    plot(time(Middle:end),q3(Middle:end))
    grid on
    
figure('name','Desired Torque vs Time')
    plot(time(Middle:end),torqueDesire(:,Middle:end))
    legend('\tau_1','\tau_2','\tau_3')
    xlabel('time (s)')
    ylabel('\tau')
    grid on
    


figure('name','Desired Torque vs Angle')
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

%%  Generate Initial value for Optimization
tic
% DoF system
nn=3; % number of joints
% DoF of Optimization 
rQ=7; % Degree of joint trajectory
rU=7; % Degree of passive torque
% B matrix
B=eye(nn);
% WeightMatrix
Weight=[ 1 1 1]';
Landa=0.98;


Time=time(Middle:end)-time(end)/2;
Q1=q1(Middle:end);
Q2=q2(Middle:end);
Q3=q3(Middle:end);
XEF=xef(Middle:end);
YEF=yef(Middle:end);
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

[Beta_UPassive1,BezireCoef_TorP1]= BezierCoeffinet(Q1,TorqueDesQ1,rU);
[Beta_UPassive2,BezireCoef_TorP2]= BezierCoeffinet(Q2,TorqueDesQ2,rU);
[Beta_UPassive3,BezireCoef_TorP3]= BezierCoeffinet(Q3,TorqueDesQ3,rU);


TorquePassiveQ1val=polyval(Beta_UPassive1,Q1);
TorquePassiveQ2val=polyval(Beta_UPassive2,Q2);
TorquePassiveQ3val=polyval(Beta_UPassive3,Q3);
TorquePassiveVal=[TorquePassiveQ1val; TorquePassiveQ2val; TorquePassiveQ3val];
    
% Omega matrix
Omega = B * (B'*B)^-1 *diag(Weight) * (B'*B)^-1 * B';

% Integral Matrix
Iu=0;
Iq=zeros(nn*(rU+1),nn*(rU+1));
Iqu=zeros(1,(rU+1)*nn);
for tt=1:length(Time)
   
    QQ_conc=zeros((rU+1)*nn,nn);% \underline{\underline{\mathcal{Q}}}^{r_u}
    for joint=1:nn
        QQ_rU_Joint=QVal(joint,tt).^(rU:-1:0)';
        QQ_conc(1+(joint-1)*(rU+1):(joint)*(rU+1),joint)= QQ_rU_Joint;
    end

    Iu = Iu + TorqueDesire(:,tt)'*Omega*TorqueDesire(:,tt);

    Iq = Iq + QQ_conc*Omega*QQ_conc';

    Iqu= Iqu+ TorqueDesire(:,tt)'*Omega*QQ_conc';
   
end

Iu=Iu*Tres;
Iq=Iq*Tres;
Iqu=Iqu*Tres;

syms s real;


DqQ_rU_Joint=[(rU:-1:1).*(s.^(rU-1:-1:0)')', 0]';
QQ=double(int(DqQ_rU_Joint*DqQ_rU_Joint',0,3/2*pi));

Psi=zeros(rU+1,rU+1,nn);
Idq=zeros(rU+1,rU+1,nn);
Idq_conc=zeros((rU+1)*nn,(rU+1)*nn);

for Joint=1:nn
    a_hat(Joint) = max(TorquePassiveVal(Joint,:)) - min(TorquePassiveVal(Joint,:));
    b_hat(Joint) = min(TorquePassiveVal(Joint,:));
    c_hat(Joint) =(max(QVal(Joint,:)) - min(QVal(Joint,:)))* 2 / 3/pi;
    d_hat(Joint) = min(QVal(Joint,:));

    for kk=1:rU+1
        Psi(kk,:,Joint)  = [zeros(1,kk-1),sym2poly( (c_hat(Joint)*s+ d_hat(Joint))^(rU+1-(kk))  ) ];
%         Psi2(kk,:,Joint) = [zeros(1,kk-1), (c_hat(Joint)^(rU+1-(kk)))* poly(-d_hat(Joint)/c_hat(Joint)*ones(1,rU+1-(kk))  ) ];
    end

    QVal_hat(Joint,:)= ( QVal(Joint,:)-d_hat(Joint))/c_hat(Joint);

    Idq(:,:,Joint) = Weight(Joint)*c_hat(Joint)* Psi(:,:,Joint)*QQ*Psi(:,:,Joint)';
    
    Idq_conc((Joint-1)*(rU+1)+1:(Joint)*(rU+1),(Joint-1)*(rU+1)+1:(Joint)*(rU+1)) = Idq(:,:,Joint);

end


% BetaOptimal=Landa*(Landa*Iq+(1-Landa)*Idq_conc)^-1*Iqu';
SVDsol=SVDBlockInvertor((Landa*Iq+(1-Landa)*Idq_conc),nn,rU+1,0.0001);
BetaOptimal=Landa*SVDsol*Iqu';
TorquePassiveQ1valOptimal=polyval(BetaOptimal(0*(rU+1)+1:1*(rU+1)),Q1);
TorquePassiveQ2valOptimal=polyval(BetaOptimal(1*(rU+1)+1:2*(rU+1)),Q2);
TorquePassiveQ3valOptimal=polyval(BetaOptimal(2*(rU+1)+1:3*(rU+1)),Q3);
TorquePassiveValOptimal=[TorquePassiveQ1valOptimal; TorquePassiveQ2valOptimal; TorquePassiveQ3valOptimal];


% show time
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
    plot(xef,yef,'linewidth',2.5,'linestyle','-','color','g')
    hold on 
    plot(RPosVal(1,:),RPosVal(2,:),'linewidth',2,'linestyle','-.','color','b')
    legend('Desired','Static Path Planing')
    hold off
    axis equal

figure('name','Passive Torques')
    subplot(3,1,1)
    plot(Q1,TorqueDesQ1,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q1,TorquePassiveQ1val,'linewidth',2,'linestyle','-.','color','r')
    plot(Q1,TorquePassiveQ1valOptimal,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    xlabel('q_1')
    ylabel('\tau_1')
    legend('Desired Torque','PassiveTorque')
    
    subplot(3,1,2)
    plot(Q2,TorqueDesQ2,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q2,TorquePassiveQ2val,'linewidth',2,'linestyle','-.','color','r')
    plot(Q2,TorquePassiveQ2valOptimal,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    xlabel('q_2')
    ylabel('\tau_2')
    legend('Desired Torque','PassiveTorque')
    
    subplot(3,1,3)
    plot(Q3,TorqueDesQ3,'linewidth',2,'linestyle','-','color','b')
    hold on
    plot(Q3,TorquePassiveQ3val,'linewidth',2,'linestyle','-.','color','r')
    plot(Q3,TorquePassiveQ3valOptimal,'linewidth',2,'linestyle','-.','color','g')
    hold off
    grid on
    xlabel('q_3')
    ylabel('\tau_3')
    legend('Desired Torque','PassiveTorque')
 toc
%% Optimization


Degree=[nn rQ rU];
Initial=[Alpha_Q1 Alpha_Q2 Alpha_Q3];
% % 
%  Xt=x;
% Xt2=x;
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
TolFun_Data=1e-5;
TolX_Data=1e-5;
TolCon_Data=1e-5;
% Algorithm='sqp';
Algorithm='interior-point';
Rand=5000*1e-10;
MinSinValue=0.005;

CostFun   = @(Alpha)CF2_TorqueCost(Alpha,Time,Degree,Tres,Weight,Landa,QQ,B,g,mL1,mL2,mL3,LL1,LL2,LL3,MinSinValue);
NonConstr = @(Alpha)CF2_NonLinearConstraint(Alpha,Time,Tres,Degree,L,XEF,YEF);


[x,fval,exitflag,output,lambda,grad,hessian] = ...
    Op_FmisCon_SQP(CostFun,NonConstr,Initial+Rand*(randn(1,3*(sum(rQ)+length(rQ)))),MaxFunEvals_Data,MaxIter_Data,TolFun_Data,TolX_Data,TolCon_Data,Algorithm);

%%
[Torque_X0,Q_X0,D1Q_X0,D2Q_X0,BetaOptimal_X0,Nothing,IntU2_X0,IntUdq_X0,IntAbsUdq_X0,IntAbsUdqDesire_X0,CostSlope_X0,Nothing,Nothing,RMSError_X0]=...
                        ShowTime(Initial,Time,Tres,Degree,Weight,Landa,[],QQ,B,XEF,YEF,m,L,g,MinSinValue,'Show','2Cycle','CostB','Initial');
[Torque_Opt,Q_Opt,D1Q_Opt,D2Q_Opt,BetaOptimal_Opt,Nothing,IntU2_Opt,IntUdq_Opt,IntAbsUdq_Opt,IntAbsUdqDesire_Opt,CostSlope_Opt,Nothing,Nothing,RMSError_Opt]=...
                        ShowTime(x      ,Time,Tres,Degree,Weight,Landa,[],QQ,B,XEF,YEF,m,L,g,MinSinValue,'Show','2Cycle','CostB','Optimized');

TotalCost_X0  = Landa* IntU2_X0    + (1-Landa)*CostSlope_X0;
TotalCost_Opt = Landa* IntU2_Opt   + (1-Landa)*CostSlope_Opt;
             

DegreeStr=sprintf('\n   rQ:%3d\n   rU:%3d\n   Cost Type:  %s\n',rQ,rU, 'Cast 2a');
Title=sprintf('%22s %10s % 11s % 11s % 15s % 15s','IntU2','C.S.','Total','RMS Err','Req. Work','Actuator Work');
Result_X0 =sprintf('%-11s %11.2e %11.2e % 10.2e % 10.2e % 15.2e % 14.2e  ','Initial:',IntU2_X0,CostSlope_X0,TotalCost_X0,RMSError_X0,sum(IntAbsUdqDesire_X0),IntAbsUdq_X0);
Result_Opt=sprintf('%-11s %11.2e %11.2e % 10.2e % 10.2e % 15.2e % 14.2e\n','Optimized:',IntU2_Opt,CostSlope_Opt,TotalCost_Opt,RMSError_Opt,sum(IntAbsUdqDesire_Opt),IntAbsUdq_Opt);
display(output.message)
disp(DegreeStr)
disp(Title)
disp(Result_X0)
disp(Result_Opt)
toc


%% Scale and shift profile

ThetaShiftScale=[];
ThetaStepscale=[];
tauShiftScale=[];

for Joint=1:nn


    ThetaStep=( (max(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)))  - min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)))) / (floor(size(Q_Opt,2)/2)));
    ThetaS=min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2))) :ThetaStep:max(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)));
    ThetaShift=ThetaS-min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)));
    ThetaShiftScale{Joint} = ThetaShift* floor(deg2rad(270) /  max(ThetaShift));
    ThetaStepscale{Joint}  = ThetaStep* floor(deg2rad(270) /  max(ThetaShift));


    tau=polyval(BetaOptimal_Opt((rU+1)*(Joint-1)+1: (Joint)*(rU+1)),ThetaS);    
    TauMin=abs(min(tau));
    tauShift=tau+TauMin*2;
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

k0=1000;
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

    [ParamA{Joint},fval,exitflag,output,population,score] = ...
        GA_FminCost(CostParam,nvars,lb,PopInitRange,PopulationSize,InitialPopulation);
    disp(output.message)
end
%%
% ParamA=[kg+600,Rg-.85,qg0+.25]
for Joint=1:nn
    k=ParamA{Joint}(1);
    R=ParamA{Joint}(2);
    q0=ParamA{Joint}(3);

    NonLinearSpring(ThetaStepscale{Joint},ThetaShiftScale{Joint},tauShiftScale{Joint},k,R,q0,.1,['Joint ',num2str(Joint)])
    % NonLinearSpring(ThetaStep,ThetaS,tauShift,k,R,q0,.1)
end
%% 