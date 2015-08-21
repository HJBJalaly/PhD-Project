
%% 2DoF

home
clear
clc
% envirment
g=9.81*1;
% parameters of robot
m=1.5;
L=1.5;
% EF motion
f=1;
A=.75;
phi=pi/2;
Tres=0.005;
time=0:Tres:4/f;


% Circle motion
f=0.5;
phi=pi/2;
Tres=0.005;
time=0:Tres:10/f;
Middle=ceil(9*length(time)/10);
xef=A*cos(2*pi*f*time)+0;
yef=A*sin(2*pi*f*time)+2;
Dxef=-(2*pi*f)*A*sin(2*pi*f*time);
Dyef= (2*pi*f)*A*cos(2*pi*f*time);
D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time);
D2yef=-(2*pi*f)^2*A*sin(2*pi*f*time);
q1=deg2rad(-60);
q2=deg2rad(48.031);
q3=deg2rad(-51.317);


% % Line motion: Horizontal
% phi=0;
% time=0:Tres:20/f;
% Middle=ceil(19*length(time)/20);
% xef=A*cos(2*pi*f*time+phi);
% yef=2.5*ones(size(time));
% Dxef=-(2*pi*f)*A*sin(2*pi*f*time+phi);
% Dyef= 2*zeros(size(time));
% D2xef=-(2*pi*f)^2*A*cos(2*pi*f*time+phi);
% D2yef= 2*zeros(size(time));
% q1=deg2rad( 42.7023); % 41.7023 deg for y=2.5m  ,  22.2 deg for y=2.0m
% q2=deg2rad(19.2663); % 18.2663 deg for y=2.5m  ,  29.468 deg for y=2.0m
% q3=deg2rad(80.3377); % 44.3377 deg for y=2.5m  ,  71.431 deg for y=2.0m



figure('name','Path (3DoF)')
plot(xef,yef,'linewidth',2.5,'linestyle','-','color','g')

Pos=[xef;yef];
DPos=[Dxef;Dyef];
D2Pos=[D2xef;D2yef];


for i=1:length(time)-1
    lPos=L*[cos(q1(i))+cos(q1(i)+q2(i));
            sin(q1(i))+sin(q1(i)+q2(i))];

    cPos=Pos(:,i+1);
    dx=cPos-lPos;

    JJ=L*[-sin(q1(i))-sin(q1(i)+q2(i)), -sin(q1(i)+q2(i));
           cos(q1(i))+cos(q1(i)+q2(i)),  cos(q1(i)+q2(i))]; 

    dq=JJ^-1*dx;
    q1(i+1)=q1(i)+dq(1);
    q2(i+1)=q2(i)+dq(2);
    
end



RPos=L*[cos(q1)+cos(q1+q2);
        sin(q1)+sin(q1+q2)];
hold all
plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','-.','color','b')
hold off
axis equal

q=[q1;q2];
Dq=differential(q,time,Tres);
D2q=differential(Dq,time,Tres);
torqueDesire=TorqueCalculator2(D2q,Dq,q,g,m,m,L,L);


figure('name','Path (3DoF)')
    plot(xef(Middle:end),yef(Middle:end),'linewidth',2.5,'linestyle','-','color','g')
    hold on 
    plot(RPos(1,Middle:end),RPos(2,Middle:end),'linewidth',2,'linestyle','-.','color','b')
    legend('Desired','Static Path Planing')
    hold off
    axis equal

figure('name','Joint Trajectories')
    subplot(2,1,1)
    plot(time(Middle:end),q1(Middle:end))
    grid on
    subplot(2,1,2)
    plot(time(Middle:end),q2(Middle:end))
    grid on
    
    
figure('name','Desired Torque vs Time')
    plot(time(Middle:end),torqueDesire(:,Middle:end))
    legend('\tau_1','\tau_2')
    xlabel('time (s)')
    ylabel('\tau')
    grid on
    

%%
% DoF system
nn=2; % number of joints
% DoF of Optimization 
rQ=8; % Degree of joint trajectory
rU=3; % Degree of passive torque
% B matrix
B=eye(nn);
% WeightMatrix
Weight=[ 3 2 2]';
Sat=[1,1,1];

% Landa for [DQ  D2q ]
SeletcStr={'DQ','D2Q'};
SelectLanda=[0 1];
Landa=[1e-7 1e-3];


Time=time(Middle:end)-time(Middle);
Q1=q1(Middle:end);
Q2=q2(Middle:end);
XEF=xef(Middle:end);
YEF=yef(Middle:end);
DXEF=Dxef(Middle:end);
DYEF=Dyef(Middle:end);
D2XEF=D2xef(Middle:end);
D2YEF=D2yef(Middle:end);
TorqueDesQ1=torqueDesire(1,Middle:end);
TorqueDesQ2=torqueDesire(2,Middle:end);
TorqueDesire=[TorqueDesQ1;TorqueDesQ2];


[Alpha_Q1,BezireCoef_q1]= BezierCoeffinet(Time,Q1,rQ);
[Alpha_Q2,BezireCoef_q2]= BezierCoeffinet(Time,Q2,rQ);
Q1val=polyval(Alpha_Q1,Time);
Q2val=polyval(Alpha_Q2,Time);
QVal=[Q1val;Q2val];

    
CoefBLS_UPassive1 = LSParamPoly(Q1',TorqueDesQ1',rU,(Landa.*SelectLanda),Sat(1));    
CoefBLS_UPassive2 = LSParamPoly(Q2',TorqueDesQ2',rU,(Landa.*SelectLanda),Sat(2));    
CoefBLS_UPassive=[CoefBLS_UPassive1;CoefBLS_UPassive2];

TorquePassiveQ1val=polyval(CoefBLS_UPassive1,Q1);
TorquePassiveQ2val=polyval(CoefBLS_UPassive2,Q2);
TorquePassiveVal=[TorquePassiveQ1val; TorquePassiveQ2val];
    

% Integral Matrix
IntU2=0;
Cost=0;
CostSub=0;
CostSlopeD1Q=0;
CostSlopeD2Q=0;
for i=1:nn
    QQ=[];
    DQ=[];
    CoefBLSI = LSParamPoly(QVal(i,:),TorqueDesire(i,:)',rU,Landa,Sat(i));    
    for j=1:length(QVal(i,:))
         QQ(j,:) = QVal(i,j).^(rU:-1:0)';
         DQ(j,:) = ([QVal(i,j).^(rU-1:-1:0) 0].*(rU:-1:0))';
         D2Q(j,:) = ([QVal(i,j).^(rU-2:-1:0) 0 0].*(rU:-1:0).*(rU-1:-1:-1))';
    end

    IntU2=IntU2 + ...
          Weight(i)* 1/2*(TorqueDesire(i,:)' - QQ*CoefBLSI )'* (TorqueDesire(i,:)' - QQ*CoefBLSI )*Tres;

    CostSlopeD1Q = CostSlopeD1Q + Weight(i)* 1/2*CoefBLSI'*(DQ'*DQ)*CoefBLSI*Tres;
    CostSlopeD2Q = CostSlopeD2Q + Weight(i)* 1/2*CoefBLSI'*(D2Q'*D2Q)*CoefBLSI*Tres;
end




BetaOptimal=[CoefBLS_UPassive1;CoefBLS_UPassive2];
TorquePassiveQ1valOptimal=polyval(BetaOptimal(0*(rU+1)+1:1*(rU+1)),Q1);
TorquePassiveQ2valOptimal=polyval(BetaOptimal(1*(rU+1)+1:2*(rU+1)),Q2val);
TorquePassiveValOptimal=[TorquePassiveQ1valOptimal; TorquePassiveQ2valOptimal];

 
figure('name','Desired Torque vs Angle')
    subplot(2,1,1)
    plot(rad2deg(q1(Middle:end)),torqueDesire(1,Middle:end),'linewidth',2)
    hold all
    plot(rad2deg(q1(Middle:end)),TorquePassiveQ1val,'linewidth',2)
    title('Initial Desired Torque-Angle Profile','FontWeight','normal','FontSize',16,'FontName','Times');
    xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('\tau_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
    subplot(2,1,2)
    plot((rad2deg(q2(Middle:end))),torqueDesire(2,Middle:end),'linewidth',2)
    hold all
    plot(rad2deg(q2(Middle:end)),TorquePassiveQ2val,'linewidth',2)
    xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('\tau_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    set(gca,'YMinorGrid','on')
    grid on
    
IntU2
CostSlopeD2Q
TotalCost_Opt = IntU2   + sum(Landa.*SelectLanda.*[CostSlopeD1Q CostSlopeD2Q])
TorqueActive=torqueDesire(:,Middle:end)-TorquePassiveVal;
IntAbsUdqActive=sum(sum(abs(TorqueActive.*Dq(:,Middle:end)),2))*Tres
IntAbsUdqActive=sum(sum(abs(torqueDesire(:,Middle:end).*Dq(:,Middle:end)),2))*Tres


%% 3DoF

clear
close
home

% envirment
g=9.81*0;

% parameters of robot
m=1;
L=1;

% EF motion
f=1;
A=1;
phi=0;
Tres=0.005;
time=0:Tres:10/f;

Start=ceil(length(time)/4);

% % Circle Motion
% Xef=A*cos(2*pi/f*time)+1.5;
% Yef=A*sin(2*pi/f*time)+1;
% 
% DXef=-(2*pi/f)*A*sin(2*pi/f*time);
% DYef= (2*pi/f)*A*cos(2*pi/f*time);
% 
% D2Xef=-(2*pi/f)^2*A*cos(2*pi/f*time);
% D2Yef=-(2*pi/f)^2*A*sin(2*pi/f*time);
% q1=deg2rad(0);
% q2=deg2rad(8.031);
% q3=deg2rad(51.317);


% Line motion: Horizontal
Xef=A*cos(2*pi*f*time+phi);
Yef=2.5*ones(size(time));
DXef=-(2*pi*f)*A*sin(2*pi*f*time+phi);
DYef= 2*zeros(size(time));
D2Xef=-(2*pi*f)^2*A*cos(2*pi*f*time+phi);
D2Yef= 2*zeros(size(time));
q1=deg2rad( 22.7023); % 41.7023 deg for y=2.5m  ,  22.2 deg for y=2.0m
q2=deg2rad(29.2663); % 18.2663 deg for y=2.5m  ,  29.468 deg for y=2.0m
q3=deg2rad(80.3377); % 44.3377 deg for y=2.5m  ,  71.431 deg for y=2.0m



figure('name','Path (3DoF)')
plot(Xef,Yef,'linewidth',2.5,'linestyle','-','color','g')

Pos=[Xef;Yef];
DPos=[DXef;DYef];
D2Pos=[D2Xef;D2Yef];

%%%%%%%%%%%%%%%%%%%%%%%%%

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
hold all
plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','-.','color','b')

%%%%%%%%%%%%%%%%%%%%%
OdeOpt= odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,7));
InitState=deg2rad([ 0 8.031 51.317 ,0 0 0, 0]);
DPos=[DXef;DYef];
D2Pos=[D2Xef;D2Yef];

Kp=10000;
Kd=200;

[T,Y] = ode15s(@(t,Y)SirDyn(t,Y,g,L,m,Pos,DPos,D2Pos,time,Kp,Kd), time,InitState,OdeOpt);

EnergyABS=Y(end,end)

q1=Y(:,1)';
q2=Y(:,2)';
q3=Y(:,3)';
RPos=L*[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
        sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];
plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','--','color','r')
xlabel('x')
ylabel('y')
hold off
axis equal
legend('Desired','Static path planing','dynamic path planing')

torqueDesire=TorqueCalculator(T,Y,g,m,m,m,L,L,L,Pos,DPos,D2Pos,time,Kp,Kd);
figure('name','Torque')
plot(time(Start:end),torqueDesire(:,Start:end))
legend('\tau_1','\tau_2','\tau_3')
xlabel('time (s)')
ylabel('\tau')
