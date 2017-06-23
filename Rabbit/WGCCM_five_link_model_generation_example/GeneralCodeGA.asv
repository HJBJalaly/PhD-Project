clear all
%% System Parameters
L_fem=.4; 
L_tib=L_fem; % According to selection of theta=c.q, femur and tibia should have equal lenght.
L_torso=.63;
Lc_fem=L_fem/2;
Lc_tib=L_tib/2;
Lc_torso=L_torso/2;
M_fem=6.8;
M_tib=3.2;
M_torso=12;
XX_fem=1;
XX_tib=1;
XX_torso=1;
g=9.8*1;
mu=.4;

%% Initialization
TT = [1   0   0   0  -1;
      0   1   0   0  -1;
     -1   0   1   0   0;
	  0  -1   0   1   0;
	  0   0   0   0   1];
c=[-1 0 -1/2 0 -1];

% Q_minus=TT*deg2rad([170 240 150 200 -10])'; % 
angl1=11;
angl2=6;
qTorso=-5;
Q_minusInit=TT*deg2rad([180-angl1+angl2 , 180+angl1+angl2 , 180-angl1+angl2-2*angl2 , 180+angl1+angl2-2*angl2 , qTorso])'; % This Q_minus
%  Q_minusInit=deg2rad([ 183.5  215.3  -20.1  -19.0   -9.5])';
% Q_minusInit=TT*deg2rad([240 , 120 , 120 , 90 , 0])' % This Q_minus
% Q_minusInit=deg2rad([183.0600  206.9524  -13.1551    -9.5970   -8.9210])';% extracted from a converged walking
% DQ_minusInit=10*TT*deg2rad([-10 10 -70 70 -1])'; % This Dq_minus
DQ_minusInit=50*deg2rad([-.1  .1  -2 1 -1 ])'; % This Dq_minus --> stable for 10 step
DQ_minusInit=1.5*deg2rad([ 5.5   -1.9  -65.0   27.8  -33.5])';
DQ_minusInit=deg2rad([ 50.82   -27.49  -91.8   16.57  -56.44])';% extracted from a converged walking
  % DQ_minus=50*deg2rad([ .1 -.2 -1 -1 -.5 ])'; 
% Impact and Relabaling

[ Q_plusInit,DQ_plusInit,V_tib2_plusInit,F2Init]=ImpactModel(Q_minusInit,1*DQ_minusInit,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
disp(['Q_plus=[ ',num2str(Q_plusInit'),']'])
disp(['DQ_plus=[ ',num2str(DQ_plusInit'),']'])
disp(['V_tib2_plus=[ ',num2str(V_tib2_plusInit'),']']);
disp(['F2=[ ',num2str(F2Init'),']']);
disp(['[F_n*mu,Ft]=[ ',num2str(F2Init(2)*mu),'>',num2str(F2Init(1)),']']);
disp('---------------------')

[KE_minusInit,PE_minusInit ] = KineticPotentialEnergy(Q_minusInit,DQ_minusInit,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[KE_plus,PE_plus ] = KineticPotentialEnergy(Q_plusInit,DQ_plusInit,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[p_tib2_minusInit ] = FoottPositon(Q_minusInit,L_fem, L_tib)';
[p_tib2_plusInit ] = FoottPositon(Q_plusInit,L_fem, L_tib)';
%% Bezier Coeffients

Ma=5;
c=[-1 0 -1/2 0 -1];
Theta_plusInit=c*Q_plusInit;
Theta_minusInit=c*Q_minusInit;

DTheta_plusInit=c*DQ_plusInit;
DTheta_minusInit=c*DQ_minusInit;
disp(['D_Theta_plus= ',num2str(DTheta_plusInit')])
disp(['D_Theta_minus= ',num2str(DTheta_minusInit')])
disp('---------------------')


AlfaInit=zeros(4,Ma+1);
AlfaInit(:,1)=Q_plusInit(1:4);
AlfaInit(:,Ma+1)=Q_minusInit(1:4);
AlfaInit(:,2)=DQ_plusInit(1:4)/(Ma*DTheta_plusInit)*(Theta_minusInit-Theta_plusInit)+AlfaInit(:,1);
AlfaInit(:,Ma)=-DQ_minusInit(1:4)/(Ma*DTheta_minusInit)*(Theta_minusInit-Theta_plusInit)+AlfaInit(:,Ma+1);

for ii=2+1:Ma-2+1
    AlfaInit(:,ii)=(1*AlfaInit(:,2)+1*AlfaInit(:,Ma-1+1))/2;
end


Hdr=zeros(4,1);
DHdr_s=zeros(4,1);
for Theta=linspace(Theta_plusInit,Theta_minusInit,50)
    [Hdr(:,end+1),DHdr_s(:,end+1)]=BezierFunction(Ma,Theta,AlfaInit,Theta_plusInit,Theta_minusInit);
end
Hdr(:,1)=[];
DHdr_s(:,1)=[];
SS=((linspace(Theta_plusInit,Theta_minusInit,50))-Theta_plusInit)/(Theta_minusInit-Theta_plusInit);
QQ=[Hdr;zeros(1,length(SS))];
DSS=[DHdr_s;zeros(1,length(SS))];
figure(2)
plot(c*QQ)
hold all
plot((c*diff(QQ')')./diff(SS))
grid on
legend('\theta','d(Hd_s)/ds')
figure(3)
subplot(2,2,1)
plot(SS,QQ(1,:))
legend('des');
ylabel('q_1')
subplot(2,2,2)
plot(SS,QQ(2,:))
legend('des');
ylabel('q_2')
subplot(2,2,3)
plot(SS,QQ(3,:))
legend('des');
ylabel('q_3')
subplot(2,2,4)
plot(SS,QQ(4,:))
legend('des');
ylabel('q_4')
    

% AnimBot3DOF(SS,QQ',1,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,1)

%% Ode Initial
LastMode='Initial';

Kp=30000;
Kd=350;

OptionsInit = odeset('RelTol',1e-3,'AbsTol',1e-3,'maxstep',1e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlfaInit,Theta_plusInit,Theta_minusInit,Ma,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));
            
[TcInit,SolCInit] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlfaInit,Theta_plusInit,Theta_minusInit,Ma,Kp,Kd),...
                    [0,2],[Q_plusInit; DQ_plusInit],OptionsInit);

Q_mnsInit =SolCInit(end,1:5 )';
DQ_mnsInit=SolCInit(end,6:10)';

q1=SolCInit(:,1)';q2=SolCInit(:,2)';q3=SolCInit(:,3)';q4=SolCInit(:,4)';q5=SolCInit(:,5)';
dq1=SolCInit(:,6)';dq2=SolCInit(:,7)';dq3=SolCInit(:,8)';dq4=SolCInit(:,9)';dq5=SolCInit(:,10)';

p_tib2_Init = [ + L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5),
                - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + L_tib*cos(q2 + q4 + q5)];

v_hip_Init =[dq1.*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) - dq2.*(L_fem*cos(q2 + q5) + L_tib*cos(q2 + q4 + q5)) + dq5.*(L_fem*cos(q1 + q5) - L_fem*cos(q2 + q5) + L_tib*cos(q1 + q3 + q5) - L_tib*cos(q2 + q4 + q5)) + L_tib*dq3.*cos(q1 + q3 + q5) - L_tib*dq4.*cos(q2 + q4 + q5)
             dq1.*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) - dq2.*(L_fem*sin(q2 + q5) + L_tib*sin(q2 + q4 + q5)) + dq5.*(L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5)) + L_tib*dq3.*sin(q1 + q3 + q5) - L_tib*dq4.*sin(q2 + q4 + q5)];

[Q_plsInit,DQ_plsInit,V_tib2_plsInit,F2Init]=ImpactModel(Q_mnsInit,DQ_mnsInit,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
                
MotionDataInit=ForceTorqueCalculator_5DoF([],[],'done');
        
T_impactInit=TcInit(end); 

[KE_mnsInit,PE_mnsInit ] = KineticPotentialEnergy(Q_mnsInit,DQ_minusInit,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[KE_plsInit,PE_plsInit ] = KineticPotentialEnergy(Q_plusInit,DQ_plusInit,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );

disp(['Q_plus=[ ',num2str(Q_plsInit'),'] (designed)'])
disp(['Q_pls=[ ',num2str(Q_plsInit'),'] (real)'])
disp(['DQ_plus=[ ',num2str(DQ_plusInit'),'] (designed)'])
disp(['DQ_pls=[ ',num2str(DQ_plsInit'),'] (real)'])
disp(['V_tib2_plus=[ ',num2str(V_tib2_plsInit'),']']);
disp(['K_mns= ',num2str(KE_mnsInit'),',  K_pls= ',num2str(KE_plsInit')])
disp(['CostAct=',num2str(sum((MotionDataInit(1,:).*MotionDataInit(1,:)+MotionDataInit(3,:).*MotionDataInit(3,:)+...
                                    MotionDataInit(2,:).*MotionDataInit(2,:)+MotionDataInit(4,:).*MotionDataInit(4,:))'.*[diff(TcInit) ;0]/2))])
disp(['Vel_hip_x_avg=',num2str(sum(v_hip_Init(1,1:end-1).*diff(TcInit'))/TcInit(end)),'m/s,    Pos_tib2_y_max=',num2str(max(p_tib2_Init(2,:))*100),'cm'])
disp('---------------------')
close all


% 
AnimBot3DOF(TcInit  ,SolCInit(:,1:5),T_impactInit,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,20)
ShowTimeOneStep(TcInit,SolCInit,MotionDataInit,c,Theta_plusInit,Theta_minusInit,AlfaInit,L_fem, L_tib, L_torso,mu,Ma);

%% GA  on only Gait
LastMode='GaitGa';

% Kp=30000;
% Kd=350;
% angl1=10;% deg
% angl2=6;% deg
% qTorso=-10;%deg
Initial=[reshape(AlfaInit(:,2+1:Ma-2+1),(Ma-3)*4,1); angl1;angl2;qTorso; DQ_minusInit];


PopulationSize=750;
PopInitRange=[-10*ones(1,(Ma-3)*4)  5  0  -45  -10*ones(1,5);
               10*ones(1,(Ma-3)*4) 45  45   0   10*ones(1,5)];

GenerationsLimit=100;
StallGenLimit=10;
nvars=(Ma-3)*4+3+5;
tic

InitialPopulation=repmat(Initial',PopulationSize,1)+1*randn(PopulationSize,nvars);

CostFun   = @(Param)CostFunctionGA(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Kp,Kd,mu);
NonConstr = @(Param)ConstraintFunctionGA(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,mu);

for i=11:20
    [Xx,fvalXx,exitflagXx,outputXx,populationXx,scoreXx] =...
        Op_GA(CostFun, NonConstr, nvars, InitialPopulation, PopInitRange, PopulationSize,GenerationsLimit,StallGenLimit);
     save(['TempMultiOde_Gait_',num2str(i),'.mat'])
end
toc
%% Ode on GA (optimized gait)
% Kp=100000;
% Kd=500;

angl1=Xx(end-7);
angl2=Xx(end-6);
qTorso=Xx(end-5);
Q_minusX=(TT*deg2rad([180-angl1+angl2 , 180+angl1+angl2 , 180-angl1+angl2-2*angl2 , 180+angl1+angl2-2*angl2 , qTorso])'); % This Q_minus;
DQ_minusX=Xx(end-4:end)';

[Q_plusX,DQ_plusX,V_tib2_plusX,F2X]=ImpactModel(Q_minusX,DQ_minusX,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);

Theta_plusX=c*Q_plusX;
Theta_minusX=c*Q_minusX;
DTheta_plusX=c*DQ_plusX;
DTheta_minusX=c*DQ_minusX;

AlfaX=zeros(4,Ma+1);

AlfaX(:,1)=Q_plusX(1:4);
AlfaX(:,Ma+1)=Q_minusX(1:4);
AlfaX(:,2)=DQ_plusX(1:4)/(Ma*DTheta_plusX)*(Theta_minusX-Theta_plusX)+AlfaX(:,1);
AlfaX(:,Ma)=-DQ_minusX(1:4)/(Ma*DTheta_minusX)*(Theta_minusX-Theta_plusX)+AlfaX(:,Ma+1);

AlfaX(:,2+1:Ma-2+1)=reshape(Xx(1:(Ma-3)*4),4,Ma-3);

OptionsX = odeset('RelTol',1e-3,'AbsTol',1e-3,'maxstep',1e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlfaX,Theta_plusX,Theta_minusX,Ma,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));
            
[TcX,SolCX] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlfaX,Theta_plusX,Theta_minusX,Ma,Kp,Kd),...
                    [0,2],[Q_plusX; DQ_plusX],OptionsX);

Q_mnsX =SolCX(end,1:5 )';
DQ_mnsX=SolCX(end,6:10)';

q1=SolCX(:,1)';q2=SolCX(:,2)';q3=SolCX(:,3)';q4=SolCX(:,4)';q5=SolCX(:,5)';
dq1=SolCX(:,6)';dq2=SolCX(:,7)';dq3=SolCX(:,8)';dq4=SolCX(:,9)';dq5=SolCX(:,10)';

p_tib2_X = [ + L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5),
                - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + L_tib*cos(q2 + q4 + q5)];

v_hip_X =[dq1.*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) - dq2.*(L_fem*cos(q2 + q5) + L_tib*cos(q2 + q4 + q5)) + dq5.*(L_fem*cos(q1 + q5) - L_fem*cos(q2 + q5) + L_tib*cos(q1 + q3 + q5) - L_tib*cos(q2 + q4 + q5)) + L_tib*dq3.*cos(q1 + q3 + q5) - L_tib*dq4.*cos(q2 + q4 + q5)
             dq1.*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) - dq2.*(L_fem*sin(q2 + q5) + L_tib*sin(q2 + q4 + q5)) + dq5.*(L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5)) + L_tib*dq3.*sin(q1 + q3 + q5) - L_tib*dq4.*sin(q2 + q4 + q5)];

[Q_plsX,DQ_plsX,V_tib2_plsX,F2X]=ImpactModel(Q_mnsX,DQ_mnsX,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);

Theta_plsX=c*Q_plsX;
Theta_mnsX=c*Q_mnsX;
DTheta_plsX=c*DQ_plsX;
DTheta_mnsX=c*DQ_mnsX;
                
MotionDataX=ForceTorqueCalculator_5DoF([],[],'done');

T_impactX=TcX(end); 


[KE_mnsX,PE_mnsX ] = KineticPotentialEnergy(Q_mnsX,DQ_mnsX,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[KE_plsX,PE_plsX ] = KineticPotentialEnergy(Q_plsX,DQ_plsX,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
disp('---------------------')
disp(['Q_pls=[ ',num2str(Q_plsX'),']'])
disp(['DQ_pls=[ ',num2str(DQ_plsX'),']'])
disp(['V_tib2_plus=[ ',num2str(V_tib2_plsX'),']']);
disp(['K_mns= ',num2str(KE_mnsX'),',  K_pls= ',num2str(KE_plsX')])
disp(['CostAct=',num2str(sum((MotionDataX(1,:).*MotionDataX(1,:)+MotionDataX(3,:).*MotionDataX(3,:)+...
                                    MotionDataX(2,:).*MotionDataX(2,:)+MotionDataX(4,:).*MotionDataX(4,:))'.*[diff(TcX) ;0]/2))])
disp(['Vel_hip_x_avg=',num2str(sum(v_hip_X(1,1:end-1).*diff(TcX'))/TcX(end)),'m/s,    Pos_tib2_y_max=',num2str(max(p_tib2_X(2,:))*100),'cm'])
disp('---------------------')

% stability
H0=eye(4,5);
H=[H0;c];
Hi=H^(-1);
qz=@(xi)Hi*[BezierFunction_hd_Fast(Ma,xi,AlfaX,Theta_plusX,Theta_minusX);xi];
subindex = @(q,r) q(r);
gamma0=@(qz)...
        [ XX_fem + XX_tib + 2*L_fem^2*M_fem + L_fem^2*M_tib + 2*L_tib^2*M_fem + L_fem^2*M_torso + 2*L_tib^2*M_tib + L_tib^2*M_torso + Lc_fem^2*M_fem + Lc_tib^2*M_tib - 2*L_fem*Lc_fem*M_fem - 2*L_tib*Lc_tib*M_tib - L_fem^2*M_tib*cos(qz(1) - qz(2)) + 4*L_fem*L_tib*M_fem*cos(qz(3)) + 2*L_fem*L_tib*M_tib*cos(qz(3)) + 2*L_fem*L_tib*M_torso*cos(qz(3)) - L_fem*L_torso*M_torso*cos(qz(1)) - 2*L_tib*Lc_fem*M_fem*cos(qz(3)) + L_fem*Lc_torso*M_torso*cos(qz(1)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3)), XX_fem + XX_tib + L_fem^2*M_tib + Lc_fem^2*M_fem + Lc_tib^2*M_tib - L_fem^2*M_tib*cos(qz(1) - qz(2)) + 2*L_fem*Lc_tib*M_tib*cos(qz(4)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)), XX_tib + 2*L_tib^2*M_fem + 2*L_tib^2*M_tib + L_tib^2*M_torso + Lc_tib^2*M_tib - 2*L_tib*Lc_tib*M_tib + 2*L_fem*L_tib*M_fem*cos(qz(3)) + L_fem*L_tib*M_tib*cos(qz(3)) + L_fem*L_tib*M_torso*cos(qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(3)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3)), XX_tib + Lc_tib^2*M_tib + L_fem*Lc_tib*M_tib*cos(qz(4)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)), 2*XX_fem + 2*XX_tib + XX_torso + 2*L_fem^2*M_fem + 2*L_fem^2*M_tib + 2*L_tib^2*M_fem + L_fem^2*M_torso + 2*L_tib^2*M_tib + L_tib^2*M_torso + L_torso^2*M_torso + 2*Lc_fem^2*M_fem + 2*Lc_tib^2*M_tib + Lc_torso^2*M_torso - 2*L_fem*Lc_fem*M_fem - 2*L_tib*Lc_tib*M_tib - 2*L_torso*Lc_torso*M_torso - 2*L_fem^2*M_tib*cos(qz(1) - qz(2)) + 4*L_fem*L_tib*M_fem*cos(qz(3)) + 2*L_fem*L_tib*M_tib*cos(qz(3)) + 2*L_fem*L_tib*M_torso*cos(qz(3)) - 2*L_fem*L_torso*M_torso*cos(qz(1)) - 2*L_tib*Lc_fem*M_fem*cos(qz(3)) + 2*L_fem*Lc_tib*M_tib*cos(qz(4)) + 2*L_fem*Lc_torso*M_torso*cos(qz(1)) - 2*L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - 2*L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - 2*L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - 2*L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - 2*L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - 2*L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + 2*L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3))];
kappa1=@(xi)c*inv([H0-BezierFunction_Dhd_Fast(Ma,xi,AlfaX,Theta_plusX,Theta_minusX)*c/(Theta_minusX-Theta_plusX);gamma0(qz(xi))])*[zeros(4,1);1];

% PE=g*(M_tib*(2*y1 - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + Lc_tib*cos(q2 + q4 + q5)) - M_fem*(cos(q1 + q5)*(L_fem - Lc_fem) - 2*y1 + L_tib*cos(q1 + q3 + q5)) + M_fem*(2*y1 - L_fem*cos(q1 + q5) + Lc_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5)) + M_torso*(2*y1 - L_fem*cos(q1 + q5) + cos(q5)*(L_torso - Lc_torso) - L_tib*cos(q1 + q3 + q5)) + M_tib*(2*y1 - cos(q1 + q3 + q5)*(L_tib - Lc_tib)))
% kappa2=-diff(PE,q5) ;
kappa2=@(qz)g*(M_fem*(sin(qz(4) + qz(5))*(L_fem - Lc_fem) + L_tib*sin(qz(4) + qz(3) + qz(5))) + M_tib*(L_fem*sin(qz(4) + qz(5)) - L_fem*sin(qz(2) + qz(5)) + L_tib*sin(qz(4) + qz(3) + qz(5)) - Lc_tib*sin(qz(2) + qz(4) + qz(5))) + M_fem*(L_fem*sin(qz(4) + qz(5)) - Lc_fem*sin(qz(2) + qz(5)) + L_tib*sin(qz(4) + qz(3) + qz(5))) + M_torso*(L_fem*sin(qz(4) + qz(5)) - sin(qz(5))*(L_torso - Lc_torso) + L_tib*sin(qz(4) + qz(3) + qz(5))) + M_tib*sin(qz(4) + qz(3) + qz(5))*(L_tib - Lc_tib));
Vzero=0;
VzeroMax=-inf;
Dxi=-(Theta_plusX-Theta_minusX)/50;
for Xi=linspace(Theta_plusX,Theta_minusX,50)
    Vzero=Vzero-kappa2(qz(Xi))/kappa1(Xi);
    if(Vzero>VzeroMax)
        VzeroMax=Vzero;
    end
end
Vzero=Vzero*Dxi;
VzeroMax=VzeroMax*Dxi;
    
DeltaDqMinus=@(Q_minus)DeltaDQ(Q_minus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
LandaHatDq=@(qm)inv([H0-Ma*(AlfaX(:,Ma+1)-AlfaX(:,Ma))*c/(Theta_minusX-Theta_plusX);gamma0(qm)])*[zeros(4,1);1];
deltaZero=gamma0(Q_plusX)*DeltaDqMinus(Q_minusX)*LandaHatDq(Q_minusX);

disp(['V_zero=    ',num2str(Vzero)])
disp(['V_zeroMax= ',num2str(VzeroMax)])
disp(['delta_zero=0<  ',num2str(deltaZero),'  <1'])
disp(['(deltaZ^2)/(1-deltaZ^2)*Vz+VzMax <0  :',num2str((deltaZero^2)/(1-deltaZero^2)*Vzero+VzeroMax)])
disp('---------------------')


% 
AnimBot3DOF(TcX  ,SolCX(:,1:5),[ T_impactX],L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,0)
ShowTimeOneStep(TcX,SolCX,MotionDataX,c,Theta_plusX,Theta_mnsX,AlfaX,L_fem, L_tib, L_torso,mu,Ma);

%% Multi Ode on GA (optimized gait)

% Kp=4000;
% Kd=130;
Stride=5;

AlphaMultiX=AlfaX;

Q_mnsMultiX=Q_minusX;
DQ_mnsMultiX=DQ_minusX;
Q_plsMultiX=Q_plusX;
DQ_plsMultiX=DQ_plusX;
Theta_mnsMultiX= Theta_minusX;
Theta_plsMultiX= Theta_plusX;

T_impactMultiX=[]; 
QqMultiX=[];
DQqMultiX=[];
TimeMultiX=[0];
MotionDataMultiX=[];
p_tib2_MulitX =[];
v_hip_MulitX =[];
    
for step=1:2*Stride
    disp('----------------------')
    disp(['step=',num2str(step)])
    
    
    OptionsMultiX = odeset('RelTol',1e-3,'AbsTol',1e-3,'maxstep',1e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlphaMultiX,Theta_plsMultiX,Theta_mnsMultiX,Ma,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));

    [TcMultiX,SolCMultiX] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlphaMultiX,Theta_plsMultiX,Theta_mnsMultiX,Ma,Kp,Kd),...
                        [0,2],[Q_plsMultiX; DQ_plsMultiX],OptionsMultiX);

    Q_mnsMultiX =SolCMultiX(end,1:5 )';
    DQ_mnsMultiX=SolCMultiX(end,6:10)';

    [ Q_plsMultiX,DQ_plsMultiX,V_tib2_pls_MulitX,F2_MulitX]=ImpactModel(Q_mnsMultiX,DQ_mnsMultiX,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
%     disp(['Q_plus=[ ',num2str(Q_pls'),']']);
    disp(['DQ_plus=[ ',num2str(DQ_plsMultiX'),']']);
    disp(['F2=[ ',num2str(F2_MulitX'),']  (max Ft=',num2str(F2_MulitX(2)*mu),')']);
%      
    Theta_plsMultiX=c*Q_plsMultiX;
    Theta_mnsMultiX=c*Q_mnsMultiX;

    DTheta_plsMultiX=c*DQ_plsMultiX;
    DTheta_mnsMultiX=c*DQ_mnsMultiX;

%     AlphaMulti=zeros(4,Ma+1);
%      AlphaMultiX(:,1)=Q_plsMultiX(1:4);
%      AlphaMultiX(:,Ma+1)=Q_mnsMultiX(1:4);
%      AlphaMultiX(:,2)=DQ_plsMultiX(1:4)/(Ma*DTheta_plsMultiX)*(Theta_mnsMultiX-Theta_plsMultiX)+AlphaMultiX(:,1);
%     AlphaMultiX(:,Ma)=-DQ_mnsMultiX(1:4)/(Ma*DTheta_mnsMultiX)*(Theta_mnsMultiX-Theta_plsMultiX)+AlphaMultiX(:,Ma+1);
%     
%     for i=2+1:Ma-2+1
%         AlphaMulti(:,i)=(0*AlphaMulti(:,2)+50*AlphaMulti(:,Ma-1+1))/50;
%     end


    MotionDataMultiX=[MotionDataMultiX, ForceTorqueCalculator_5DoF([],[],'done')];

    q1=SolCMultiX(:,1)';q2=SolCMultiX(:,2)';q3=SolCMultiX(:,3)';q4=SolCMultiX(:,4)';q5=SolCMultiX(:,5)';
    dq1=SolCMultiX(:,6)';dq2=SolCMultiX(:,7)';dq3=SolCMultiX(:,8)';dq4=SolCMultiX(:,9)';dq5=SolCMultiX(:,10)';

    p_tib2_MulitX =[p_tib2_MulitX,  [ + L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5),
                                    - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + L_tib*cos(q2 + q4 + q5)]];

    v_hip_MulitX =[v_hip_MulitX , [dq1.*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) - dq2.*(L_fem*cos(q2 + q5) + L_tib*cos(q2 + q4 + q5)) + dq5.*(L_fem*cos(q1 + q5) - L_fem*cos(q2 + q5) + L_tib*cos(q1 + q3 + q5) - L_tib*cos(q2 + q4 + q5)) + L_tib*dq3.*cos(q1 + q3 + q5) - L_tib*dq4.*cos(q2 + q4 + q5)
                                 dq1.*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) - dq2.*(L_fem*sin(q2 + q5) + L_tib*sin(q2 + q4 + q5)) + dq5.*(L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5)) + L_tib*dq3.*sin(q1 + q3 + q5) - L_tib*dq4.*sin(q2 + q4 + q5)]];

    QqMultiX =[QqMultiX ;SolCMultiX(:,1:5)];
    DQqMultiX=[DQqMultiX ;SolCMultiX(:,6:10)];
    TimeMultiX=[TimeMultiX; TimeMultiX(end)+TcMultiX];
    T_impactMultiX=[T_impactMultiX TimeMultiX(end)];
    
end


% stabolity
qzM=@(xi)Hi*[BezierFunction_hd_Fast(Ma,xi,AlphaMultiX,Theta_plsMultiX,Theta_mnsMultiX);xi];
gamma0M=@(qz)...
        [ XX_fem + XX_tib + 2*L_fem^2*M_fem + L_fem^2*M_tib + 2*L_tib^2*M_fem + L_fem^2*M_torso + 2*L_tib^2*M_tib + L_tib^2*M_torso + Lc_fem^2*M_fem + Lc_tib^2*M_tib - 2*L_fem*Lc_fem*M_fem - 2*L_tib*Lc_tib*M_tib - L_fem^2*M_tib*cos(qz(1) - qz(2)) + 4*L_fem*L_tib*M_fem*cos(qz(3)) + 2*L_fem*L_tib*M_tib*cos(qz(3)) + 2*L_fem*L_tib*M_torso*cos(qz(3)) - L_fem*L_torso*M_torso*cos(qz(1)) - 2*L_tib*Lc_fem*M_fem*cos(qz(3)) + L_fem*Lc_torso*M_torso*cos(qz(1)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3)), XX_fem + XX_tib + L_fem^2*M_tib + Lc_fem^2*M_fem + Lc_tib^2*M_tib - L_fem^2*M_tib*cos(qz(1) - qz(2)) + 2*L_fem*Lc_tib*M_tib*cos(qz(4)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)), XX_tib + 2*L_tib^2*M_fem + 2*L_tib^2*M_tib + L_tib^2*M_torso + Lc_tib^2*M_tib - 2*L_tib*Lc_tib*M_tib + 2*L_fem*L_tib*M_fem*cos(qz(3)) + L_fem*L_tib*M_tib*cos(qz(3)) + L_fem*L_tib*M_torso*cos(qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(3)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3)), XX_tib + Lc_tib^2*M_tib + L_fem*Lc_tib*M_tib*cos(qz(4)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)), 2*XX_fem + 2*XX_tib + XX_torso + 2*L_fem^2*M_fem + 2*L_fem^2*M_tib + 2*L_tib^2*M_fem + L_fem^2*M_torso + 2*L_tib^2*M_tib + L_tib^2*M_torso + L_torso^2*M_torso + 2*Lc_fem^2*M_fem + 2*Lc_tib^2*M_tib + Lc_torso^2*M_torso - 2*L_fem*Lc_fem*M_fem - 2*L_tib*Lc_tib*M_tib - 2*L_torso*Lc_torso*M_torso - 2*L_fem^2*M_tib*cos(qz(1) - qz(2)) + 4*L_fem*L_tib*M_fem*cos(qz(3)) + 2*L_fem*L_tib*M_tib*cos(qz(3)) + 2*L_fem*L_tib*M_torso*cos(qz(3)) - 2*L_fem*L_torso*M_torso*cos(qz(1)) - 2*L_tib*Lc_fem*M_fem*cos(qz(3)) + 2*L_fem*Lc_tib*M_tib*cos(qz(4)) + 2*L_fem*Lc_torso*M_torso*cos(qz(1)) - 2*L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - 2*L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - 2*L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - 2*L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - 2*L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - 2*L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + 2*L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3))];
kappa1M=@(xi)c*inv([H0-BezierFunction_Dhd_Fast(Ma,xi,AlphaMultiX,Theta_plsMultiX,Theta_mnsMultiX)*c/(Theta_mnsMultiX-Theta_plsMultiX);gamma0M(qz(xi))])*[zeros(4,1);1];
VzeroM=0;
VzeroMaxM=-inf;
Dxi=-(Theta_plsMultiX-Theta_mnsMultiX)/50;
for Xi=linspace(Theta_plsMultiX,Theta_mnsMultiX,50)
    VzeroM=VzeroM-kappa2(qzM(Xi))/kappa1M(Xi);
    if(VzeroM>VzeroMaxM)
        VzeroMaxM=VzeroM;
    end
end
VzeroM=VzeroM*Dxi;
VzeroMaxM=VzeroMaxM*Dxi;
    
DeltaDqMinusM=@(Q_mnsMulti)DeltaDQ(Q_mnsMulti,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
LandaHatDqM=@(qm)inv([H0-Ma*(AlphaMultiX(:,Ma+1)-AlphaMultiX(:,Ma))*c/(Theta_mnsMultiX-Theta_plsMultiX);gamma0M(qm)])*[zeros(4,1);1];
deltaZeroM=gamma0M(Q_plsMultiX)*DeltaDqMinusM(Q_mnsMultiX)*LandaHatDqM(Q_mnsMultiX);

disp(['V_zero= ',num2str(VzeroM),',  V_zMax= ',num2str(VzeroMaxM),',  delta_zero=0< ',num2str(deltaZeroM),' <1'])
disp(['(deltaZ^2)/(1-deltaZ^2)*Vz+VzMax <0  :',num2str((deltaZeroM^2)/(1-deltaZeroM^2)*VzeroM+VzeroMaxM)])
disp(['Vel_hip_x_avg=',num2str(sum(v_hip_MulitX(1,end-length(TcMultiX)+2:end).*diff(TcMultiX'))/TcMultiX(end)),'m/s,    Pos_tib2_y_max=',num2str(max(p_tib2_MulitX(2,:))*100),'cm'])
disp(['ang1=',num2str(angl1),', ang2=',num2str(angl2)])
disp('---------------------')

TimeMultiX(1)=[];

ShowTimeMultiStep(TimeMultiX,T_impactMultiX,QqMultiX,DQqMultiX,MotionDataMultiX,p_tib2_MulitX,v_hip_MulitX,mu)
AnimBot3DOF( TimeMultiX, QqMultiX,T_impactMultiX,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,0)

%% Optimization damper-spring on Multi GA (optimized gait + Sequential Optimization)
Weight=[1 1];
SampleRate=10;   
rM=3;
rB=0;
rD=0-1;
Landa=0.00000005;
Gamma=0.0005;

ii=length(T_impactMultiX)/2;
indx1=find(TimeMultiX==T_impactMultiX(2*(ii)-2),1)+1;
indx2=find(TimeMultiX==T_impactMultiX(2*(ii)-1),1);
indx3=find(TimeMultiX==T_impactMultiX(2*(ii)-0),1);

TimeStrideMultiX=TimeMultiX(indx1:indx3);

Ur1MultiX=MotionDataMultiX(1,indx1:indx2)';
Ur1MultiX=[Ur1MultiX;MotionDataMultiX(2,indx2+1:indx3)'];
Ur3MultiX=MotionDataMultiX(3,indx1:indx2)';
Ur3MultiX=[Ur3MultiX;MotionDataMultiX(4,indx2+1:indx3)'];

Q1MultiX=QqMultiX(indx1:indx2,1);
Q1MultiX=[Q1MultiX;QqMultiX(indx2+1:indx3,2)];
Q3MultiX=QqMultiX(indx1:indx2,3);
Q3MultiX=[Q3MultiX;QqMultiX(indx2+1:indx3,4)];

Qhat13MultiX=Q1MultiX+Q3MultiX;

DtQ1MultiX=DQqMultiX(indx1:indx2,1);
DtQ1MultiX=[DtQ1MultiX;DQqMultiX(indx2+1:indx3,2)];
DtQ3MultiX=DQqMultiX(indx1:indx2  ,3);
DtQ3MultiX=[DtQ3MultiX;DQqMultiX(indx2+1:indx3,4)];

[CostActuationMultiX,CostD2QMultiX,CostParaRegMultiX,BetaOptimalMultiX,ThetaOptimalMultiX,EtaOptimalMultiX,TorqueMonoSpring1MultiX,TorqueMonoSpring3MultiX,TorqueDamper1MultiX,TorqueDamper3MultiX,TorqueBiSpring13MultiX]=...
                    OptimalParam(TimeStrideMultiX,Q1MultiX,Q3MultiX,Qhat13MultiX,DtQ1MultiX,DtQ3MultiX,Ur1MultiX,Ur3MultiX,rM,rB,rD,Landa,Gamma,Weight,SampleRate,1,0);
                
Ua1MultiX=(Ur1MultiX -TorqueMonoSpring1MultiX -TorqueDamper1MultiX -TorqueBiSpring13MultiX );
Ua3MultiX=(Ur3MultiX -TorqueMonoSpring3MultiX -TorqueDamper3MultiX -TorqueBiSpring13MultiX );


Work_rMultiX=sum(abs(Ur1MultiX.*DtQ1MultiX).*[diff(TimeStrideMultiX) ;0])+sum(abs(Ur3MultiX.*DtQ3MultiX).*[diff(TimeStrideMultiX) ;0]);
Work_aMultiX=sum(abs(Ua1MultiX.*DtQ1MultiX).*[diff(TimeStrideMultiX) ;0])+sum(abs(Ua3MultiX.*DtQ3MultiX).*[diff(TimeStrideMultiX) ;0]);

CostReqMultiX=sum((Ur1MultiX.*Ur1MultiX+Ur3MultiX.*Ur3MultiX).*[diff(TimeStrideMultiX);0])/2;

Title=sprintf('%22s  % 12s % 12s'  ,    'Cost','Req. Work','Act. Work');
Result_InitMultiX=sprintf('%-11s %11.2f %10.2f %10.2f',   'Initial:',CostReqMultiX, Work_rMultiX,Work_rMultiX);
Result_OptMultiX =sprintf('%-11s %11.2f %10.2f %10.2f \n',   'Optimized:',CostActuationMultiX, Work_rMultiX,Work_aMultiX);

disp(Title)
disp(Result_InitMultiX)
disp(Result_OptMultiX)
%% Simultinous GA
clear *Xx
close all
LastMode='SimGaitGa';
% Kp=30000;
% Kd=350;
% angl1=10;% deg
% angl2=6;% deg
% qTorso=-10;%deg

Initial=[reshape(AlfaInit(:,2+1:Ma-2+1),(Ma-3)*4,1); angl1;angl2;qTorso; DQ_minusInit];


PopulationSize=750;
PopInitRange=[-10*ones(1,(Ma-3)*4)  5  0  -45  -10*ones(1,5);
               10*ones(1,(Ma-3)*4) 45  45   0   10*ones(1,5)];

GenerationsLimit=100;
StallGenLimit=10;
nvars=(Ma-3)*4+3+5;
tic

rM=3;
rB=3;
rD=0;
for i=11:20
    InitialPopulation=repmat(Initial',PopulationSize,1)+1*randn(PopulationSize,nvars);

    CostFunSim   = @(Param)CostFunctionSimGA(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Kp,Kd,mu,rM,rB,rD,Landa,Gamma,Weight,SampleRate);
    NonConstr = @(Param)ConstraintFunctionGA(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,mu);

    [XxSim,fvalXxSim,exitflagXxSim,outputXxSim,populationXxSim,scoreXxSim] =...
        Op_GA(CostFunSim, NonConstr, nvars, InitialPopulation, PopInitRange, PopulationSize,GenerationsLimit,StallGenLimit);
    save(['TempMultiOde_SimGait_MBD_',num2str(i),'.mat'])
end


%% Ode on Sim GA
angl1Sim=XxSim(end-7);
angl2Sim=XxSim(end-6);
qTorsoSim=XxSim(end-5);
Q_minusXSim=(TT*deg2rad([180-angl1Sim+angl2Sim , 180+angl1Sim+angl2Sim , 180-angl1Sim+angl2Sim-2*angl2Sim , 180+angl1Sim+angl2Sim-2*angl2Sim , qTorsoSim])'); % This Q_minus;
DQ_minusXSim=XxSim(end-4:end)';

[Q_plusXSim,DQ_plusXSim,V_tib2_plusXSim,F2XSim]=ImpactModel(Q_minusXSim,DQ_minusXSim,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);

Theta_plusXSim=c*Q_plusXSim;
Theta_minusXSim=c*Q_minusXSim;
DTheta_plusXSim=c*DQ_plusXSim;
DTheta_minusXSim=c*DQ_minusXSim;

AlfaXSim=zeros(4,Ma+1);

AlfaXSim(:,1)=Q_plusXSim(1:4);
AlfaXSim(:,Ma+1)=Q_minusXSim(1:4);
AlfaXSim(:,2)=DQ_plusXSim(1:4)/(Ma*DTheta_plusXSim)*(Theta_minusXSim-Theta_plusXSim)+AlfaXSim(:,1);
AlfaXSim(:,Ma)=-DQ_minusXSim(1:4)/(Ma*DTheta_minusXSim)*(Theta_minusXSim-Theta_plusXSim)+AlfaXSim(:,Ma+1);

AlfaXSim(:,2+1:Ma-2+1)=reshape(XxSim(1:(Ma-3)*4),4,Ma-3);


OptionsXSim = odeset('RelTol',1e-3,'AbsTol',1e-3,'maxstep',1e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlfaXSim,Theta_plusXSim,Theta_minusXSim,Ma,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));
            
[TcXSim,SolCXSim] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlfaXSim,Theta_plusXSim,Theta_minusXSim,Ma,Kp,Kd),...
                    [0,2],[Q_plusXSim; DQ_plusXSim],OptionsXSim);

Q_mnsXSim =SolCXSim(end,1:5 )';
DQ_mnsXSim=SolCXSim(end,6:10)';

q1=SolCXSim(:,1)';q2=SolCXSim(:,2)';q3=SolCXSim(:,3)';q4=SolCXSim(:,4)';q5=SolCXSim(:,5)';
dq1=SolCXSim(:,6)';dq2=SolCXSim(:,7)';dq3=SolCXSim(:,8)';dq4=SolCXSim(:,9)';dq5=SolCXSim(:,10)';

p_tib2_XSim = [ + L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5),
                - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + L_tib*cos(q2 + q4 + q5)];

v_hip_XSim =[dq1.*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) - dq2.*(L_fem*cos(q2 + q5) + L_tib*cos(q2 + q4 + q5)) + dq5.*(L_fem*cos(q1 + q5) - L_fem*cos(q2 + q5) + L_tib*cos(q1 + q3 + q5) - L_tib*cos(q2 + q4 + q5)) + L_tib*dq3.*cos(q1 + q3 + q5) - L_tib*dq4.*cos(q2 + q4 + q5)
             dq1.*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) - dq2.*(L_fem*sin(q2 + q5) + L_tib*sin(q2 + q4 + q5)) + dq5.*(L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5)) + L_tib*dq3.*sin(q1 + q3 + q5) - L_tib*dq4.*sin(q2 + q4 + q5)];


[Q_plsXSim,DQ_plsXSim,V_tib2_plsXSim,F2XSim]=ImpactModel(Q_mnsXSim,DQ_mnsXSim,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);

Theta_plsXSim=c*Q_plsXSim;
Theta_mnsXSim=c*Q_mnsXSim;
DTheta_plsXSim=c*DQ_plsXSim;
DTheta_mnsXSim=c*DQ_mnsXSim;
                
MotionDataXSim=ForceTorqueCalculator_5DoF([],[],'done');

T_impactXSim=[TcXSim(end)]; 

[KE_mnsXSim,PE_mnsXSim ] = KineticPotentialEnergy(Q_mnsXSim,DQ_mnsXSim,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[KE_plsXSim,PE_plsXSim ] = KineticPotentialEnergy(Q_plsXSim,DQ_plsXSim,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );

disp(['Q_pls=[ ',num2str(Q_plsXSim'),']'])
disp(['DQ_pls=[ ',num2str(DQ_plsXSim'),']'])
disp(['V_tib2_plus=[ ',num2str(V_tib2_plsXSim'),']']);
disp(['K_mns= ',num2str(KE_mnsXSim'),',  K_pls= ',num2str(KE_plsXSim')])
disp(['CostAct=',num2str(sum((MotionDataXSim(1,:).*MotionDataXSim(1,:)+MotionDataXSim(3,:).*MotionDataXSim(3,:)+...
                                    MotionDataXSim(2,:).*MotionDataXSim(2,:)+MotionDataXSim(4,:).*MotionDataXSim(4,:))'.*[diff(TcXSim) ;0]/2))])
disp(['Vel_hip_x_avg=',num2str(mean(v_hip_XSim(1,:))),'m/s,    Pos_tib2_y_max=',num2str(max(p_tib2_XSim(2,:))*100),'cm'])

% stability
H0=eye(4,5);
H=[H0;c];
Hi=H^(-1);
qz=@(xi)Hi*[BezierFunction_hd_Fast(Ma,xi,AlfaXSim,Theta_plusXSim,Theta_minusXSim);xi];
subindex = @(q,r) q(r);
gamma0=@(qz)...
        [ XX_fem + XX_tib + 2*L_fem^2*M_fem + L_fem^2*M_tib + 2*L_tib^2*M_fem + L_fem^2*M_torso + 2*L_tib^2*M_tib + L_tib^2*M_torso + Lc_fem^2*M_fem + Lc_tib^2*M_tib - 2*L_fem*Lc_fem*M_fem - 2*L_tib*Lc_tib*M_tib - L_fem^2*M_tib*cos(qz(1) - qz(2)) + 4*L_fem*L_tib*M_fem*cos(qz(3)) + 2*L_fem*L_tib*M_tib*cos(qz(3)) + 2*L_fem*L_tib*M_torso*cos(qz(3)) - L_fem*L_torso*M_torso*cos(qz(1)) - 2*L_tib*Lc_fem*M_fem*cos(qz(3)) + L_fem*Lc_torso*M_torso*cos(qz(1)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3)), XX_fem + XX_tib + L_fem^2*M_tib + Lc_fem^2*M_fem + Lc_tib^2*M_tib - L_fem^2*M_tib*cos(qz(1) - qz(2)) + 2*L_fem*Lc_tib*M_tib*cos(qz(4)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)), XX_tib + 2*L_tib^2*M_fem + 2*L_tib^2*M_tib + L_tib^2*M_torso + Lc_tib^2*M_tib - 2*L_tib*Lc_tib*M_tib + 2*L_fem*L_tib*M_fem*cos(qz(3)) + L_fem*L_tib*M_tib*cos(qz(3)) + L_fem*L_tib*M_torso*cos(qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(3)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3)), XX_tib + Lc_tib^2*M_tib + L_fem*Lc_tib*M_tib*cos(qz(4)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)), 2*XX_fem + 2*XX_tib + XX_torso + 2*L_fem^2*M_fem + 2*L_fem^2*M_tib + 2*L_tib^2*M_fem + L_fem^2*M_torso + 2*L_tib^2*M_tib + L_tib^2*M_torso + L_torso^2*M_torso + 2*Lc_fem^2*M_fem + 2*Lc_tib^2*M_tib + Lc_torso^2*M_torso - 2*L_fem*Lc_fem*M_fem - 2*L_tib*Lc_tib*M_tib - 2*L_torso*Lc_torso*M_torso - 2*L_fem^2*M_tib*cos(qz(1) - qz(2)) + 4*L_fem*L_tib*M_fem*cos(qz(3)) + 2*L_fem*L_tib*M_tib*cos(qz(3)) + 2*L_fem*L_tib*M_torso*cos(qz(3)) - 2*L_fem*L_torso*M_torso*cos(qz(1)) - 2*L_tib*Lc_fem*M_fem*cos(qz(3)) + 2*L_fem*Lc_tib*M_tib*cos(qz(4)) + 2*L_fem*Lc_torso*M_torso*cos(qz(1)) - 2*L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - 2*L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - 2*L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - 2*L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - 2*L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - 2*L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + 2*L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3))];
kappa1=@(xi)c*inv([H0-BezierFunction_Dhd_Fast(Ma,xi,AlfaXSim,Theta_plusXSim,Theta_minusXSim)*c/(Theta_minusXSim-Theta_plusXSim);gamma0(qz(xi))])*[zeros(4,1);1];

% PE=g*(M_tib*(2*y1 - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + Lc_tib*cos(q2 + q4 + q5)) - M_fem*(cos(q1 + q5)*(L_fem - Lc_fem) - 2*y1 + L_tib*cos(q1 + q3 + q5)) + M_fem*(2*y1 - L_fem*cos(q1 + q5) + Lc_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5)) + M_torso*(2*y1 - L_fem*cos(q1 + q5) + cos(q5)*(L_torso - Lc_torso) - L_tib*cos(q1 + q3 + q5)) + M_tib*(2*y1 - cos(q1 + q3 + q5)*(L_tib - Lc_tib)))
% kappa2=-diff(PE,q5) ;
kappa2=@(qz)g*(M_fem*(sin(qz(4) + qz(5))*(L_fem - Lc_fem) + L_tib*sin(qz(4) + qz(3) + qz(5))) + M_tib*(L_fem*sin(qz(4) + qz(5)) - L_fem*sin(qz(2) + qz(5)) + L_tib*sin(qz(4) + qz(3) + qz(5)) - Lc_tib*sin(qz(2) + qz(4) + qz(5))) + M_fem*(L_fem*sin(qz(4) + qz(5)) - Lc_fem*sin(qz(2) + qz(5)) + L_tib*sin(qz(4) + qz(3) + qz(5))) + M_torso*(L_fem*sin(qz(4) + qz(5)) - sin(qz(5))*(L_torso - Lc_torso) + L_tib*sin(qz(4) + qz(3) + qz(5))) + M_tib*sin(qz(4) + qz(3) + qz(5))*(L_tib - Lc_tib));
Vzero=0;
VzeroMax=-inf;
Dxi=-(Theta_plusXSim-Theta_minusXSim)/50;
for Xi=linspace(Theta_plusXSim,Theta_minusXSim,50)
    Vzero=Vzero-kappa2(qz(Xi))/kappa1(Xi);
    if(Vzero>VzeroMax)
        VzeroMax=Vzero;
    end
end
Vzero=Vzero*Dxi;
VzeroMax=VzeroMax*Dxi;
    
DeltaDqMinus=@(Q_minus)DeltaDQ(Q_minus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
LandaHatDq=@(qm)inv([H0-Ma*(AlfaXSim(:,Ma+1)-AlfaXSim(:,Ma))*c/(Theta_minusXSim-Theta_plusXSim);gamma0(qm)])*[zeros(4,1);1];
deltaZero=gamma0(Q_plusXSim)*DeltaDqMinus(Q_minusXSim)*LandaHatDq(Q_minusXSim);

disp(['V_zero=    ',num2str(Vzero),',    V_zeroMax= ',num2str(VzeroMax)])
disp(['delta_zero=0<  ',num2str(deltaZero),'  <1'])
disp(['(deltaZ^2)/(1-deltaZ^2)*Vz+VzMax <0  :',num2str((deltaZero^2)/(1-deltaZero^2)*Vzero+VzeroMax)])

% 
AnimBot3DOF(TcXSim  ,SolCXSim(:,1:5), T_impactXSim,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,0)
ShowTimeOneStep(TcXSim,SolCXSim,MotionDataXSim,c,Theta_plusXSim,Theta_mnsXSim,AlfaXSim,L_fem, L_tib, L_torso,mu,Ma);


Ur1XSim=MotionDataXSim(1,:)';
Ur1XSim=[Ur1XSim;MotionDataXSim(2,:)'];
Ur3XSim=MotionDataXSim(3,:)';
Ur3XSim=[Ur3XSim;MotionDataXSim(4,:)'];

Q1XSim=SolCXSim(:,1);
Q1XSim=[Q1XSim;SolCXSim(:,2)];
Q3XSim=SolCXSim(:,3);
Q3XSim=[Q3XSim;SolCXSim(:,4)];

Qhat13XSim=Q1XSim+Q3XSim;

DtQ1XSim=SolCXSim(:,6);
DtQ1XSim=[DtQ1XSim;SolCXSim(:,7)];
DtQ3XSim=SolCXSim(:,8);
DtQ3XSim=[DtQ3XSim;SolCXSim(:,9)];

[CostActuationXSim,CostD2QXSim,CostParaRegXSim,BetaOptimalXSim,ThetaOptimalXSim,EtaOptimalXSim,TorqueMonoSpring1XSim,TorqueMonoSpring3XSim,TorqueDamper1XSim,TorqueDamper3XSim,TorqueBiSpring13XSim]=...
                    OptimalParam([TcXSim;TcXSim+TcXSim(end)],Q1XSim,Q3XSim,Qhat13XSim,DtQ1XSim,DtQ3XSim,Ur1XSim,Ur3XSim,rM,rB,rD,Landa,Gamma,Weight,SampleRate,1,0);
                
Ua1XSim=(Ur1XSim -TorqueMonoSpring1XSim -TorqueDamper1XSim -TorqueBiSpring13XSim );
Ua3XSim=(Ur3XSim -TorqueMonoSpring3XSim -TorqueDamper3XSim -TorqueBiSpring13XSim );


Work_rXSim=sum(abs(Ur1XSim.*DtQ1XSim).*[diff([TcXSim;TcXSim+TcXSim(end)]) ;0])+sum(abs(Ur3XSim.*DtQ3XSim).*[diff([TcXSim;TcXSim+TcXSim(end)]) ;0]);
Work_a=sum(abs(Ua1XSim.*DtQ1XSim).*[diff([TcXSim;TcXSim+TcXSim(end)]) ;0])+sum(abs(Ua3XSim.*DtQ3XSim).*[diff([TcXSim;TcXSim+TcXSim(end)]) ;0]);

CostReq=sum((Ur1XSim.*Ur1XSim+Ur3XSim.*Ur3XSim).*[diff([TcXSim;TcXSim+TcXSim(end)]);0])/2;

Title=sprintf('%22s  % 12s % 12s'  ,    'Cost','Req. Work','Act. Work');
Result_Init_XSim=sprintf('%-11s %11.2f %10.2f %10.2f',   'Initial:',CostReq, Work_rXSim,Work_rXSim);
Result_Opt_XSim =sprintf('%-11s %11.2f %10.2f %10.2f',   'Optimized:',CostActuationXSim, Work_rXSim,Work_a);

disp(Title)
disp(Result_Init_XSim)
disp(Result_Opt_XSim)


%% Multi Ode on Sim GA

% Kp=4000;
% Kd=130;

Stride=5;
AlphaMultiXSim=AlfaXSim;
Q_mnsMultiXSim=Q_minusXSim;
DQ_mnsMultiXSim=DQ_minusXSim;
Q_plsMultiXSim=Q_plusXSim;
DQ_plsMultiXSim=DQ_plusXSim;
Theta_mnsMultiXSim= Theta_minusXSim;
Theta_plsMultiXSim= Theta_plusXSim;

T_impactMultiXSim=[]; 
QqMultiXSim=[];
DQqMultiXSim=[];
TimeMultiXSim=[0];
MotionDataMultiXSim=[];
p_tib2_MulitXSim =[];
v_hip_MulitXSim =[];
    
for step=1:2*Stride
    disp('----------------------')
    disp(['step=',num2str(step)])
    
    OptionsMultiXSim = odeset('RelTol',1e-3,'AbsTol',1e-3,'maxstep',1e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlphaMultiXSim,Theta_plsMultiXSim,Theta_mnsMultiXSim,Ma,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));

    [TcMultiXSim,SolCMultiXSim] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlphaMultiXSim,Theta_plsMultiXSim,Theta_mnsMultiXSim,Ma,Kp,Kd),...
                        [0,2],[Q_plsMultiXSim; DQ_plsMultiXSim],OptionsMultiXSim);
    
    Q_mnsMultiXSim=SolCMultiXSim(end,1:5 )';
    DQ_mnsMultiXSim=SolCMultiXSim(end,6:10)';

    [Q_plsMultiXSim,DQ_plsMultiXSim,V_tib2_pls_MulitXSim,F2_MulitXSim]=ImpactModel(Q_mnsMultiXSim,DQ_mnsMultiXSim,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
%     disp(['Q_plus=[ ',num2str(Q_pls'),']']);
    disp(['DQ_plus=[ ',num2str(DQ_plsMultiXSim'),']']);
    disp(['F2=[ ',num2str(F2_MulitXSim'),']  (max Ft=',num2str(F2_MulitXSim(2)*mu),')']);

    Theta_plsMultiXSim=c*Q_plsMultiXSim;
    Theta_mnsMultiXSim=c*Q_mnsMultiXSim;

    DTheta_plsMultiXSim=c*DQ_plsMultiXSim;
    DTheta_mnsMultiXSim=c*DQ_mnsMultiXSim;

%     AlphaMulti=zeros(4,Ma+1);
%     AlphaMultiXSim(:,1)=Q_plsMultiXSim(1:4);
%     AlphaMultiXSim(:,Ma+1)=Q_mnsMultiXSim(1:4);
%     AlphaMultiXSim(:,2)=DQ_plsMultiXSim(1:4)/(Ma*DTheta_plsMultiXSim)*(Theta_mnsMultiXSim-Theta_plsMultiXSim)+AlphaMultiXSim(:,1);
%     AlphaMultiXSim(:,Ma)=-DQ_mnsMultiXSim(1:4)/(Ma*DTheta_mnsMultiXSim)*(Theta_mnsMultiXSim-Theta_plsMultiXSim)+AlphaMultiXSim(:,Ma+1);
%     
%     for i=2+1:Ma-2+1
%         AlphaMulti(:,i)=(0*AlphaMulti(:,2)+50*AlphaMulti(:,Ma-1+1))/50;
%     end


    MotionDataMultiXSim=[MotionDataMultiXSim, ForceTorqueCalculator_5DoF([],[],'done')];

    q1=SolCMultiXSim(:,1)';q2=SolCMultiXSim(:,2)';q3=SolCMultiXSim(:,3)';q4=SolCMultiXSim(:,4)';q5=SolCMultiXSim(:,5)';
    dq1=SolCMultiXSim(:,6)';dq2=SolCMultiXSim(:,7)';dq3=SolCMultiXSim(:,8)';dq4=SolCMultiXSim(:,9)';dq5=SolCMultiXSim(:,10)';

    p_tib2_MulitXSim =[p_tib2_MulitXSim,  [ + L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5),
                                    - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + L_tib*cos(q2 + q4 + q5)]];


    v_hip_MulitXSim =[v_hip_MulitXSim , [dq1.*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) - dq2.*(L_fem*cos(q2 + q5) + L_tib*cos(q2 + q4 + q5)) + dq5.*(L_fem*cos(q1 + q5) - L_fem*cos(q2 + q5) + L_tib*cos(q1 + q3 + q5) - L_tib*cos(q2 + q4 + q5)) + L_tib*dq3.*cos(q1 + q3 + q5) - L_tib*dq4.*cos(q2 + q4 + q5)
                                         dq1.*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) - dq2.*(L_fem*sin(q2 + q5) + L_tib*sin(q2 + q4 + q5)) + dq5.*(L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5)) + L_tib*dq3.*sin(q1 + q3 + q5) - L_tib*dq4.*sin(q2 + q4 + q5)]];
        

    QqMultiXSim =[QqMultiXSim ;SolCMultiXSim(:,1:5)];
    DQqMultiXSim=[DQqMultiXSim ;SolCMultiXSim(:,6:10)];
    TimeMultiXSim=[TimeMultiXSim; TimeMultiXSim(end)+TcMultiXSim];
    T_impactMultiXSim=[T_impactMultiXSim TimeMultiXSim(end)];
        
end


% stability
qzM=@(xi)Hi*[BezierFunction_hd_Fast(Ma,xi,AlphaMultiXSim,Theta_plsMultiXSim,Theta_mnsMultiXSim);xi];
gamma0M=@(qz)...
        [ XX_fem + XX_tib + 2*L_fem^2*M_fem + L_fem^2*M_tib + 2*L_tib^2*M_fem + L_fem^2*M_torso + 2*L_tib^2*M_tib + L_tib^2*M_torso + Lc_fem^2*M_fem + Lc_tib^2*M_tib - 2*L_fem*Lc_fem*M_fem - 2*L_tib*Lc_tib*M_tib - L_fem^2*M_tib*cos(qz(1) - qz(2)) + 4*L_fem*L_tib*M_fem*cos(qz(3)) + 2*L_fem*L_tib*M_tib*cos(qz(3)) + 2*L_fem*L_tib*M_torso*cos(qz(3)) - L_fem*L_torso*M_torso*cos(qz(1)) - 2*L_tib*Lc_fem*M_fem*cos(qz(3)) + L_fem*Lc_torso*M_torso*cos(qz(1)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3)), XX_fem + XX_tib + L_fem^2*M_tib + Lc_fem^2*M_fem + Lc_tib^2*M_tib - L_fem^2*M_tib*cos(qz(1) - qz(2)) + 2*L_fem*Lc_tib*M_tib*cos(qz(4)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)), XX_tib + 2*L_tib^2*M_fem + 2*L_tib^2*M_tib + L_tib^2*M_torso + Lc_tib^2*M_tib - 2*L_tib*Lc_tib*M_tib + 2*L_fem*L_tib*M_fem*cos(qz(3)) + L_fem*L_tib*M_tib*cos(qz(3)) + L_fem*L_tib*M_torso*cos(qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(3)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3)), XX_tib + Lc_tib^2*M_tib + L_fem*Lc_tib*M_tib*cos(qz(4)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)), 2*XX_fem + 2*XX_tib + XX_torso + 2*L_fem^2*M_fem + 2*L_fem^2*M_tib + 2*L_tib^2*M_fem + L_fem^2*M_torso + 2*L_tib^2*M_tib + L_tib^2*M_torso + L_torso^2*M_torso + 2*Lc_fem^2*M_fem + 2*Lc_tib^2*M_tib + Lc_torso^2*M_torso - 2*L_fem*Lc_fem*M_fem - 2*L_tib*Lc_tib*M_tib - 2*L_torso*Lc_torso*M_torso - 2*L_fem^2*M_tib*cos(qz(1) - qz(2)) + 4*L_fem*L_tib*M_fem*cos(qz(3)) + 2*L_fem*L_tib*M_tib*cos(qz(3)) + 2*L_fem*L_tib*M_torso*cos(qz(3)) - 2*L_fem*L_torso*M_torso*cos(qz(1)) - 2*L_tib*Lc_fem*M_fem*cos(qz(3)) + 2*L_fem*Lc_tib*M_tib*cos(qz(4)) + 2*L_fem*Lc_torso*M_torso*cos(qz(1)) - 2*L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - 2*L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - 2*L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - 2*L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - 2*L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - 2*L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + 2*L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3))];
kappa1M=@(xi)c*inv([H0-BezierFunction_Dhd_Fast(Ma,xi,AlphaMultiXSim,Theta_plsMultiXSim,Theta_mnsMultiXSim)*c/(Theta_mnsMultiXSim-Theta_plsMultiXSim);gamma0M(qz(xi))])*[zeros(4,1);1];
VzeroM=0;
VzeroMaxM=-inf;
Dxi=-(Theta_plsMultiXSim-Theta_mnsMultiXSim)/50;
for Xi=linspace(Theta_plsMultiXSim,Theta_mnsMultiXSim,50)
    VzeroM=VzeroM-kappa2(qzM(Xi))/kappa1M(Xi);
    if(VzeroM>VzeroMaxM)
        VzeroMaxM=VzeroM;
    end
end
VzeroM=VzeroM*Dxi;
VzeroMaxM=VzeroMaxM*Dxi;
    
DeltaDqMinusM=@(VarQ_mnsMulti)DeltaDQ(VarQ_mnsMulti,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
LandaHatDqM=@(Var_qm)inv([H0-Ma*(AlphaMultiXSim(:,Ma+1)-AlphaMultiXSim(:,Ma))*c/(Theta_mnsMultiXSim-Theta_plsMultiXSim);gamma0M(Var_qm)])*[zeros(4,1);1];
deltaZeroM=gamma0M(Q_plsMultiXSim)*DeltaDqMinusM(Q_mnsMultiXSim)*LandaHatDqM(Q_mnsMultiXSim);

disp(['V_zero= ',num2str(VzeroM),',  V_zMax= ',num2str(VzeroMaxM),',  delta_zero=0< ',num2str(deltaZeroM),' <1'])
disp(['(deltaZ^2)/(1-deltaZ^2)*Vz+VzMax <0  :',num2str((deltaZeroM^2)/(1-deltaZeroM^2)*VzeroM+VzeroMaxM)])
disp(['Vel_hip_x_avg=',num2str(sum(v_hip_MulitXSim(1,end-length(TcMultiXSim)+2:end).*diff(TcMultiXSim'))/TcMultiXSim(end)),'m/s,    Pos_tib2_y_max=',num2str(max(p_tib2_MulitXSim(2,:))*100),'cm'])
disp(['ang1=',num2str(angl1Sim),', ang2=',num2str(angl2Sim)])
disp('---------------------')

TimeMultiXSim(1)=[];

ShowTimeMultiStep(TimeMultiXSim,[0 T_impactMultiXSim],QqMultiXSim,DQqMultiXSim,MotionDataMultiXSim,p_tib2_MulitXSim,v_hip_MulitXSim,mu)
AnimBot3DOF( TimeMultiXSim, QqMultiXSim,T_impactMultiXSim,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,0)

%% Optimization damper-spring on Multi Ode on Sim GA

ii=length(T_impactMultiXSim)/2;
indx1=find(TimeMultiXSim==T_impactMultiXSim(2*(ii)-2),1)+1;
indx2=find(TimeMultiXSim==T_impactMultiXSim(2*(ii)-1),1);
indx3=find(TimeMultiXSim==T_impactMultiXSim(2*(ii)-0),1);

TimeStrideMultiXSim=TimeMultiXSim(indx1:indx3);

Ur1MultiXSim=MotionDataMultiXSim(1,indx1:indx2)';
Ur1MultiXSim=[Ur1MultiXSim;MotionDataMultiXSim(2,indx2+1:indx3)'];
Ur3MultiXSim=MotionDataMultiXSim(3,indx1:indx2)';
Ur3MultiXSim=[Ur3MultiXSim;MotionDataMultiXSim(4,indx2+1:indx3)'];

Q1MultiXSim=QqMultiXSim(indx1:indx2,1);
Q1MultiXSim=[Q1MultiXSim;QqMultiXSim(indx2+1:indx3,2)];
Q3MultiXSim=QqMultiXSim(indx1:indx2,3);
Q3MultiXSim=[Q3MultiXSim;QqMultiXSim(indx2+1:indx3,4)];

Qhat13MultiXSim=Q1MultiXSim+Q3MultiXSim;

DtQ1MultiXSim=DQqMultiXSim(indx1:indx2,1);
DtQ1MultiXSim=[DtQ1MultiXSim;DQqMultiXSim(indx2+1:indx3,2)];
DtQ3MultiXSim=DQqMultiXSim(indx1:indx2  ,3);
DtQ3MultiXSim=[DtQ3MultiXSim;DQqMultiXSim(indx2+1:indx3,4)];

[CostActuationMultiXSim,CostD2QMultiXSim,CostParaRegMultiXSim,BetaOptimalMultiXSim,ThetaOptimalMultiXSim,EtaOptimalMultiXSim,TorqueMonoSpring1MultiXSim,TorqueMonoSpring3MultiXSim,TorqueDamper1MultiXSim,TorqueDamper3MultiXSim,TorqueBiSpring13MultiXSim]=...
                    OptimalParam(TimeStrideMultiXSim,Q1MultiXSim,Q3MultiXSim,Qhat13MultiXSim,DtQ1MultiXSim,DtQ3MultiXSim,Ur1MultiXSim,Ur3MultiXSim,rM,rB,rD,Landa,Gamma,Weight,SampleRate,1,1);
                
Ua1MultiXSim=(Ur1MultiXSim -TorqueMonoSpring1MultiXSim -TorqueDamper1MultiXSim -TorqueBiSpring13MultiXSim );
Ua3MultiXSim=(Ur3MultiXSim -TorqueMonoSpring3MultiXSim -TorqueDamper3MultiXSim -TorqueBiSpring13MultiXSim );


Work_r_MultiXSim=sum(abs(Ur1MultiXSim.*DtQ1MultiXSim).*[diff(TimeStrideMultiXSim) ;0])+sum(abs(Ur3MultiXSim.*DtQ3MultiXSim).*[diff(TimeStrideMultiXSim) ;0]);
Work_a_MultiXSim=sum(abs(Ua1MultiXSim.*DtQ1MultiXSim).*[diff(TimeStrideMultiXSim) ;0])+sum(abs(Ua3MultiXSim.*DtQ3MultiXSim).*[diff(TimeStrideMultiXSim) ;0]);

CostReq_MultiXSim=sum((Ur1MultiXSim.*Ur1MultiXSim+Ur3MultiXSim.*Ur3MultiXSim).*[diff(TimeStrideMultiXSim);0])/2;

Title=sprintf('%22s  % 12s % 12s'  ,    'Cost','Req. Work','Act. Work');
Result_Init_MultiXSim=sprintf('%-11s %11.2f %10.2f %10.2f',   'Initial:',CostReq_MultiXSim, Work_r_MultiXSim,Work_r_MultiXSim);
Result_Opt_MultiXSim =sprintf('%-11s %11.2f %10.2f %10.2f \n',   'Optimized:',CostActuationMultiXSim, Work_r_MultiXSim,Work_a_MultiXSim);

disp(Title)
disp(Result_Init_MultiXSim)
disp(Result_Opt_MultiXSim)
