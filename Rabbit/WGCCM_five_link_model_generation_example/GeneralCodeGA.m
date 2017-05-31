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
Q_minus=TT*deg2rad([180-angl1+angl2 , 180+angl1+angl2 , 180-angl1+angl2-2*angl2 , 180+angl1+angl2-2*angl2 , qTorso])'; % This Q_minus
%  Q_minus=deg2rad([ 183.5  215.3  -20.1  -19.0   -9.5])';
% Q_minus=TT*deg2rad([240 , 120 , 120 , 90 , 0])' % This Q_minus
% Q_minus=deg2rad([183.0600  206.9524  -13.1551    -9.5970   -8.9210])';% extracted from a converged walking
% DQ_minus=10*TT*deg2rad([-10 10 -70 70 -1])'; % This Dq_minus
DQ_minus=50*deg2rad([-.1  .1  -2 1 -1 ])'; % This Dq_minus --> stable for 10 step
DQ_minus=1.5*deg2rad([ 5.5   -1.9  -65.0   27.8  -33.5])';
DQ_minus=deg2rad([ 50.82   -27.49  -91.8   16.57  -56.44])';% extracted from a converged walking
  % DQ_minus=50*deg2rad([ .1 -.2 -1 -1 -.5 ])'; 
% Impact and Relabaling

[ Q_plus,DQ_plus,V_tib2_plus,F2]=ImpactModel(Q_minus,1*DQ_minus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
disp(['Q_plus=[ ',num2str(Q_plus'),']'])
disp(['DQ_plus=[ ',num2str(DQ_plus'),']'])
disp(['V_tib2_plus=[ ',num2str(V_tib2_plus'),']']);
disp(['F2=[ ',num2str(F2'),']']);
disp(['[F_n*mu,Ft]=[ ',num2str(F2(2)*mu),'>',num2str(F2(1)),']']);
disp('---------------------')

[KE_minus,PE_minus ] = KineticPotentialEnergy(Q_minus,DQ_minus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[KE_plus,PE_plus ] = KineticPotentialEnergy(Q_plus,DQ_plus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[p_tib2_minus ] = FoottPositon(Q_minus,L_fem, L_tib)';
[p_tib2_plus ] = FoottPositon(Q_plus,L_fem, L_tib)';
%% Bezier Coeffients

Ma=4;
c=[-1 0 -1/2 0 -1];
Theta_plus=c*Q_plus;
Theta_minus=c*Q_minus;

DTheta_plus=c*DQ_plus;
DTheta_minus=c*DQ_minus;
disp(['D_Theta_plus= ',num2str(DTheta_plus')])
disp(['D_Theta_minus= ',num2str(DTheta_minus')])
disp('---------------------')


Alfa=zeros(4,Ma+1);
Alfa(:,1)=Q_plus(1:4);
Alfa(:,Ma+1)=Q_minus(1:4);
Alfa(:,2)=DQ_plus(1:4)/(Ma*DTheta_plus)*(Theta_minus-Theta_plus)+Alfa(:,1);
Alfa(:,Ma)=-DQ_minus(1:4)/(Ma*DTheta_minus)*(Theta_minus-Theta_plus)+Alfa(:,Ma+1);

for i=2+1:Ma-2+1
    Alfa(:,i)=(1*Alfa(:,2)+1*Alfa(:,Ma-1+1))/2;
end


Hdr=zeros(4,1);
DHdr_s=zeros(4,1);
for Theta=linspace(Theta_plus,Theta_minus,50)
    [Hdr(:,end+1),DHdr_s(:,end+1)]=BezierFunction(Ma,Theta,Alfa,Theta_plus,Theta_minus);
end
Hdr(:,1)=[];
DHdr_s(:,1)=[];
SS=((linspace(Theta_plus,Theta_minus,50))-Theta_plus)/(Theta_minus-Theta_plus);
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
Kp=30000;
Kd=350;

Options = odeset('RelTol',2e-3,'AbsTol',2e-3,'maxstep',2e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Ma,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));
            
[Tc,SolC] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Ma,Kp,Kd),...
                    [0,2],[Q_plus; DQ_plus],Options);

Q_mns =SolC(end,1:5 )';
DQ_mns=SolC(end,6:10)';

[Q_pls,DQ_pls,V_tib2_pls,F2]=ImpactModel(Q_mns,DQ_mns,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
disp(['Q_plus=[ ',num2str(Q_pls'),'] (designed)'])
disp(['Q_pls=[ ',num2str(Q_pls'),'] (real)'])
disp(['DQ_plus=[ ',num2str(DQ_plus'),'] (designed)'])
disp(['DQ_pls=[ ',num2str(DQ_pls'),'] (real)'])
disp(['V_tib2_plus=[ ',num2str(V_tib2_pls'),']']);
                
MotionData=ForceTorqueCalculator_5DoF([],[],'done');
        
T_impact=Tc(end); 

[KE_mns,PE_mns ] = KineticPotentialEnergy(Q_mns,DQ_minus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[KE_pls,PE_pls ] = KineticPotentialEnergy(Q_plus,DQ_plus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
disp(['K_mns= ',num2str(KE_mns'),',  K_pls= ',num2str(KE_pls')])
disp(['CostAct=',num2str(sum((MotionData(1,:).*MotionData(1,:)+MotionData(3,:).*MotionData(3,:)+...
                                    MotionData(2,:).*MotionData(2,:)+MotionData(4,:).*MotionData(4,:))'.*[diff(Tc) ;0]/2))])
disp('---------------------')
close all


% 
AnimBot3DOF(Tc  ,SolC(:,1:5),T_impact,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,20)
ShowTimeOneStep(Tc,SolC,MotionData,c,Theta_plus,Theta_minus,Alfa,L_fem, L_tib, L_torso,mu,Ma);

%% GA 

% Kp=30000;
% Kd=350;
% angl1=10;% deg
% angl2=6;% deg
% qTorso=-10;%deg

Initial=[reshape(Alfa(:,2+1:Ma-2+1),(Ma-3)*4,1); angl1;angl2;qTorso; DQ_minus];


% WeightMatrix
Weight=[ 10 1 1]';

PopulationSize=30;
PopInitRange=[-10*ones(1,(Ma-3)*4)  5  0  -45  -100*ones(1,5);
               10*ones(1,(Ma-3)*4) 45  45   0     100*ones(1,5)];

GenerationsLimit=100;
StallGenLimit=10;
nvars=(Ma-3)*4+3+5;
tic

InitialPopulation=repmat(Initial',PopulationSize,1)+0.5*randn(PopulationSize,nvars);

CostFun   = @(Param)CostFunctionGA(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Kp,Kd);
NonConstr = @(Param)ConstraintFunctionGA(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g);

[Xx,fval,exitflag,output,population,score] =...
    Op_GA(CostFun, NonConstr, nvars, InitialPopulation, PopInitRange, PopulationSize,GenerationsLimit,StallGenLimit);

%% Ode on GA
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

% for i=2+1:Ma-2+1
%     AlfaX(:,i)=(0*AlfaX(:,2)+50*AlfaX(:,Ma))/50;
% end

Options = odeset('RelTol',1e-3,'AbsTol',1e-3,'maxstep',1e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlfaX,Theta_plusX,Theta_minusX,Ma,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));
            
[TcX,SolCX] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlfaX,Theta_plusX,Theta_minusX,Ma,Kp,Kd),...
                    [0,2],[Q_plusX; DQ_plusX],Options);

Q_mnsX =SolCX(end,1:5 )';
DQ_mnsX=SolCX(end,6:10)';
[Q_plsX,DQ_plsX,V_tib2_plsX,F2X]=ImpactModel(Q_mnsX,DQ_mnsX,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);

Theta_plsX=c*Q_plsX;
Theta_mnsX=c*Q_mnsX;
DTheta_plsX=c*DQ_plsX;
DTheta_mnsX=c*DQ_mnsX;

disp(['Q_pls=[ ',num2str(Q_plsX'),']'])
disp(['DQ_pls=[ ',num2str(DQ_plsX'),']'])
disp(['V_tib2_plus=[ ',num2str(V_tib2_plsX'),']']);
                
MotionDataX=ForceTorqueCalculator_5DoF([],[],'done');

T_impactX=[TcX(end)]; 

[KE_mnsX,PE_mnsX ] = KineticPotentialEnergy(Q_mnsX,DQ_mnsX,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[KE_plsX,PE_plsX ] = KineticPotentialEnergy(Q_plsX,DQ_plsX,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
disp(['K_mns= ',num2str(KE_mnsX'),',  K_pls= ',num2str(KE_plsX')])
disp(['CostAct=',num2str(sum((MotionDataX(1,:).*MotionDataX(1,:)+MotionDataX(3,:).*MotionDataX(3,:)+...
                                    MotionDataX(2,:).*MotionDataX(2,:)+MotionDataX(4,:).*MotionDataX(4,:))'.*[diff(TcX) ;0]/2))])
disp('---------------------')
close all

% stabolity
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
LandaHatDq=@(qm)inv([H0-Ma*(AlfaX(:,Ma+1)-Alfa(:,Ma))*c/(Theta_minusX-Theta_plusX);gamma0(qm)])*[zeros(4,1);1];
deltaZero=gamma0(Q_plusX)*DeltaDqMinus(Q_minusX)*LandaHatDq(Q_minusX);

disp(['V_zero=    ',num2str(Vzero)])
disp(['V_zeroMax= ',num2str(VzeroMax)])
disp(['delta_zero=0<  ',num2str(deltaZero),'  <1'])
disp(['(deltaZ^2)/(1-deltaZ^2)*Vz+VzMax <0  :',num2str((deltaZero^2)/(1-deltaZero^2)*Vzero+VzeroMax)])
disp('---------------------')


% 
AnimBot3DOF(TcX  ,SolCX(:,1:5),[ T_impactX],L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,2)
ShowTimeOneStep(TcX,SolCX,MotionDataX,c,Theta_plusX,Theta_mnsX,AlfaX,L_fem, L_tib, L_torso,mu,Ma);

%% Multi Ode on GA

% Kp=4000;
% Kd=130;


AlphaMulti=AlfaX;

Q_mnsMulti=Q_minusX;
DQ_mnsMulti=DQ_minusX;
Q_plsMulti=Q_plusX;
DQ_plsMulti=DQ_plusX;
Theta_mnsMulti= Theta_minusX;
Theta_plsMulti= Theta_plusX;

T_impactMulti=[]; 
QqMulti=[];
DQqMulti=[];
TimeMutli=[0];
MotionDataMulti=[];
p_tib2_Mulit =[];
v_hip_Mulit =[];
    
for step=1:2*60
    step
    
    
    Options = odeset('RelTol',1e-3,'AbsTol',1e-3,'maxstep',1e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlphaMulti,Theta_plsMulti,Theta_mnsMulti,Ma,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));

    [TcMulti,SolCMulti] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlphaMulti,Theta_plsMulti,Theta_mnsMulti,Ma,Kp,Kd),...
                        [0,2],[Q_plsMulti; DQ_plsMulti],Options);

    Q_mnsMulti =SolCMulti(end,1:5 )';
    DQ_mnsMulti=SolCMulti(end,6:10)';

    [ Q_plsMulti,DQ_plsMulti,V_tib2_pls_Mulit,F2_Mulit]=ImpactModel(Q_mnsMulti,DQ_mnsMulti,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
%     disp(['Q_plus=[ ',num2str(Q_pls'),']']);
    disp(['DQ_plus=[ ',num2str(DQ_plsMulti'),']']);
%     disp(['V_tib2_plus=[ ',num2str(V_tib2_pls'),']']);
%                 
    Theta_plsMulti=c*Q_plsMulti;
    Theta_mnsMulti=c*Q_mnsMulti;

    DTheta_plsMulti=c*DQ_plsMulti;
    DTheta_mnsMulti=c*DQ_mnsMulti;

%     AlphaMulti=zeros(4,Ma+1);
    AlphaMulti(:,1)=Q_plsMulti(1:4);
    AlphaMulti(:,Ma+1)=Q_mnsMulti(1:4);
    AlphaMulti(:,2)=DQ_plsMulti(1:4)/(Ma*DTheta_plsMulti)*(Theta_mnsMulti-Theta_plsMulti)+AlphaMulti(:,1);
    AlphaMulti(:,Ma)=-DQ_mnsMulti(1:4)/(Ma*DTheta_mnsMulti)*(Theta_mnsMulti-Theta_plsMulti)+AlphaMulti(:,Ma+1);
    
%     for i=2+1:Ma-2+1
%         AlphaMulti(:,i)=(0*AlphaMulti(:,2)+50*AlphaMulti(:,Ma-1+1))/50;
%     end


    MotionDataMulti=[MotionDataMulti, ForceTorqueCalculator_5DoF([],[],'done')];

    q1=SolCMulti(:,1)';
    q2=SolCMulti(:,2)';
    q3=SolCMulti(:,3)';
    q4=SolCMulti(:,4)';
    q5=SolCMulti(:,5)';
    dq1=SolCMulti(:,6)';
    dq2=SolCMulti(:,7)';
    dq3=SolCMulti(:,8)';
    dq4=SolCMulti(:,9)';
    dq5=SolCMulti(:,10)';

    p_tib2_Mulit =[p_tib2_Mulit,  [ + L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5),
                                    - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + L_tib*cos(q2 + q4 + q5)]];


    v_hip_Mulit =[v_hip_Mulit , [dq1.*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) - dq2.*(L_fem*cos(q2 + q5) + L_tib*cos(q2 + q4 + q5)) + dq5.*(L_fem*cos(q1 + q5) - L_fem*cos(q2 + q5) + L_tib*cos(q1 + q3 + q5) - L_tib*cos(q2 + q4 + q5)) + L_tib*dq3.*cos(q1 + q3 + q5) - L_tib*dq4.*cos(q2 + q4 + q5)
                                 dq1.*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) - dq2.*(L_fem*sin(q2 + q5) + L_tib*sin(q2 + q4 + q5)) + dq5.*(L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5)) + L_tib*dq3.*sin(q1 + q3 + q5) - L_tib*dq4.*sin(q2 + q4 + q5)]];
        

    QqMulti =[QqMulti ;SolCMulti(:,1:5)];
    DQqMulti=[DQqMulti ;SolCMulti(:,6:10)];
    TimeMutli=[TimeMutli; TimeMutli(end)+TcMulti];
    T_impactMulti=[T_impactMulti TimeMutli(end)];
    
    
end


% stabolity
qzM=@(xi)Hi*[BezierFunction_hd_Fast(Ma,xi,AlphaMulti,Theta_plsMulti,Theta_mnsMulti);xi];
gamma0M=@(qz)...
        [ XX_fem + XX_tib + 2*L_fem^2*M_fem + L_fem^2*M_tib + 2*L_tib^2*M_fem + L_fem^2*M_torso + 2*L_tib^2*M_tib + L_tib^2*M_torso + Lc_fem^2*M_fem + Lc_tib^2*M_tib - 2*L_fem*Lc_fem*M_fem - 2*L_tib*Lc_tib*M_tib - L_fem^2*M_tib*cos(qz(1) - qz(2)) + 4*L_fem*L_tib*M_fem*cos(qz(3)) + 2*L_fem*L_tib*M_tib*cos(qz(3)) + 2*L_fem*L_tib*M_torso*cos(qz(3)) - L_fem*L_torso*M_torso*cos(qz(1)) - 2*L_tib*Lc_fem*M_fem*cos(qz(3)) + L_fem*Lc_torso*M_torso*cos(qz(1)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3)), XX_fem + XX_tib + L_fem^2*M_tib + Lc_fem^2*M_fem + Lc_tib^2*M_tib - L_fem^2*M_tib*cos(qz(1) - qz(2)) + 2*L_fem*Lc_tib*M_tib*cos(qz(4)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)), XX_tib + 2*L_tib^2*M_fem + 2*L_tib^2*M_tib + L_tib^2*M_torso + Lc_tib^2*M_tib - 2*L_tib*Lc_tib*M_tib + 2*L_fem*L_tib*M_fem*cos(qz(3)) + L_fem*L_tib*M_tib*cos(qz(3)) + L_fem*L_tib*M_torso*cos(qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(3)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3)), XX_tib + Lc_tib^2*M_tib + L_fem*Lc_tib*M_tib*cos(qz(4)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)), 2*XX_fem + 2*XX_tib + XX_torso + 2*L_fem^2*M_fem + 2*L_fem^2*M_tib + 2*L_tib^2*M_fem + L_fem^2*M_torso + 2*L_tib^2*M_tib + L_tib^2*M_torso + L_torso^2*M_torso + 2*Lc_fem^2*M_fem + 2*Lc_tib^2*M_tib + Lc_torso^2*M_torso - 2*L_fem*Lc_fem*M_fem - 2*L_tib*Lc_tib*M_tib - 2*L_torso*Lc_torso*M_torso - 2*L_fem^2*M_tib*cos(qz(1) - qz(2)) + 4*L_fem*L_tib*M_fem*cos(qz(3)) + 2*L_fem*L_tib*M_tib*cos(qz(3)) + 2*L_fem*L_tib*M_torso*cos(qz(3)) - 2*L_fem*L_torso*M_torso*cos(qz(1)) - 2*L_tib*Lc_fem*M_fem*cos(qz(3)) + 2*L_fem*Lc_tib*M_tib*cos(qz(4)) + 2*L_fem*Lc_torso*M_torso*cos(qz(1)) - 2*L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - 2*L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - 2*L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - 2*L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - 2*L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - 2*L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + 2*L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3))];
kappa1M=@(xi)c*inv([H0-BezierFunction_Dhd_Fast(Ma,xi,AlphaMulti,Theta_plsMulti,Theta_mnsMulti)*c/(Theta_mnsMulti-Theta_plsMulti);gamma0M(qz(xi))])*[zeros(4,1);1];
VzeroM=0;
VzeroMaxM=-inf;
Dxi=-(Theta_plsMulti-Theta_mnsMulti)/50;
for Xi=linspace(Theta_plsMulti,Theta_mnsMulti,50)
    VzeroM=VzeroM-kappa2(qzM(Xi))/kappa1M(Xi);
    if(VzeroM>VzeroMaxM)
        VzeroMaxM=VzeroM;
    end
end
VzeroM=VzeroM*Dxi;
VzeroMaxM=VzeroMaxM*Dxi;
    
DeltaDqMinusM=@(Q_mnsMulti)DeltaDQ(Q_mnsMulti,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
LandaHatDqM=@(qm)inv([H0-Ma*(AlphaMulti(:,Ma+1)-AlphaMulti(:,Ma))*c/(Theta_mnsMulti-Theta_plsMulti);gamma0M(qm)])*[zeros(4,1);1];
deltaZeroM=gamma0M(Q_plsMulti)*DeltaDqMinusM(Q_mnsMulti)*LandaHatDqM(Q_mnsMulti);

disp(['V_zero=    ',num2str(VzeroM)])
disp(['V_zeroMax= ',num2str(VzeroMaxM)])
disp(['delta_zero=0<  ',num2str(deltaZeroM),'  <1'])
disp(['(deltaZ^2)/(1-deltaZ^2)*Vz+VzMax <0  :',num2str((deltaZeroM^2)/(1-deltaZeroM^2)*VzeroM+VzeroMaxM)])
disp('---------------------')



TimeMutli(1)=[];

AnimBot3DOF( TimeMutli, QqMulti,T_impactMulti,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,2)
ShowTimeMultiStep(TimeMutli,T_impactMulti,QqMulti,DQqMulti,MotionDataMulti,p_tib2_Mulit,v_hip_Mulit,mu)

%% %% Optimization damper-spring on GA
Weight=[1 1];
SampleRate=10;   
rM=4-0;
rB=4-0;
rD=0+0-1;
Landa=0.00000005;
Gamma=0.0005;

i=length(T_impactMulti)/2;
indx1=find(TimeMutli==T_impactMulti(2*(i)-2),1)+1;
indx2=find(TimeMutli==T_impactMulti(2*(i)-1),1);
indx3=find(TimeMutli==T_impactMulti(2*(i)-0),1);

TimeStride=TimeMutli(indx1:indx3);

Ur1=MotionDataMulti(1,indx1:indx2)';
Ur1=[Ur1;MotionDataMulti(2,indx2+1:indx3)'];
Ur3=MotionDataMulti(3,indx1:indx2)';
Ur3=[Ur3;MotionDataMulti(4,indx2+1:indx3)'];

Q1=QqMulti(indx1:indx2,1);
Q1=[Q1;QqMulti(indx2+1:indx3,2)];
Q3=QqMulti(indx1:indx2,3);
Q3=[Q3;QqMulti(indx2+1:indx3,4)];

Qhat13=Q1+Q3;

DtQ1=DQqMulti(indx1:indx2,1);
DtQ1=[DtQ1;DQqMulti(indx2+1:indx3,2)];
DtQ3=DQqMulti(indx1:indx2  ,3);
DtQ3=[DtQ3;DQqMulti(indx2+1:indx3,4)];

[CostActuation,CostD2Q,CostParaReg,BetaOptimal,ThetaOptimal,EtaOptimal,TorqueMonoSpring1,TorqueMonoSpring3,TorqueDamper1,TorqueDamper3,TorqueBiSpring13]=...
                    OptimalParam(TimeStride,Q1,Q3,Qhat13,DtQ1,DtQ3,Ur1,Ur3,rM,rB,rD,Landa,Gamma,Weight,SampleRate,1,0);
                
Ua1=(Ur1 -TorqueMonoSpring1 -TorqueDamper1 -TorqueBiSpring13 );
Ua3=(Ur3 -TorqueMonoSpring3 -TorqueDamper3 -TorqueBiSpring13 );


Work_r=sum(abs(Ur1.*DtQ1).*[diff(TimeStride) ;0])+sum(abs(Ur3.*DtQ3).*[diff(TimeStride) ;0]);
Work_a=sum(abs(Ua1.*DtQ1).*[diff(TimeStride) ;0])+sum(abs(Ua3.*DtQ3).*[diff(TimeStride) ;0]);

CostReq=sum((Ur1.*Ur1+Ur3.*Ur3).*[diff(TimeStride);0])/2;

Title=sprintf('%22s  % 12s % 12s'  ,    'Cost','Req. Work','Act. Work');
Result_Init=sprintf('%-11s %11.2f %10.2f %10.2f',   'Initial:',CostReq, Work_r,Work_r);
Result_Opt =sprintf('%-11s %11.2f %10.2f %10.2f \n',   'Optimized:',CostActuation, Work_r,Work_a);

disp(Title)
disp(Result_Init)
disp(Result_Opt)

%% Simultinous GA
% Kp=30000;
% Kd=350;
% angl1=10;% deg
% angl2=6;% deg
% qTorso=-10;%deg

Initial=[reshape(Alfa(:,2+1:Ma-2+1),(Ma-3)*4,1); angl1;angl2;qTorso; DQ_minus];


PopulationSize=30;
PopInitRange=[-10*ones(1,(Ma-3)*4)  5  0  -45  -100*ones(1,5);
               10*ones(1,(Ma-3)*4) 45  45   0     100*ones(1,5)];

GenerationsLimit=100;
StallGenLimit=10;
nvars=(Ma-3)*4+3+5;
tic

InitialPopulation=repmat(Initial',PopulationSize,1)+0.5*randn(PopulationSize,nvars);

CostFunSim   = @(Param)CostFunctionSimGA(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Kp,Kd,rM,rB,rD,Landa,Gamma,Weight,SampleRate);
NonConstr = @(Param)ConstraintFunctionGA(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g);

[XxSim,fval,exitflag,output,population,score] =...
    Op_GA(CostFunSim, NonConstr, nvars, InitialPopulation, PopInitRange, PopulationSize,GenerationsLimit,StallGenLimit);

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

% for i=2+1:Ma-2+1
%     AlfaX(:,i)=(0*AlfaX(:,2)+50*AlfaX(:,Ma))/50;
% end

Options = odeset('RelTol',5e-3,'AbsTol',5e-3,'maxstep',5e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlfaXSim,Theta_plusXSim,Theta_minusXSim,Ma,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));
            
[TcXSim,SolCXSim] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlfaXSim,Theta_plusXSim,Theta_minusXSim,Ma,Kp,Kd),...
                    [0,2],[Q_plusXSim; DQ_plusXSim],Options);

Q_mnsXSim =SolCXSim(end,1:5 )';
DQ_mnsXSim=SolCXSim(end,6:10)';
[Q_plsXSim,DQ_plsXSim,V_tib2_plsXSim,F2XSim]=ImpactModel(Q_mnsXSim,DQ_mnsXSim,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);

Theta_plsXSim=c*Q_plsXSim;
Theta_mnsXSim=c*Q_mnsXSim;
DTheta_plsXSim=c*DQ_plsXSim;
DTheta_mnsXSim=c*DQ_mnsXSim;

disp(['Q_pls=[ ',num2str(Q_plsXSim'),']'])
disp(['DQ_pls=[ ',num2str(DQ_plsXSim'),']'])
disp(['V_tib2_plus=[ ',num2str(V_tib2_plsXSim'),']']);
                
MotionDataXSim=ForceTorqueCalculator_5DoF([],[],'done');

T_impactXSim=[TcXSim(end)]; 

[KE_mnsXSim,PE_mnsXSim ] = KineticPotentialEnergy(Q_mnsXSim,DQ_mnsXSim,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[KE_plsXSim,PE_plsXSim ] = KineticPotentialEnergy(Q_plsXSim,DQ_plsXSim,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
disp(['K_mns= ',num2str(KE_mnsXSim'),',  K_pls= ',num2str(KE_plsXSim')])
disp(['CostAct=',num2str(sum((MotionDataXSim(1,:).*MotionDataXSim(1,:)+MotionDataXSim(3,:).*MotionDataXSim(3,:)+...
                                    MotionDataXSim(2,:).*MotionDataXSim(2,:)+MotionDataXSim(4,:).*MotionDataXSim(4,:))'.*[diff(TcXSim) ;0]/2))])
disp('---------------------')
close all

% stabolity
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

disp(['V_zero=    ',num2str(Vzero)])
disp(['V_zeroMax= ',num2str(VzeroMax)])
disp(['delta_zero=0<  ',num2str(deltaZero),'  <1'])
disp(['(deltaZ^2)/(1-deltaZ^2)*Vz+VzMax <0  :',num2str((deltaZero^2)/(1-deltaZero^2)*Vzero+VzeroMax)])
disp('---------------------')


% 
AnimBot3DOF(TcXSim  ,SolCXSim(:,1:5),[ T_impactXSim],L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,2)
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
Result_Init=sprintf('%-11s %11.2f %10.2f %10.2f',   'Initial:',CostReq, Work_rXSim,Work_rXSim);
Result_Opt =sprintf('%-11s %11.2f %10.2f %10.2f \n',   'Optimized:',CostActuationXSim, Work_rXSim,Work_a);

disp(Title)
disp(Result_Init)
disp(Result_Opt)