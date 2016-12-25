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

%% Initialization
TT = [1   0   0   0  -1;
      0   1   0   0  -1;
     -1   0   1   0   0;
	  0  -1   0   1   0;
	  0   0   0   0   1];
c=[-1 0 -1/2 0 -1];

% Q_minus=TT*deg2rad([170 240 150 200 -10])'; % 
angl1=15;
angl2=10;
Q_minus=TT*deg2rad([180-angl1+angl2 , 180+angl1+angl2 , 180-angl1+angl2-2*angl2 , 180+angl1+angl2-2*angl2 , -10])'; % This Q_minus
Q_minus=deg2rad([ 183.5  215.3  -20.1  -19.0   -9.5])';
% Q_minus=TT*deg2rad([240 , 120 , 120 , 90 , 0])' % This Q_minus
% DQ_minus=10*TT*deg2rad([-10 10 -70 70 -1])'; % This Dq_minus
DQ_minus=50*deg2rad([-.1  .1  -2 1 -1 ])'; % This Dq_minus --> stable for 10 step
DQ_minus=1*deg2rad([ 5.5   -1.9  -65.0   27.8  -33.5])';
  % DQ_minus=50*deg2rad([ .1 -.2 -1 -1 -.5 ])'; 
% Impact and Relabaling

[ Q_plus,DQ_plus,V_tib2_plus,F2]=ImpactModel(Q_minus,DQ_minus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
disp(['Q_plus=[ ',num2str(Q_plus'),']'])
disp(['DQ_plus=[ ',num2str(DQ_plus'),']'])
disp(['V_tib2_plus=[ ',num2str(V_tib2_plus'),']']);
disp('---------------------')

[KE_minus,PE_minus ] = KineticPotentialEnergy(Q_minus,DQ_minus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[KE_plus,PE_plus ] = KineticPotentialEnergy(Q_plus,DQ_plus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[p_tib2_minus ] = FoottPositon(Q_minus,L_fem, L_tib)';
[p_tib2_plus ] = FoottPositon(Q_plus,L_fem, L_tib)';
%% Bezier Coeffients

Ma=3;
c=[-1 0 -1/2 0 -1];
Theta_plus=c*Q_plus;
Theta_minus=c*Q_minus;

DTheta_plus=c*DQ_plus;
DTheta_minus=c*DQ_minus;
disp(['D_Theta_plus= ',num2str(DTheta_plus')])
disp(['D_Theta_minus= ',num2str(DTheta_minus')])


Alfa=zeros(4,4);
Alfa(:,1)=Q_plus(1:4);
Alfa(:,4)=Q_minus(1:4);
Alfa(:,2)=DQ_plus(1:4)/(Ma*DTheta_plus)*(Theta_minus-Theta_plus)+Alfa(:,1);
Alfa(:,3)=-DQ_minus(1:4)/(Ma*DTheta_minus)*(Theta_minus-Theta_plus)+Alfa(:,4);



Hdr=zeros(4,1);
DHdr=zeros(4,1);
for Theta=linspace(Theta_plus,Theta_minus,50)
    [Hdr(:,end+1),DHdr(:,end+1)]=BezierFunction(Theta,Alfa,Theta_plus,Theta_minus);
end
Hdr(:,1)=[];
DHdr(:,1)=[];
SS=((linspace(Theta_plus,Theta_minus,50))-Theta_plus)/(Theta_minus-Theta_plus);
QQ=[Hdr;zeros(1,length(SS))];
DQQ=[DHdr;zeros(1,length(SS))];
figure(2)
plot(c*QQ)
hold all
plot(c*DQQ/(Theta_minus-Theta_plus))
grid on
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
%% Ode
Kp=30000;
Kd=350;

Options = odeset('RelTol',2e-3,'AbsTol',2e-3,'maxstep',2e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));
            
[Tc,SolC] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Kp,Kd),...
                    [0,2],[Q_plus; DQ_plus],Options);

Q_mns =SolC(end,1:5 )';
DQ_mns=SolC(end,6:10)';

[Q_pls,DQ_pls,V_tib2_pls,F2]=ImpactModel(Q_mns,DQ_mns,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
disp(['Q_pls=[ ',num2str(Q_pls'),']'])
disp(['DQ_pls=[ ',num2str(DQ_pls'),']'])
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
AnimBot3DOF(Tc  ,SolC(:,1:5),T_impact,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,2)
ShowTimeOneStep(Tc,SolC,MotionData,c,Theta_plus,Theta_minus,Alfa,L_fem, L_tib, L_torso);
%% Multi Ode 

% Kp=4000;
% Kd=130;


Q_pls=Q_plus;
DQ_pls=DQ_plus;
Alpha=Alfa;
Theta_mns= Theta_minus;
Theta_pls= Theta_plus;

T_impactMulti=[]; 
QqMulti=[];
DQqMulti=[];
TimeMutli=[0];
MotionDataMulti=[];
p_tib2 =[];
v_hip =[];
    
for step=1:2*8
    step
    
    
    Options = odeset('RelTol',1e-2,'AbsTol',1e-2,'maxstep',1e-2,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alpha,Theta_pls,Theta_mns,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));

    [TcMulti,SolCMulti] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alpha,Theta_pls,Theta_mns,Kp,Kd),...
                        [0,2],[Q_pls; DQ_pls],Options);

    Q_mns =SolCMulti(end,1:5 )';
    DQ_mns=SolCMulti(end,6:10)';

    [ Q_pls,DQ_pls,V_tib2_pls,F2]=ImpactModel(Q_mns,DQ_mns,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
%     disp(['Q_plus=[ ',num2str(Q_pls'),']']);
%     disp(['DQ_plus=[ ',num2str(DQ_pls'),']']);
%     disp(['V_tib2_plus=[ ',num2str(V_tib2_pls'),']']);
%                 
    Theta_pls=c*Q_pls;
    Theta_mns=c*Q_mns;

    DTheta_pls=c*DQ_pls;
    DTheta_mns=c*DQ_mns;

    Alpha=zeros(4,4);
    Alpha(:,1)=Q_pls(1:4);
    Alpha(:,4)=Q_mns(1:4);
    Alpha(:,2)=DQ_pls(1:4)/(Ma*DTheta_pls)*(Theta_mns-Theta_pls)+Alpha(:,1);
    Alpha(:,3)=-DQ_mns(1:4)/(Ma*DTheta_mns)*(Theta_mns-Theta_pls)+Alpha(:,4);

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

    p_tib2 =[p_tib2,  [ + L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5),
                      - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + L_tib*cos(q2 + q4 + q5)]];


    v_hip =[v_hip , [dq1.*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) - dq2.*(L_fem*cos(q2 + q5) + L_tib*cos(q2 + q4 + q5)) + dq5.*(L_fem*cos(q1 + q5) - L_fem*cos(q2 + q5) + L_tib*cos(q1 + q3 + q5) - L_tib*cos(q2 + q4 + q5)) + L_tib*dq3.*cos(q1 + q3 + q5) - L_tib*dq4.*cos(q2 + q4 + q5)
        dq1.*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) - dq2.*(L_fem*sin(q2 + q5) + L_tib*sin(q2 + q4 + q5)) + dq5.*(L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5)) + L_tib*dq3.*sin(q1 + q3 + q5) - L_tib*dq4.*sin(q2 + q4 + q5)]];
        

    QqMulti =[QqMulti ;SolCMulti(:,1:5)];
    DQqMulti=[DQqMulti ;SolCMulti(:,6:10)];
    TimeMutli=[TimeMutli; TimeMutli(end)+TcMulti];
    T_impactMulti=[T_impactMulti TimeMutli(end)];
    
    
end
TimeMutli(1)=[];

AnimBot3DOF( TimeMutli, QqMulti,T_impactMulti,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,1)
ShowTimeMultiStep(TimeMutli,T_impactMulti,QqMulti,DQqMulti,MotionDataMulti,p_tib2,v_hip)

%% Optimization damper-spring

Weight=[1 1];
SampleRate=10;   
rM=3-0;
rB=3-0;
rD=0+0-0;
Landa=0.0000;
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

[BetaOptimal,ThetaOptimal,EtaOptimal,CostActuation,CostD2Q,CostParaReg,TorqueMonoSpring1,TorqueMonoSpring3,TorqueDamper1,TorqueDamper3,TorqueBiSpring13]=...
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

%% Rabit With Spring-Damper

Options = odeset('RelTol',2e-3,'AbsTol',2e-3,'maxstep',2e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculatorWithSpringDamper_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Kp,Kd,rM,BetaOptimal,rD,EtaOptimal,rB,ThetaOptimal),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));

[TcSD,SolSD] = ode15s(@(t,x)RabitDynamicWithSpringDamper(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Kp,Kd,rM,BetaOptimal,rD,EtaOptimal,rB,ThetaOptimal),...
                    [0,2],[Q_plus; DQ_plus],Options);

Q_mnsSD =SolSD(end,1:5 )';
DQ_mnsSD=SolSD(end,6:10)';

[Q_plsSD,DQ_plsSD,V_tib2_plsSD,F2]=ImpactModel(Q_mnsSD,DQ_mnsSD,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
disp(['Q_pls=[ ',num2str(Q_plsSD'),']'])
disp(['DQ_pls=[ ',num2str(DQ_plsSD'),']'])
disp(['V_tib2_plus=[ ',num2str(V_tib2_plsSD'),']']);
                
MotionDataSD=ForceTorqueCalculatorWithSpringDamper_5DoF([],[],'done');

        
T_impact=TcSD(end); 

[KE_mns,PE_mns ] = KineticPotentialEnergy(Q_mnsSD,DQ_mnsSD,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[KE_pls,PE_pls ] = KineticPotentialEnergy(Q_plus,DQ_plus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
disp(['K_mns= ',num2str(KE_mns'),',  K_pls= ',num2str(KE_pls')])
disp(['CostAct=',num2str(sum((MotionDataSD(1,:).*MotionDataSD(1,:)+MotionDataSD(3,:).*MotionDataSD(3,:)+...
                                    MotionDataSD(2,:).*MotionDataSD(2,:)+MotionDataSD(4,:).*MotionDataSD(4,:))'.*[diff(TcSD) ;0]/2))])
disp('---------------------')

 
AnimBot3DOF(TcSD  ,SolSD(:,1:5),T_impact,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,2)
ShowTimeOneStep(TcSD,SolSD,MotionDataSD,c,Theta_plus,Theta_minus,Alfa,L_fem, L_tib, L_torso);
