
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
DQ_minus=deg2rad([ 5.5   -1.9  -65.0   27.8  -33.5])';
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
%% Bezier Coeffients Initialization
Ma=4;
c=[-1 0 -1/2 0 -1];
Theta_plus=c*Q_plus;
Theta_minus=c*Q_minus;

DTheta_plus=c*DQ_plus;
DTheta_minus=c*DQ_minus;
disp(['D_Theta_plus= ',num2str(DTheta_plus')])
disp(['D_Theta_minus= ',num2str(DTheta_minus')])


Alfa=zeros(4,Ma+1);
Alfa(:,1)=Q_plus(1:4);
Alfa(:,Ma+1)=Q_minus(1:4);
Alfa(:,2)=DQ_plus(1:4)/(Ma*DTheta_plus)*(Theta_minus-Theta_plus)+Alfa(:,1);
Alfa(:,Ma)=-DQ_minus(1:4)/(Ma*DTheta_minus)*(Theta_minus-Theta_plus)+Alfa(:,Ma+1);

for i=2+1:Ma-2+1
    Alfa(:,i)=(Alfa(:,2)+Alfa(:,Ma-1+1))/2;
end


Hdr=zeros(4,1);
DHdr=zeros(4,1);
for Theta=linspace(Theta_plus,Theta_minus,50)
    [Hdr(:,end+1),DHdr(:,end+1)]=BezierFunction(Theta,Alfa,Theta_plus,Theta_minus);
end
Hdr(:,1)=[];
DHdr(:,1)=[];
SS=((linspace(Theta_plus,Theta_minus,50))-Theta_plus)/(Theta_minus-Theta_plus);
QQ=[Hdr;ones(1,length(SS))*Q_plus(5)];
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

%% Optimization

Kp=30000;
Kd=350;

Initial=[reshape(Alfa(:,2+1:Ma-2+1),(Ma-3)*4,1); Q_minus; DQ_minus];
Initial=[reshape(Alfa(:,2+1:Ma-2+1),(Ma-3)*4,1)];

% WeightMatrix
Weight=[ 10 1 1]';

tic
MaxFunEvals_Data=5000;
MaxIter_Data=1000;
TolFun_Data=1e-5;
TolX_Data=1e-5;
TolCon_Data=1e-5;
Algorithm='sqp';
Algorithm='interior-point';

CostFun   = @(Param)CostFunction(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Kp,Kd);
NonConstr = @(Param)ConstraintFunction(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Kp,Kd);


[xx,fval,exitflag,output,lambda,grad,hessian] = ...
    Op_FmisCon_SQP(CostFun,NonConstr,Initial,MaxFunEvals_Data,MaxIter_Data,TolFun_Data,TolX_Data,TolCon_Data,Algorithm);




%% Ode

% Q_mnsX=xx(end-9:end-5);
Q_mnsX=deg2rad([ 183.5  215.3  -20.1  -19.0   -9.5])';
% DQ_mnsX=xx(end-4:end);
DQ_mnsX=deg2rad([ 5.5   -1.9  -65.0   27.8  -33.5])';


[Q_plsX,DQ_plsX,V_tib2_plus,F2]=ImpactModel(Q_mnsX,DQ_mnsX,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);


Theta_plsX=c*Q_plsX;
Theta_mnsX=c*Q_mnsX;
DTheta_plsX=c*DQ_plsX;
DTheta_mnsX=c*DQ_mnsX;


AlfaX=zeros(4,Ma+1);

AlfaX(:,0+1)=Q_plsX(1:4);
AlfaX(:,Ma+1)=Q_mnsX(1:4);
AlfaX(:,1+1)=DQ_plsX(1:4)/(Ma*DTheta_plsX)*(Theta_mnsX-Theta_plsX)+AlfaX(:,0+1);
AlfaX(:,Ma-1+1)=-DQ_mnsX(1:4)/(Ma*DTheta_mnsX)*(Theta_mnsX-Theta_plsX)+AlfaX(:,Ma+1);

Alfa(:,2+1:Ma-2+1)=reshape(xx(1:(Ma-3)*4),4,Ma-3);


Options = odeset('RelTol',2e-3,'AbsTol',2e-3,'maxstep',2e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlfaX,Theta_plsX,Theta_mnsX,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));
            
[Tc,SolC] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,AlfaX,Theta_plsX,Theta_mnsX,Kp,Kd),...
                    [0,2],[Q_plsX; DQ_plsX],Options);

Q_mns =SolC(end,1:5 )';
DQ_mns=SolC(end,6:10)';

[Q_pls,DQ_pls,V_tib2_pls,F2]=ImpactModel(Q_mns,DQ_mns,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
disp(['Q_pls=[ ',num2str(Q_pls'),']'])
disp(['DQ_pls=[ ',num2str(DQ_pls'),']'])
disp(['V_tib2_plus=[ ',num2str(V_tib2_pls'),']']);
                
Force=ForceTorqueCalculator_5DoF([],[],'done');

q1=SolC(:,1)';
q2=SolC(:,2)';
q3=SolC(:,3)';
q4=SolC(:,4)';
q5=SolC(:,5)';
dq1=SolC(:,6)';
dq2=SolC(:,7)';
dq3=SolC(:,8)';
dq4=SolC(:,9)';
dq5=SolC(:,10)';

p_tib2 = [ + L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5)
           - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + L_tib*cos(q2 + q4 + q5)];

[p_tib2_mnsX ] = FoottPositon(Q_mns,L_fem, L_tib)';
[p_tib2_plsX ] = FoottPositon(Q_plsX,L_fem, L_tib)';
       
       
v_hip =[dq1.*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) - dq2.*(L_fem*cos(q2 + q5) + L_tib*cos(q2 + q4 + q5)) + dq5.*(L_fem*cos(q1 + q5) - L_fem*cos(q2 + q5) + L_tib*cos(q1 + q3 + q5) - L_tib*cos(q2 + q4 + q5)) + L_tib*dq3.*cos(q1 + q3 + q5) - L_tib*dq4.*cos(q2 + q4 + q5)
        dq1.*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) - dq2.*(L_fem*sin(q2 + q5) + L_tib*sin(q2 + q4 + q5)) + dq5.*(L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5)) + L_tib*dq3.*sin(q1 + q3 + q5) - L_tib*dq4.*sin(q2 + q4 + q5)];
        
T_impact=[Tc(end)]; 

[KE_mns,PE_mns ] = KineticPotentialEnergy(Q_mns,DQ_minus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
[KE_pls,PE_pls ] = KineticPotentialEnergy(Q_plsX,DQ_plus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g );
disp(['K_mns= ',num2str(KE_mns'),',  K_pls= ',num2str(KE_pls')])
disp('---------------------')


% 
AnimBot3DOF([Tc ] ,[SolC(:,1:5)],T_impact,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,2)

figure(6)
    plot(Tc, Force(1:4,:)')
    legend('u_1','u_2','u_3','u_4')
    grid on
    xlabel('Time (s)')
    ylabel('Torque (N.m)')

    

figure(7)
    Thetas=c*[q1;q2;q3;q4;q5];
    Hds=zeros(4,1);
    DHds=zeros(4,1);
    Hdr=zeros(4,1);
    DHdr=zeros(4,1);
    Thetar=linspace(Theta_plsX,Theta_mnsX,length(Tc));
    for i=1:length(Thetar)
        [Hds(:,end+1),DHds(:,end+1)]=BezierFunction(Thetas(i),AlfaX,Theta_plsX,Theta_mnsX);
        [Hdr(:,end+1),DHdr(:,end+1)]=BezierFunction(Thetar(i),Alfa,Theta_plsX,Theta_mnsX);
    end
    Hds(:,1)=[];
    DHds(:,1)=[];
    Hdr(:,1)=[];
    DHdr(:,1)=[];
    subplot(3,2,1)
    plot(Tc,q1)
    hold all
    plot(Tc,Hds(1,:),'-.')
    plot(Tc,Hdr(1,:),'-.r')
    legend1=legend('sim','ref','des');
    set(legend1,'Orientation','horizontal');
    ylabel('q_1')
    hold off
    subplot(3,2,2)
    plot(Tc,q2)
    hold all
    plot(Tc,Hds(2,:),'-.')
    plot(Tc,Hdr(2,:),'-.r')
    legend1=legend('sim','ref','des');
    set(legend1,'Orientation','horizontal');
    hold off
    ylabel('q_2')
    subplot(3,2,3)
    plot(Tc,q3)
    hold all
    plot(Tc,Hds(3,:),'-.')
    plot(Tc,Hdr(3,:),'-.r')
    legend1=legend('sim','ref','des');
    set(legend1,'Orientation','horizontal');
    hold off
    ylabel('q_3')
    subplot(3,2,4)
    plot(Tc,q4)
    hold all
    plot(Tc,Hds(4,:),'-.')
    plot(Tc,Hdr(4,:),'-.r')
    legend1=legend('sim','ref','des');
    set(legend1,'Orientation','horizontal');
    ylabel('q_4')
    hold off
    subplot(3,2,5)
    plot(Tc,q5)
    ylabel('q_5')
    subplot(3,2,6)
    plot(Tc,Thetas)
    ylabel('\theta')
    

figure(8)
    plot(Tc,Force(5,:))
    hold all
    plot(Tc,Force(6,:))
    grid on
    legend('F_t','F_n')
    hold off
    xlabel('Time (s)')
    ylabel('GRF (N)')

    
figure(9)
    subplot(2,1,1)
    plot(Tc,p_tib2(1,:))
    hold all
    plot(Tc,p_tib2(2,:))
    grid on
    legend('p^x_t_2','p^y_t_2')
    hold off
    xlabel('Time (s)')
    ylabel('p_t_i_b_2 (m)')
    subplot(2,1,2)
    plot(Tc,v_hip(1,:))
    hold all
    plot(Tc,v_hip(2,:))
    grid on
    legend('v^x_h','v^y_h')
    hold off
    xlabel('Time (s)')
    ylabel('v_h_i_p (m/s)')
    
%% Multi Ode

% Kp=4000;
% Kd=130;


Q_pls=Q_plus;
DQ_pls=DQ_plus;
Alpha=Alfa;
Theta_mns= Theta_minus;
Theta_pls= Theta_plus;

T_impact=[]; 
Qq=[];
DQq=[];
Time=[0];
Force=[];
p_tib2 =[];
v_hip =[];
    
for step=1:2*10
    step
    
    
    Options = odeset('RelTol',1e-2,'AbsTol',1e-2,'maxstep',1e-2,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alpha,Theta_pls,Theta_mns,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));

    [Tc,SolC] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alpha,Theta_pls,Theta_mns,Kp,Kd),...
                        [0,2],[Q_pls; DQ_pls],Options);

    Q_mns =SolC(end,1:5 )';
    DQ_mns=SolC(end,6:10)';

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

    Force=[Force, ForceTorqueCalculator_5DoF([],[],'done')];

    q1=SolC(:,1)';
    q2=SolC(:,2)';
    q3=SolC(:,3)';
    q4=SolC(:,4)';
    q5=SolC(:,5)';
    dq1=SolC(:,6)';
    dq2=SolC(:,7)';
    dq3=SolC(:,8)';
    dq4=SolC(:,9)';
    dq5=SolC(:,10)';

    p_tib2 =[p_tib2,  [ + L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5),
                      - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + L_tib*cos(q2 + q4 + q5)]];


    v_hip =[v_hip , [dq1.*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) + dq5.*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) + L_tib*dq3.*cos(q1 + q3 + q5),
                     dq1.*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) + dq5.*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) + L_tib*dq3.*sin(q1 + q3 + q5)]];
        

    Qq =[Qq ;SolC(:,1:5)];
    DQq=[DQq ;SolC(:,6:10)];
    Time=[Time; Time(end)+Tc];
    T_impact=[T_impact Time(end)];
    
    
end
Time(1)=[];

AnimBot3DOF( Time, Qq,T_impact,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso,1)

%%

figure(3)
    plot(Time,Force(1:4,:)','linewidth',2)
    legend('u_1','u_2','u_3','u_4')
    grid on
    xlabel('Time (s)')
    ylabel('Torque (N.m)')
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','k','linestyle','-.','linewidth',1)   
    end
    hold off
    
    

figure(4)
    subplot(3,2,1)
    plot(Time,Qq(:,1))
    ylabel('q_1')
    grid on
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
    end
    hold off
    
    subplot(3,2,2)
    plot(Time,Qq(:,2))
    ylabel('q_2')
    grid on
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
    end
    hold off
    
    subplot(3,2,3)
    plot(Time,Qq(:,3))
    ylabel('q_3')
    grid on
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
    end
    hold off
    
    subplot(3,2,4)
    plot(Time,Qq(:,4))
    ylabel('q_4')
    grid on
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
    end
    hold off
    
    subplot(3,2,5)
    plot(Time,Qq(:,5))
    ylabel('q_5')
    grid on
    hold on
    h=get(gca,'YLIM')';
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.')
        drawnow
    end
    drawnow
    hold off
    
    subplot(3,2,6)
    plot(0,0)
    

figure(5)
    plot(Time,Force(5,:))
    hold all
    plot(Time,Force(6,:))
    grid on
    legend('F_t','F_n')
    hold off
    xlabel('Time (s)')
    ylabel('GRF (N)')
        h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.','linewidth',2)   
    end
    hold off


figure(6)
    subplot(2,1,1)
    plot(Time,p_tib2(1,:))
    hold all
    plot(Time,p_tib2(2,:))
    grid on
    legend('p^x_t_2','p^y_t_2')
    hold off
    xlabel('Time (s)')
    ylabel('p_t_i_b_2 (m)')
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.','linewidth',2)   
    end
    hold off
    
    subplot(2,1,2)
    plot(Time,v_hip(1,:))
    hold all
    plot(Time,v_hip(2,:))
    grid on
    legend('v^x_h','v^y_h')
    hold off
    xlabel('Time (s)')
    ylabel('v_h_i_p (m/s)')
    h=get(gca,'YLIM')';
    hold on
    for i=1:length(T_impact)
        a=line([T_impact(i) T_impact(i)],h);
        set(a,'color','r','linestyle','-.','linewidth',2)   
    end
    hold off
    
figure(7)
    subplot(3,2,1)
    plot(Qq(:,1),DQq(:,1))
    grid on
    hold on
    for i=1:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,1),DQq(index,1),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('q_1')
    ylabel('dq_1')
    
    subplot(3,2,2)
    plot(Qq(:,2),DQq(:,2))
    grid on
    hold on
    for i=1:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,2),DQq(index,2),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('q_2')
    ylabel('dq_2')
    
    
    subplot(3,2,3)
    plot(Qq(:,3),DQq(:,3))
    grid on
    hold on
    for i=1:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,3),DQq(index,3),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('q_3')
    ylabel('dq_3')
    
    subplot(3,2,4)
    plot(Qq(:,4),DQq(:,4))
    grid on
    hold on
    for i=1:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,4),DQq(index,4),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('q_4')
    ylabel('dq_4')
    
    subplot(3,2,5)
    plot(Qq(:,5),DQq(:,5))
    grid on
    hold on
    for i=1:length(T_impact)
        index=find(Time==T_impact(i),1);
        plot(Qq(index,5),DQq(index,5),'color','r','linestyle','none','marker','*')
    end
    hold off
    xlabel('q_5')
    ylabel('dq_5')

figure(8)
    
    indx1=1;
    for i=1:length(T_impact)/2
        indx2=find(Time==T_impact(2*(i)-1),1);
        indx3=find(Time==T_impact(2*(i)-0),1);
    
        subplot(3,2,1)
        plot(Qq(indx1:indx2,1),Force(1,indx1:indx2))
        hold on
        plot(Qq(indx2+1:indx3,2),Force(2,indx2+1:indx3),'r')
        grid on
%         hold off
        xlabel('q_1 without relabling (rad)')
        ylabel('u_1 without relabling (N.m)')

        subplot(3,2,2)
        plot(Qq(indx1:indx2  ,2),Force(2,indx1:indx2))
        hold on
        plot(Qq(indx2+1:indx3,1),Force(1,indx2+1:indx3),'r')
        grid on
%         hold off
        xlabel('q_2 without relabling (rad)')
        ylabel('u_2 without relabling (N.m)')


        subplot(3,2,3)
        plot(Qq(indx1:indx2  ,3),Force(3,indx1:indx2))
        hold on
        plot(Qq(indx2+1:indx3,4),Force(4,indx2+1:indx3),'r')
        grid on
%         hold off
        xlabel('q_3 without relabling (rad)')
        ylabel('u_3 without relabling (N.m)')

        subplot(3,2,4)
        plot(Qq(indx1:indx2  ,4),Force(4,indx1:indx2))
        hold on
        plot(Qq(indx2+1:indx3,3),Force(3,indx2+1:indx3),'r')
        grid on
%         hold off
        xlabel('q_4 without relabling (rad)')
        ylabel('u_4 without relabling (N.m)')

        subplot(3,2,5)
        plot(Qq(indx1:indx2  ,5),Force(5,indx1:indx2))
        hold on
        plot(Qq(indx2+1:indx3,5),Force(5,indx2+1:indx3),'r')
        grid on
%         hold off
        xlabel('q_5 without relabling (rad)')
        ylabel('u_5 without relabling (N.m)')

       indx1=indx2+1;
    end
    subplot(3,2,1)
    hold off
    subplot(3,2,2)
    hold off
    subplot(3,2,3)
    hold off
    subplot(3,2,4)
    hold off
    subplot(3,2,5)
    hold off
    subplot(3,2,6)
    hold off

figure(9)    
    i=length(T_impact)/2;
    indx1=find(Time==T_impact(2*(i)-2),1)+1;
    indx2=find(Time==T_impact(2*(i)-1),1);
    indx3=find(Time==T_impact(2*(i)-0),1);
    
    subplot(3,2,1)
    plot(Qq(indx1,1),Force(1,indx1),'marker','*')
    hold on
    plot(Qq(indx1:indx2,1),Force(1,indx1:indx2))
    plot(Qq(indx2+1,2),Force(2,indx2+1),'r','marker','*')
    plot(Qq(indx2+1:indx3,2),Force(2,indx2+1:indx3),'r')
    grid on
    hold off
    xlabel('q_1 without relabling (rad)')
    ylabel('u_1 without relabling (N.m)')

    subplot(3,2,2)
    plot(Qq(indx1  ,2),Force(2,indx1),'marker','*')
    hold on
    plot(Qq(indx1:indx2  ,2),Force(2,indx1:indx2))
    plot(Qq(indx2+1,1),Force(1,indx2+1),'r','marker','*')
    plot(Qq(indx2+1:indx3,1),Force(1,indx2+1:indx3),'r')
    grid on
    hold off
    xlabel('q_2 without relabling (rad)')
    ylabel('u_2 without relabling (N.m)')


    subplot(3,2,3)
    plot(Qq(indx1 ,3),Force(3,indx1),'marker','*')
    hold on
    plot(Qq(indx1:indx2  ,3),Force(3,indx1:indx2))
    plot(Qq(indx2+1,4),Force(4,indx2+1),'r','marker','*')
    plot(Qq(indx2+1:indx3,4),Force(4,indx2+1:indx3),'r')
    grid on
    hold off
    xlabel('q_3 without relabling (rad)')
    ylabel('u_3 without relabling (N.m)')

    subplot(3,2,4)
    plot(Qq(indx1  ,4),Force(4,indx1),'marker','*')
    hold on
    plot(Qq(indx1:indx2  ,4),Force(4,indx1:indx2))
    plot(Qq(indx2+1,3),Force(3,indx2+1),'r','marker','*')
    plot(Qq(indx2+1:indx3,3),Force(3,indx2+1:indx3),'r')
    grid on
    hold off
    xlabel('q_4 without relabling (rad)')
    ylabel('u_4 without relabling (N.m)')

    subplot(3,2,5)
    plot(Qq(indx1 ,5),Force(5,indx1),'marker','*')
    hold on
    plot(Qq(indx1:indx2  ,5),Force(5,indx1:indx2))
    plot(Qq(indx2+1,5),Force(5,indx2+1),'r','marker','*')
    plot(Qq(indx2+1:indx3,5),Force(5,indx2+1:indx3),'r')
    grid on
    hold off
    xlabel('q_5 without relabling (rad)')
    ylabel('u_5 without relabling (N.m)')

   