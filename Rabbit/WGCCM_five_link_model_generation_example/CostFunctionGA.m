function Cost=CostFunctionGA(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Kp,Kd,mu)

%Param=[ vector(Alpha(:,3:M-1)) ; Q_mius ; DQ_mius]

angl1=Param(end-7);
angl2=Param(end-6);
qTorso=Param(end-5);

TT = [1   0   0   0  -1;
      0   1   0   0  -1;
     -1   0   1   0   0;
	  0  -1   0   1   0;
	  0   0   0   0   1];

Q_minus=(TT*deg2rad([180-angl1+angl2 , 180+angl1+angl2 , 180-angl1+angl2-2*angl2 , 180+angl1+angl2-2*angl2 , qTorso])'); % This Q_minus;
DQ_minus=Param(end-4:end)';

c=[-1 0 -1/2 0 -1];

[Q_plus,DQ_plus]=ImpactModel(Q_minus,DQ_minus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
Theta_plus=c*Q_plus;
Theta_minus=c*Q_minus;
DTheta_minus=c*DQ_minus;
DTheta_plus=c*DQ_plus;

Alfa=zeros(4,Ma+1);
Alfa(:,1)=Q_plus(1:4);
Alfa(:,2)=DQ_plus(1:4)/(Ma*DTheta_plus)*(Theta_minus-Theta_plus)+Alfa(:,1);
Alfa(:,Ma+1)=Q_minus(1:4);
Alfa(:,Ma)=-DQ_minus(1:4)/(Ma*DTheta_minus)*(Theta_minus-Theta_plus)+Alfa(:,Ma+1);
Alfa(:,2+1:Ma-2+1)=reshape(Param(1:(Ma-3)*4),4,Ma-3);

for i=1:4

    Options = odeset('RelTol',1e-3,'AbsTol',1e-3,'maxstep',1e-3,...
                     'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Ma,Kp,Kd),...
                     'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));

    [Tc,SolC] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Ma,Kp,Kd),...
                        [0,2],[Q_plus; DQ_plus],Options);

    p_tib2_minus = FoottPositon(SolC(end,1:5 )',L_fem, L_tib)';
    p_tib2_plus = FoottPositon(Q_plus,L_fem, L_tib)';

    if(Tc(end)>1.99 || any(abs(SolC(end,1))>2*pi))
        Cost=2e65;
        return
    end

    if(p_tib2_minus(1)<(-.75*p_tib2_plus(1)))
        Cost=2e50;
        return
    end

    if(any(SolC(:,5)>0))
        Cost=1e50;
        return
    end
    
    Q_minus =SolC(end,1:5 )';
    DQ_minus=SolC(end,6:10)';
    [Q_plus,DQ_plus]=ImpactModel(Q_minus,DQ_minus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
    Theta_plus=c*Q_plus;
    Theta_minus=c*Q_minus;
end

Force=ForceTorqueCalculator_5DoF([],[],'done');

Cost=sum((Force(1,:).*Force(1,:)+Force(3,:).*Force(3,:)+...
              Force(2,:).*Force(2,:)+Force(4,:).*Force(4,:))'.*[diff(Tc) ;0]/2);

if( any( abs(Force(5,:))>(abs(Force(6,:))*mu))) % Ft> Fn*mu
    Cost=Cost+5e6;
    return
end

q1=SolC(:,1)';q2=SolC(:,2)';q3=SolC(:,3)';q4=SolC(:,4)';q5=SolC(:,5)';
dq1=SolC(:,6)';dq2=SolC(:,7)';dq3=SolC(:,8)';dq4=SolC(:,9)';dq5=SolC(:,10)';

v_hip =[dq1.*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) - dq2.*(L_fem*cos(q2 + q5) + L_tib*cos(q2 + q4 + q5)) + dq5.*(L_fem*cos(q1 + q5) - L_fem*cos(q2 + q5) + L_tib*cos(q1 + q3 + q5) - L_tib*cos(q2 + q4 + q5)) + L_tib*dq3.*cos(q1 + q3 + q5) - L_tib*dq4.*cos(q2 + q4 + q5)
             dq1.*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) - dq2.*(L_fem*sin(q2 + q5) + L_tib*sin(q2 + q4 + q5)) + dq5.*(L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5)) + L_tib*dq3.*sin(q1 + q3 + q5) - L_tib*dq4.*sin(q2 + q4 + q5)];

Vel_hip_x_avg=sum(v_hip(1,1:end-1).*diff(Tc'))/Tc(end);
if( ~(Vel_hip_x_avg>1 && Vel_hip_x_avg <1.2))
    Cost=Cost+300000;
end
    
end