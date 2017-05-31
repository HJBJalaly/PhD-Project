function Cost=CostFunctionGA(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Kp,Kd)

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

% for i=2+1:Ma-2+1
%     Alfa(:,i)=(0*Alfa(:,2)+50*Alfa(:,Ma-1+1))/50;
% end

Options = odeset('RelTol',5e-3,'AbsTol',5e-3,'maxstep',5e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Ma,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));
            
[Tc,SolC] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Ma,Kp,Kd),...
                    [0,2],[Q_plus; DQ_plus],Options);

Force=ForceTorqueCalculator_5DoF([],[],'done');

Cost=sum((Force(1,:).*Force(1,:)+Force(3,:).*Force(3,:)+...
              Force(2,:).*Force(2,:)+Force(4,:).*Force(4,:))'.*[diff(Tc) ;0]/2);


p_tib2_minus = FoottPositon(SolC(end,1:5 )',L_fem, L_tib)';

if(p_tib2_minus(1)<0)
    Cost=Cost+20000;
end
          
end