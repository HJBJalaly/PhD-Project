function Cost=CostFunctionSimGA(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Kp,Kd,rM,rB,rD,Landa,Gamma,Weight,SampleRate)

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

%%
Options = odeset('RelTol',5e-3,'AbsTol',5e-3,'maxstep',5e-3,...
                 'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Ma,Kp,Kd),...
                 'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));
            
[Tc,SolC] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Ma,Kp,Kd),...
                    [0,2],[Q_plus; DQ_plus],Options);


p_tib2_minus = FoottPositon(SolC(end,1:5 )',L_fem, L_tib)';

if(p_tib2_minus(1)<0)
    Cost=2e50;
    return
end

% else:

%%
Force=ForceTorqueCalculator_5DoF([],[],'done');

Ur1=Force(1,:)';
Ur1=[Ur1;Force(2,:)'];
Ur3=Force(3,:)';
Ur3=[Ur3;Force(4,:)'];

Q1=SolC(:,1);
Q1=[Q1;SolC(:,2)];
Q3=SolC(:,3);
Q3=[Q3;SolC(:,4)];

Qhat13=Q1+Q3;

DtQ1=SolC(:,6);
DtQ1=[DtQ1;SolC(:,7)];
DtQ3=SolC(:,8);
DtQ3=[DtQ3;SolC(:,9)];

Cost=OptimalParam([Tc ;Tc+Tc(end)],Q1,Q3,Qhat13,DtQ1,DtQ3,Ur1,Ur3,rM,rB,rD,Landa,Gamma,Weight,SampleRate,0,0);
                
% Ua1=(Ur1 -TorqueMonoSpring1 -TorqueDamper1 -TorqueBiSpring13 );
% Ua3=(Ur3 -TorqueMonoSpring3 -TorqueDamper3 -TorqueBiSpring13 );

% 
% Work_r=sum(abs(Ur1.*DtQ1).*[diff(TimeStride) ;0])+sum(abs(Ur3.*DtQ3).*[diff(TimeStride) ;0]);
% Work_a=sum(abs(Ua1.*DtQ1).*[diff(TimeStride) ;0])+sum(abs(Ua3.*DtQ3).*[diff(TimeStride) ;0]);

% CostReq=sum((Ur1.*Ur1+Ur3.*Ur3).*[diff(TimeStride);0])/2;
% 
% Cost
%  Cost=sum((Force(1,:).*Force(1,:)+Force(3,:).*Force(3,:)+...
%                Force(2,:).*Force(2,:)+Force(4,:).*Force(4,:))'.*[diff(Tc) ;0]/2)


end