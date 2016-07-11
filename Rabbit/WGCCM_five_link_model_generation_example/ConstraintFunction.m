function [Cneq,Ceq]=ConstraintFunction(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Kp,Kd)


% Param=[ vector(Alpha(:,3:M-1)) ; Q_mius ; DQ_mius]

% Q_minus=Param(end-9:end-5);
Q_minus=deg2rad([ 183.5  215.3  -20.1  -19.0   -9.5])';
% DQ_minus=Param(end-4:end);
DQ_minus=deg2rad([ 5.5   -1.9  -65.0   27.8  -33.5])';


[Q_plus,DQ_plus,V_tib2_plus,F2]=ImpactModel(Q_minus,DQ_minus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);


c=[-1 0 -1/2 0 -1];
Theta_plus=c*Q_plus;
Theta_minus=c*Q_minus;
DTheta_plus=c*DQ_plus;
DTheta_minus=c*DQ_minus;


Alfa=zeros(4,Ma+1);

Alfa(:,0+1)=Q_plus(1:4);
Alfa(:,Ma+1)=Q_minus(1:4);
Alfa(:,1+1)=DQ_plus(1:4)/(Ma*DTheta_plus)*(Theta_minus-Theta_plus)+Alfa(:,0+1);
Alfa(:,Ma-1+1)=-DQ_minus(1:4)/(Ma*DTheta_minus)*(Theta_minus-Theta_plus)+Alfa(:,Ma+1);

Alfa(:,2+1:Ma-2+1)=reshape(Param(1:(Ma-3)*4),4,Ma-3);
    

MotionData=ForceTorqueCalculator_5DoF([],[],'done');
[p_tib2_plus ]  = FoottPositon(Q_plus,L_fem, L_tib)';
[p_tib2_minus ] = FoottPositon(Q_minus,L_fem, L_tib)'
 
Q_mns =MotionData(7:11,end);

p_hip_y_mins=- L_fem*cos(Q_mns(1)+ Q_mns(5)) - L_tib*cos(Q_mns(1) + Q_mns(3) + Q_mns(5));
p_hip_y_plus=- L_fem*cos(Q_plus(1)+ Q_plus(5)) - L_tib*cos(Q_plus(1) + Q_plus(3) + Q_plus(5));

RR=[ 0 1 0 0 0;
     1 0 0 0 0;
     0 0 0 1 0;
     0 0 1 0 0;
     0 0 0 0 1];

Cneq=[ Q_minus(1:2)-3*pi/2;
      -Q_minus(1:2)+pi/2;
       Q_minus(3:4);
      -Q_minus(3:4)-pi/2;
       Q_minus(5)-pi/4;
      -Q_minus(5);
       Q_minus(1)-Q_minus(2);
       DQ_minus-5;
      -DQ_minus-5;
   10*(p_tib2_plus(1)-p_tib2_minus(1))
      -p_hip_y_mins;
      -p_hip_y_plus;
      -V_tib2_plus(2)];

Ceq=[ RR*Q_mns-Q_plus;
      p_tib2_minus(2)];
Cneq=[-p_hip_y_mins;
      -p_hip_y_plus;];
Ceq=[];
% 
% [Q_pls,DQ_pls,V_tib2_pls,F2]=ImpactModel(Q_mns,DQ_mns,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
                

end