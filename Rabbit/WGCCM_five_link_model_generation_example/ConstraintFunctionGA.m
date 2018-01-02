function [Cneq,Ceq]=ConstraintFunctionGA(Param,Ma,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,mu)


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

[Q_plus,DQ_plus,V_tib2_plus,F2,DeltaDqMinus]=ImpactModelDeltaDq(Q_minus,DQ_minus,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso);
% [p_tib2_minus ] = FoottPositon(Q_minus,L_fem, L_tib)';
% [p_tib2_plus ] = FoottPositon(Q_plus,L_fem, L_tib)';

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

% q3 and q4

Hdr=zeros(4,1);
p_tib2_plus=[];
for Theta=linspace(Theta_plus,Theta_minus,10)
    Hdr(:,end+1)=BezierFunction_hd_Fast(Ma,Theta,Alfa,Theta_plus,Theta_minus);
    p_tib2_plus(:,end+1) = FoottPositon([Hdr(:,end) ;deg2rad( qTorso)],L_fem, L_tib);
end
Hdr(:,1)=[];


%stabilty
H0=eye(4,5);
H=[H0;c];
Hi=H^(-1);
qz=@(xi)Hi*[BezierFunction_hd_Fast(Ma,xi,Alfa,Theta_plus,Theta_minus);xi];
gamma0=@(qz)...
        [ XX_fem + XX_tib + 2*L_fem^2*M_fem + L_fem^2*M_tib + 2*L_tib^2*M_fem + L_fem^2*M_torso + 2*L_tib^2*M_tib + L_tib^2*M_torso + Lc_fem^2*M_fem + Lc_tib^2*M_tib - 2*L_fem*Lc_fem*M_fem - 2*L_tib*Lc_tib*M_tib - L_fem^2*M_tib*cos(qz(1) - qz(2)) + 4*L_fem*L_tib*M_fem*cos(qz(3)) + 2*L_fem*L_tib*M_tib*cos(qz(3)) + 2*L_fem*L_tib*M_torso*cos(qz(3)) - L_fem*L_torso*M_torso*cos(qz(1)) - 2*L_tib*Lc_fem*M_fem*cos(qz(3)) + L_fem*Lc_torso*M_torso*cos(qz(1)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3)), XX_fem + XX_tib + L_fem^2*M_tib + Lc_fem^2*M_fem + Lc_tib^2*M_tib - L_fem^2*M_tib*cos(qz(1) - qz(2)) + 2*L_fem*Lc_tib*M_tib*cos(qz(4)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)), XX_tib + 2*L_tib^2*M_fem + 2*L_tib^2*M_tib + L_tib^2*M_torso + Lc_tib^2*M_tib - 2*L_tib*Lc_tib*M_tib + 2*L_fem*L_tib*M_fem*cos(qz(3)) + L_fem*L_tib*M_tib*cos(qz(3)) + L_fem*L_tib*M_torso*cos(qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(3)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3)), XX_tib + Lc_tib^2*M_tib + L_fem*Lc_tib*M_tib*cos(qz(4)) - L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)), 2*XX_fem + 2*XX_tib + XX_torso + 2*L_fem^2*M_fem + 2*L_fem^2*M_tib + 2*L_tib^2*M_fem + L_fem^2*M_torso + 2*L_tib^2*M_tib + L_tib^2*M_torso + L_torso^2*M_torso + 2*Lc_fem^2*M_fem + 2*Lc_tib^2*M_tib + Lc_torso^2*M_torso - 2*L_fem*Lc_fem*M_fem - 2*L_tib*Lc_tib*M_tib - 2*L_torso*Lc_torso*M_torso - 2*L_fem^2*M_tib*cos(qz(1) - qz(2)) + 4*L_fem*L_tib*M_fem*cos(qz(3)) + 2*L_fem*L_tib*M_tib*cos(qz(3)) + 2*L_fem*L_tib*M_torso*cos(qz(3)) - 2*L_fem*L_torso*M_torso*cos(qz(1)) - 2*L_tib*Lc_fem*M_fem*cos(qz(3)) + 2*L_fem*Lc_tib*M_tib*cos(qz(4)) + 2*L_fem*Lc_torso*M_torso*cos(qz(1)) - 2*L_fem*Lc_tib*M_tib*cos(qz(1) - qz(2) - qz(4)) - 2*L_fem*Lc_fem*M_fem*cos(qz(1) - qz(2)) - 2*L_tib*Lc_tib*M_tib*cos(qz(1) - qz(2) + qz(3) - qz(4)) - 2*L_fem*L_tib*M_tib*cos(qz(1) - qz(2) + qz(3)) - 2*L_tib*Lc_fem*M_fem*cos(qz(1) - qz(2) + qz(3)) - 2*L_tib*L_torso*M_torso*cos(qz(1) + qz(3)) + 2*L_tib*Lc_torso*M_torso*cos(qz(1) + qz(3))];
kappa1=@(xi)c*inv([H0-BezierFunction_Dhd_Fast(Ma,xi,Alfa,Theta_plus,Theta_minus)*c/(Theta_minus-Theta_plus);gamma0(qz(xi))])*[zeros(4,1);1];

% PE=g*(M_tib*(2*y1 - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + Lc_tib*cos(q2 + q4 + q5)) - M_fem*(cos(q1 + q5)*(L_fem - Lc_fem) - 2*y1 + L_tib*cos(q1 + q3 + q5)) + M_fem*(2*y1 - L_fem*cos(q1 + q5) + Lc_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5)) + M_torso*(2*y1 - L_fem*cos(q1 + q5) + cos(q5)*(L_torso - Lc_torso) - L_tib*cos(q1 + q3 + q5)) + M_tib*(2*y1 - cos(q1 + q3 + q5)*(L_tib - Lc_tib)))
% kappa2=-diff(PE,q5) ;
kappa2=@(qz)g*(M_fem*(sin(qz(4) + qz(5))*(L_fem - Lc_fem) + L_tib*sin(qz(4) + qz(3) + qz(5))) + M_tib*(L_fem*sin(qz(4) + qz(5)) - L_fem*sin(qz(2) + qz(5)) + L_tib*sin(qz(4) + qz(3) + qz(5)) - Lc_tib*sin(qz(2) + qz(4) + qz(5))) + M_fem*(L_fem*sin(qz(4) + qz(5)) - Lc_fem*sin(qz(2) + qz(5)) + L_tib*sin(qz(4) + qz(3) + qz(5))) + M_torso*(L_fem*sin(qz(4) + qz(5)) - sin(qz(5))*(L_torso - Lc_torso) + L_tib*sin(qz(4) + qz(3) + qz(5))) + M_tib*sin(qz(4) + qz(3) + qz(5))*(L_tib - Lc_tib));
Vzero=0;
VzeroMax=-inf;
Dxi=-(Theta_plus-Theta_minus)/25;
for Xi=linspace(Theta_plus,Theta_minus,25)
    Vzero=Vzero-kappa2(qz(Xi))/kappa1(Xi);
    if(Vzero>VzeroMax)
        VzeroMax=Vzero;
    end
end
Vzero=Vzero*Dxi;
VzeroMax=VzeroMax*Dxi;
    
LandaHatDq=@(qm)inv([H0-Ma*(Alfa(:,Ma+1)-Alfa(:,Ma))*c/(Theta_minus-Theta_plus);gamma0(qm)])*[zeros(4,1);1];
deltaZero=gamma0(Q_plus)*DeltaDqMinus*LandaHatDq(Q_minus);

if(deltaZero>1)
    disp('2')
end
if( ((deltaZero^2)/(1-deltaZero^2)*Vzero+VzeroMax )>0)
    disp('3')
end

% for i=2+1:Ma-2+1
%     Alfa(:,i)=(0*Alfa(:,2)+50*Alfa(:,Ma-1+1))/50;
% end

% Options = odeset('RelTol',5e-3,'AbsTol',5e-3,'maxstep',5e-3,...
%                  'OutputFcn',@(t,x,flag)ForceTorqueCalculator_5DoF(t,x,flag,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Kp,Kd),...
%                  'Events',@(t,x)EventTouchDown(t,x,L_fem, L_tib));
%             
% [Tc,SolC] = ode15s(@(t,x)RabitDynamic(t,x,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g,Alfa,Theta_plus,Theta_minus,Kp,Kd),...
%                     [0,2],[Q_plus; DQ_plus],Options);
% 
% Force=ForceTorqueCalculator_5DoF([],[],'done');
% 
% Cost=sum((Force(1,:).*Force(1,:)+Force(3,:).*Force(3,:)+...
%               Force(2,:).*Force(2,:)+Force(4,:).*Force(4,:))'.*[diff(Tc) ;0]/2);
% 

Cneq=[ angl1-45;    % in deg
      -angl1;       % in deg
       angl2-45;    % in deg
      -angl2;       % in deg
       qTorso;      % in deg
       -qTorso-45;  % in deg
%        Q_minus(1)-Q_minus(2);
%        DQ_minus-5;
%       -DQ_minus-5;
%       p_tib2_plus(1);
%      -p_tib2_minus(1);
%       -p_hip_y_mins;
%       -p_hip_y_plus;
      -V_tib2_plus(2); % vertical velocity of foot after the impact moment
      -F2(2);
      -deltaZero^2; %stability condition
       deltaZero^2-1'; %stability condition
       (deltaZero^2)/(1-deltaZero^2)*Vzero+VzeroMax; %stability condition
       Hdr(3,:)'+0.035;%knee joint
       Hdr(4,:)'+0.035; % knee joint
      -p_tib2_plus(2,:)'-0.002; % position of swing foot over one step
       DQ_minus-10; % velocity
      -DQ_minus-10;
       abs(F2(1))-abs(F2(2))*mu;]; % non-slippinf condtion at impact

Ceq=[ ];

end