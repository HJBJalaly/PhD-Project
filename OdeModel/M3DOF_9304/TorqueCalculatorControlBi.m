function status=TorqueCalculatorControlBi(t,Y,flag,g,L,m,Qref,DQref,Time,Kp,Kd,PassiveParam,rU,rB)

persistent Torque;

if(strcmp(flag,'init'))
    Torque=[];
    t=t(1);
elseif(strcmp(flag,'done'))
    status=Torque;
    return;
else
    status=0;
end

% t
% q1=min(max(Y(1),-pi),pi);
% q2=min(max(Y(2),-pi),pi);
% q3=min(max(Y(3),-pi),pi);

q1=Y(1);
q2=Y(2);
q3=Y(3);

D1q1=Y(4);
D1q2=Y(5);
D1q3=Y(6);
D1Q=[D1q1 D1q2 D1q3]';

mL1=m;
mL2=m;
mL3=m;

LL1=L;
LL2=L;
LL3=L;


%% Dynamic
MM=[mL3*LL1*LL3*cos(q2+q3)+LL1*LL2*cos(q2)*mL2+2*LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL3*LL1^2+LL1^2*mL2+mL1*LL1^2/3+mL2*LL2^2/3,mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6;
    mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6;
    mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6,mL3*LL3^2/3];


CC= [-mL3*LL1*LL3*sin(q2+q3)*(D1q2+D1q3)/2-LL2*(LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL1*LL3*(D1q1+D1q2+D1q3)*sin(q2+q3)/2-LL2*(sin(q2)*D1q1*LL1*(mL2+2*mL3)+LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL3*(LL1*sin(q2+q3)+LL2*sin(q3))*(D1q1+D1q2+D1q3)/2;
      mL3*LL1*D1q1*LL3*sin(q2+q3)/2+(sin(q2)*D1q1*LL1*(mL2+2*mL3)-mL3*LL3*sin(q3)*D1q3)*LL2/2,-mL3*LL2*LL3*sin(q3)*D1q3/2,-mL3*LL3*sin(q3)*(D1q1+D1q2+D1q3)*LL2/2;
      mL3*(LL1*D1q1*sin(q2+q3)+LL2*sin(q3)*(D1q1+D1q2))*LL3/2,LL3*LL2*sin(q3)*(D1q1+D1q2)*mL3/2,0];


GG= [g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2+2*cos(q1)*LL1*mL3+2*cos(q1)*LL1*mL2+cos(q1)*LL1*mL1)/2;
     g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2)/2;
     g*mL3*LL3*cos(q1+q2+q3)/2];


%%  FeedForward Law Control
TorquePassiveQ1Uni=polyval(PassiveParam(0*(rU+1)+1:1*(rU+1)),q1);
TorquePassiveQ2Uni=polyval(PassiveParam(1*(rU+1)+1:2*(rU+1)),q2);
TorquePassiveQ3Uni=polyval(PassiveParam(2*(rU+1)+1:3*(rU+1)),q3);
TorquePassiveUni=[TorquePassiveQ1Uni; TorquePassiveQ2Uni; TorquePassiveQ3Uni];

TorquePassiveQ12Bi=polyval(PassiveParam(3*(rU+1)+0*(rB+1)+1:3*(rU+1)+1*(rB+1)),q1+q2);
TorquePassiveQ23Bi=polyval(PassiveParam(3*(rU+1)+1*(rB+1)+1:3*(rU+1)+2*(rB+1)),q2+q3);
TorquePassiveBi=[ TorquePassiveQ12Bi; TorquePassiveQ23Bi];

% TorqueActiveQ1=interp1(Time,FF_Torque(1,:),t,'spline');
% TorqueActiveQ2=interp1(Time,FF_Torque(2,:),t,'spline');
% TorqueActiveQ3=interp1(Time,FF_Torque(3,:),t,'spline');
% TorqueActive=[TorqueActiveQ1; TorqueActiveQ2; TorqueActiveQ3];

Q1ref=interp1(Time,Qref(1,:),t,'spline');
Q2ref=interp1(Time,Qref(2,:),t,'spline');
Q3ref=interp1(Time,Qref(3,:),t,'spline');
DQ1ref=interp1(Time,DQref(1,:),t,'spline');
DQ2ref=interp1(Time,DQref(2,:),t,'spline');
DQ3ref=interp1(Time,DQref(3,:),t,'spline');

Torque= [Torque,[Kp*([Q1ref Q2ref Q3ref]'-[q1 q2 q3]')+Kd*([DQ1ref DQ2ref DQ3ref]'-D1Q); TorquePassiveUni;[TorquePassiveBi;0]+[0;TorquePassiveBi]]];
    
 

end