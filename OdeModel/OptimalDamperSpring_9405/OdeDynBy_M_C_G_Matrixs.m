

%% Ode in M, D and C Matrix for my robot for 3 DOF :
home 
clear
syms mL1 mL2 mL3  LL1 LL2 LL3 g real positive
syms q1 q2 q3 real
syms D1q1 D1q2 D1q3 real



MM=[mL3*LL1*LL3*cos(q2+q3)+LL1*LL2*cos(q2)*mL2+2*LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL3*LL1^2+LL1^2*mL2+mL1*LL1^2/3+mL2*LL2^2/3,mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6;
    mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6;
    mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6,mL3*LL3^2/3];


CC= [-mL3*LL1*LL3*sin(q2+q3)*(D1q2+D1q3)/2-LL2*(LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL1*LL3*(D1q1+D1q2+D1q3)*sin(q2+q3)/2-LL2*(sin(q2)*D1q1*LL1*(mL2+2*mL3)+LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL3*(LL1*sin(q2+q3)+LL2*sin(q3))*(D1q1+D1q2+D1q3)/2;
      mL3*LL1*D1q1*LL3*sin(q2+q3)/2+(sin(q2)*D1q1*LL1*(mL2+2*mL3)-mL3*LL3*sin(q3)*D1q3)*LL2/2,-mL3*LL2*LL3*sin(q3)*D1q3/2,-mL3*LL3*sin(q3)*(D1q1+D1q2+D1q3)*LL2/2;
      mL3*(LL1*D1q1*sin(q2+q3)+LL2*sin(q3)*(D1q1+D1q2))*LL3/2,LL3*LL2*sin(q3)*(D1q1+D1q2)*mL3/2,0];


GG= [g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2+2*cos(q1)*LL1*mL3+2*cos(q1)*LL1*mL2+cos(q1)*LL1*mL1)/2;
     g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2)/2;
     g*mL3*LL3*cos(q1+q2+q3)/2];


whos

%% Ode in M, D and C Matrix for my robot for 2 DOF :
home 
clear
syms mL1 mL2   LL1 LL2 g real positive
syms q1 q2 real
syms D1q1 D1q2 real



MM=[mL2*LL1*LL2*cos(q2)+mL2*LL1^2+mL2*LL2^2/3+mL1*LL1^2/3, mL2*LL2*(3*LL1*cos(q2)+2*LL2)/6;
    mL2*LL2*(3*LL1*cos(q2)+2*LL2)/6, mL2*LL2^2/3];


CC= [-mL2*LL2*sin(q2)*D1q2*LL1/2, -mL2*LL2*sin(q2)*(D1q1+D1q2)*LL1/2; 
      mL2*LL2*LL1*D1q1*sin(q2)/2, 0];


GG= [g*(mL2*LL2*cos(q1+q2)+LL1*cos(q1)*mL1+2*LL1*cos(q1)*mL2)/2;
     g*mL2*LL2*cos(q1+q2)/2];


whos


%% Lagrange for 3 DOF :
home 
clear
syms mL1 mL2 mL3  LL1 LL2 LL3 g real positive
syms q1 q2 q3 real
syms D1q1 D1q2 D1q3 real

LL=mL2*LL1*D1q1*LL2*D1q2*cos(q2)/2 + mL1*LL1^2*D1q1^2/6 + mL2*LL2^2*D1q1*D1q2/3 + ...
   mL2*LL1^2*D1q1^2/2 + mL2*LL2^2*D1q1^2/6 + mL2*LL2^2*D1q2^2/6 + mL2*LL1*D1q1^2*LL2*cos(q2)/2 + ...
   mL3*LL1*D1q1^2*LL2*cos(q2) + mL3*LL1*D1q1^2*LL3*cos(q2+q3)/2 + mL3*LL2*D1q1^2*LL3*cos(q3)/2 + ...
   mL3*LL2^2*D1q1*D1q2 + mL3*LL3^2*D1q2*D1q3/3 + mL3*LL3^2*D1q1*D1q3/3 + mL3*LL3^2*D1q1*D1q2/3 + ...
   mL3*LL1^2*D1q1^2/2 + mL3*LL2^2*D1q1^2/2 + mL3*LL2^2*D1q2^2/2 + mL3*LL3^2*D1q1^2/6 + ...
   mL3*LL3^2*D1q2^2/6 + mL3*LL3^2*D1q3^2/6 + mL3*LL1*D1q1*LL2*D1q2*cos(q2) + ...
   mL3*LL1*D1q1*LL3*D1q2*cos(q2+q3)/2 + mL3*LL1*D1q1*LL3*D1q3*cos(q2+q3)/2 + ...
   mL3*LL2*D1q1*LL3*D1q2*cos(q3) + mL3*LL2*D1q1*LL3*D1q3*cos(q3)/2 + mL3*LL2*D1q2*LL3*D1q3*cos(q3)/2 + ...
   mL3*LL2*D1q2^2*LL3*cos(q3)/2 - g*(mL3*LL3*sin(q1+q2+q3) + LL2*(2*mL3+mL2)*sin(q1+q2) + ...
   sin(q1)*LL1*(mL1 + 2*mL2 + 2*mL3))/2;

