function Cost=CF1_TorqueCostWithoutSlope(Coef,Time,Degree,Tres,Select,Weight,Landa,g,mL1,mL2,mL3,LL1,LL2,LL3)

% Torque=zeros(3,length(Time));

CoefP_q1=Coef(1:(Degree+1));
CoefP_q2=Coef((Degree+1)+1:2*(Degree+1));
CoefP_q3=Coef(2*(Degree+1)+1:3*(Degree+1));

D1CoefP_q1=CoefP_q1(1:end-1).*(Degree:-1:1);
D1CoefP_q2=CoefP_q2(1:end-1).*(Degree:-1:1);
D1CoefP_q3=CoefP_q3(1:end-1).*(Degree:-1:1);

D2CoefP_q1=D1CoefP_q1(1:end-1).*(Degree-1:-1:1);
D2CoefP_q2=D1CoefP_q2(1:end-1).*(Degree-1:-1:1);
D2CoefP_q3=D1CoefP_q3(1:end-1).*(Degree-1:-1:1);

Q1=polyval(CoefP_q1,Time);
Q2=polyval(CoefP_q2,Time);
Q3=polyval(CoefP_q3,Time);

D1Q1=polyval(D1CoefP_q1,Time);
D1Q2=polyval(D1CoefP_q2,Time);
D1Q3=polyval(D1CoefP_q3,Time);

D2Q1=polyval(D2CoefP_q1,Time);
D2Q2=polyval(D2CoefP_q2,Time);
D2Q3=polyval(D2CoefP_q3,Time);

for i=1:length(Time)
    q1=Q1(i);
    q2=Q2(i);
    q3=Q3(i);
    D1q1=D1Q1(i);
    D1q2=D1Q2(i);
    D1q3=D1Q3(i);
    D2q1=D2Q1(i);
    D2q2=D2Q2(i);
    D2q3=D2Q3(i);

    
    MM=[mL3*LL1*LL3*cos(q2+q3)+LL1*LL2*cos(q2)*mL2+2*LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL3*LL1^2+LL1^2*mL2+mL1*LL1^2/3+mL2*LL2^2/3,mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6;
        mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6;
        mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6,mL3*LL3^2/3];

    CC= [-mL3*LL1*LL3*sin(q2+q3)*(D1q2+D1q3)/2-LL2*(LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL1*LL3*(D1q1+D1q2+D1q3)*sin(q2+q3)/2-LL2*(sin(q2)*D1q1*LL1*(mL2+2*mL3)+LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL3*(LL1*sin(q2+q3)+LL2*sin(q3))*(D1q1+D1q2+D1q3)/2;
          mL3*LL1*D1q1*LL3*sin(q2+q3)/2+(sin(q2)*D1q1*LL1*(mL2+2*mL3)-mL3*LL3*sin(q3)*D1q3)*LL2/2,-mL3*LL2*LL3*sin(q3)*D1q3/2,-mL3*LL3*sin(q3)*(D1q1+D1q2+D1q3)*LL2/2;
          mL3*(LL1*D1q1*sin(q2+q3)+LL2*sin(q3)*(D1q1+D1q2))*LL3/2,LL3*LL2*sin(q3)*(D1q1+D1q2)*mL3/2,0];

    GG= [g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2+2*cos(q1)*LL1*mL3+2*cos(q1)*LL1*mL2+cos(q1)*LL1*mL1)/2;
         g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2)/2;
         g*mL3*LL3*cos(q1+q2+q3)/2];


     Torque(:,i) = MM*[D2q1;D2q2;D2q3] + CC*[D1q1;D1q2;D1q3] + GG;
    
end

IntU2=sum(sum(Torque.^2,2).*Weight);
IntAbsUdq=sum(sum(abs(Torque.*[D1Q1;D1Q2;D1Q3]),2).*Weight);
IntUdq=sum(((sum((Torque.*[D1Q1;D1Q2;D1Q3]),2)).^2).*Weight);

CostArea=sum(Select.*[ IntU2 IntAbsUdq IntUdq])*Tres;
CostSlope=CostArea;
Cost=((Landa)*CostArea+(1-Landa)*CostSlope)/sum(Weight);