function [Cneq,Ceq]=CF1_NonLinearConstraintWithSplopeLimit( Coef,time,Tres,Degree,Weight,Xef,Yef,g,mL1,mL2,mL3,LL1,LL2,LL3)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Pos=[Xef;Yef];

CoefP_q1=Coef(1:(Degree+1));
CoefP_q2=Coef((Degree+1)+1:2*(Degree+1));
CoefP_q3=Coef(2*(Degree+1)+1:3*(Degree+1));

D1CoefP_q1=CoefP_q1(1:end-1).*(Degree:-1:1);
D1CoefP_q2=CoefP_q2(1:end-1).*(Degree:-1:1);
D1CoefP_q3=CoefP_q3(1:end-1).*(Degree:-1:1);

D2CoefP_q1=D1CoefP_q1(1:end-1).*(Degree-1:-1:1);
D2CoefP_q2=D1CoefP_q2(1:end-1).*(Degree-1:-1:1);
D2CoefP_q3=D1CoefP_q3(1:end-1).*(Degree-1:-1:1);

Q1=polyval(CoefP_q1,time);
Q2=polyval(CoefP_q2,time);
Q3=polyval(CoefP_q3,time);

D1Q1=polyval(D1CoefP_q1,time);
D1Q2=polyval(D1CoefP_q2,time);
D1Q3=polyval(D1CoefP_q3,time);

D2Q1=polyval(D2CoefP_q1,time);
D2Q2=polyval(D2CoefP_q2,time);
D2Q3=polyval(D2CoefP_q3,time);

for i=1:length(time)
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


CostSlope=0;
Qq=[Q1;Q2;Q3];
for Joint=1:3
    
    Till=floor(size(Qq,2)/1);
    
    ThetaShift=Qq(Joint, 1:Till)-min(Qq(Joint, 1:Till));
    ThetaShiftScale = ThetaShift* (deg2rad(270) /  max(ThetaShift));
    
    tau=Torque(Joint, 1:Till)-min(Torque(Joint, 1:Till));
    tauShiftScale= tau /max(tau);
    
    DTa=diff(tauShiftScale)./diff(ThetaShiftScale);
	CostSlope(Joint)=max((DTa.^2)*Weight(Joint));
    

end


RPos=[LL1*cos(Q1)+LL2*cos(Q1+Q2)+LL3*cos(Q1+Q2+Q3);
      LL1*sin(Q1)+LL2*sin(Q1+Q2)+LL3*sin(Q1+Q2+Q3)];

    
Cneq =  [sum(sum((RPos-Pos).^2))*Tres-.0001;
         max(CostSlope)                     ];
    
Ceq  =  [RPos(:,1)-Pos(:,1);
         RPos(:,end)-Pos(:,end);
         Q1(1)-Q1(end);
         Q2(1)-Q2(end);
         Q3(1)-Q3(end);
         D1Q1(1)-D1Q1(end);
         D1Q2(1)-D1Q2(end);
         D1Q3(1)-D1Q3(end);
         D2Q1(1)-D2Q1(end);
         D2Q2(1)-D2Q2(end);
         D2Q3(1)-D2Q3(end)];    
end

