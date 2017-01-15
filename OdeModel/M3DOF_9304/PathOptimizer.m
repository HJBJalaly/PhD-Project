function Cost=PathOptimizer(Alpha,Beta,Theta,Time,Degree,Tres,Weight,g,mL1,mL2,mL3,LL1,LL2,LL3)
   
% Torque=zeros(3,length(Time));
nn=Degree(1);
rQ=Degree(2);
rM=Degree(3);
rB=Degree(4);

Alpha_Q1=Alpha(1:(rQ+1));
Alpha_Q2=Alpha((rQ+1)+1:2*(rQ+1));
Alpha_Q3=Alpha(2*(rQ+1)+1:3*(rQ+1));

Alpha_D1Q1=Alpha_Q1(1:end-1).*(rQ:-1:1);
Alpha_D1Q2=Alpha_Q2(1:end-1).*(rQ:-1:1);
Alpha_D1Q3=Alpha_Q3(1:end-1).*(rQ:-1:1);

Alpha_D2Q1=Alpha_D1Q1(1:end-1).*(rQ-1:-1:1);
Alpha_D2Q2=Alpha_D1Q2(1:end-1).*(rQ-1:-1:1);
Alpha_D2Q3=Alpha_D1Q3(1:end-1).*(rQ-1:-1:1);

Q1=polyval(Alpha_Q1,Time);
Q2=polyval(Alpha_Q2,Time);
Q3=polyval(Alpha_Q3,Time);
QJ=[Q1;Q2;Q3];

Qhat1=Q1+Q2;
Qhat2=Q2+Q3;
QhatJ=[Qhat1;Qhat2];


D1Q1=polyval(Alpha_D1Q1,Time);
D1Q2=polyval(Alpha_D1Q2,Time);
D1Q3=polyval(Alpha_D1Q3,Time);

D2Q1=polyval(Alpha_D2Q1,Time);
D2Q2=polyval(Alpha_D2Q2,Time);
D2Q3=polyval(Alpha_D2Q3,Time);

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


     TorqueDesire(:,i) = MM*[D2q1;D2q2;D2q3] + CC*[D1q1;D1q2;D1q3] + GG;
    
end

rMp=rM+1;
rBp=rB+1;
for i=1:nn-1
    TorqueMonoOptimal(i,:)=polyval(Beta((i-1)*(rMp)+1:i*(rMp)),QJ(i,:));
    TorqueBicepsOptimal(i,:)=polyval(Theta((i-1)*(rBp)+1:(i)*(rBp)),QhatJ(i,:));
end
i=nn;
TorqueMonoOptimal(i,:)=polyval(Beta((i-1)*(rMp)+1:i*(rMp)),QJ(i,:));

TorqueActive=TorqueDesire-TorqueMonoOptimal-[TorqueBicepsOptimal;zeros(size(Q1))]-[zeros(size(Q1));TorqueBicepsOptimal];
    
Cost=0;
for i=1:nn
    Cost=Cost + ...
          Weight(i)*( 1/2*(TorqueActive(i,:)')'*TorqueActive(i,:)')*Tres;
end

end