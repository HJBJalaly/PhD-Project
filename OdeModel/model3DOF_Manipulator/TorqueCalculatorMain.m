function Torque=TorqueCalculatorMain(D2q,Dq,q,g,mL1,mL2,mL3,LL1,LL2,LL3,Time,Kp,Kd)

Torque=zeros(size(q));

for i=1:size(q,2)
    q1=q(1,i);
    q2=q(2,i);
    D1q1=Dq(1,i);
    D1q2=Dq(2,i);
    

    MM=[mL2*LL1*LL2*cos(q2)+mL2*LL1^2+mL2*LL2^2/3+mL1*LL1^2/3, mL2*LL2*(3*LL1*cos(q2)+2*LL2)/6;
        mL2*LL2*(3*LL1*cos(q2)+2*LL2)/6, mL2*LL2^2/3];


    CC= [-mL2*LL2*sin(q2)*D1q2*LL1/2, -mL2*LL2*sin(q2)*(D1q1+D1q2)*LL1/2; 
          mL2*LL2*LL1*D1q1*sin(q2)/2, 0];


    GG= [g*(mL2*LL2*cos(q1+q2)+LL1*cos(q1)*mL1+2*LL1*cos(q1)*mL2)/2;
         g*mL2*LL2*cos(q1+q2)/2];

    Torque(:,i) = MM*D2q(:,i) + CC*Dq(:,i) + GG;
    
end