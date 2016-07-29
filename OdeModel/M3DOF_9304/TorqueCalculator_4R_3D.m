function Torque=TorqueCalculator_4R_3D(D2q,Dq,q,g,mL1,mL2,mL3,mL4,LL1,LL2,LL3,LL4)

Torque=zeros(size(q));

for i=1:size(q,2)
    q1=q(1,i);
    q2=q(2,i);
    q3=q(3,i);
    q4=q(4,i);
    D1q1=Dq(1,i);
    D1q2=Dq(2,i);
    D1q3=Dq(3,i);
    D1q4=Dq(4,i);
    

    MM= [ (63*mL1)/8000 + (LL3^2*mL3*sin(q2 + q3)^2)/4 + LL3^2*mL4*sin(q2 + q3)^2 + (LL2^2*mL2*sin(q2)^2)/4 + LL2^2*mL3*sin(q2)^2 + LL2^2*mL4*sin(q2)^2 + (LL4^2*mL4*sin(q2 + q3 + q4)^2)/4 - LL2*LL3*mL3*sin(q3/2)^2 - 2*LL2*LL3*mL4*sin(q3/2)^2 - LL3*LL4*mL4*sin(q4/2)^2 - LL2*LL4*mL4*sin(q3/2 + q4/2)^2 + LL3*LL4*mL4*sin(q2 + q3 + q4/2)^2 + LL2*LL4*mL4*sin(q2 + q3/2 + q4/2)^2 + LL2*LL3*mL3*sin(q2 + q3/2)^2 + 2*LL2*LL3*mL4*sin(q2 + q3/2)^2,                                                                                                                                                                                                                 0,                                                                                                                                                                              0,                                                                                         0;
                                                                                                                                                                                                                                                                                                                                                                                                                                            0, (73*mL1)/20000 + (LL1^2*mL1)/4 + (LL2^2*mL2)/4 + LL2^2*mL3 + LL2^2*mL4 + (LL3^2*mL3)/4 + LL3^2*mL4 + (LL4^2*mL4)/4 + LL2*LL4*mL4*cos(q3 + q4) + LL2*LL3*mL3*cos(q3) + 2*LL2*LL3*mL4*cos(q3) + LL3*LL4*mL4*cos(q4), (13*mL1)/8000 + (LL1^2*mL1)/6 + (LL3^2*mL3)/4 + LL3^2*mL4 + (LL4^2*mL4)/4 + (LL2*LL4*mL4*cos(q3 + q4))/2 + (LL2*LL3*mL3*cos(q3))/2 + LL2*LL3*mL4*cos(q3) + LL3*LL4*mL4*cos(q4), mL1/2500 + (LL1^2*mL1)/12 + (LL4^2*mL4)/4 + (LL2*LL4*mL4*cos(q3 + q4))/2 + (LL3*LL4*mL4*cos(q4))/2;
                                                                                                                                                                                                                                                                                                                                                                                                                                            0,                                    (13*mL1)/8000 + (LL1^2*mL1)/6 + (LL3^2*mL3)/4 + LL3^2*mL4 + (LL4^2*mL4)/4 + (LL2*LL4*mL4*cos(q3 + q4))/2 + (LL2*LL3*mL3*cos(q3))/2 + LL2*LL3*mL4*cos(q3) + LL3*LL4*mL4*cos(q4),                                                                                (13*mL1)/8000 + (LL1^2*mL1)/6 + (LL3^2*mL3)/4 + LL3^2*mL4 + (LL4^2*mL4)/4 + LL3*LL4*mL4*cos(q4),                                (mL1*LL1^2)/12 + (mL4*LL4^2)/4 + (LL3*mL4*cos(q4)*LL4)/2 + mL1/2500;
                                                                                                                                                                                                                                                                                                                                                                                                                                            0,                                                                                                                mL1/2500 + (LL1^2*mL1)/12 + (LL4^2*mL4)/4 + (LL2*LL4*mL4*cos(q3 + q4))/2 + (LL3*LL4*mL4*cos(q4))/2,                                                                                                            (mL1*LL1^2)/12 + (mL4*LL4^2)/4 + (LL3*mL4*cos(q4)*LL4)/2 + mL1/2500,                                                          (mL1*LL1^2)/12 + (mL4*LL4^2)/4 + mL1/2500 ];
 

    CC= [ D1q4*((LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4))/8 - (LL2*LL4*mL4*sin(q3 + q4))/4 - (LL3*LL4*mL4*sin(q4))/4 + (LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4))/4 + (LL2*LL4*mL4*sin(2*q2 + q3 + q4))/4) + D1q3*((LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4))/8 + (LL3^2*mL3*sin(2*q2 + 2*q3))/8 + (LL3^2*mL4*sin(2*q2 + 2*q3))/2 - (LL2*LL4*mL4*sin(q3 + q4))/4 - (LL2*LL3*mL3*sin(q3))/4 - (LL2*LL3*mL4*sin(q3))/2 + (LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4))/2 + (LL2*LL3*mL3*sin(2*q2 + q3))/4 + (LL2*LL3*mL4*sin(2*q2 + q3))/2 + (LL2*LL4*mL4*sin(2*q2 + q3 + q4))/4) + D1q2*((LL2^2*mL2*sin(2*q2))/8 + (LL2^2*mL3*sin(2*q2))/2 + (LL2^2*mL4*sin(2*q2))/2 + (LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4))/8 + (LL3^2*mL3*sin(2*q2 + 2*q3))/8 + (LL3^2*mL4*sin(2*q2 + 2*q3))/2 + (LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4))/2 + (LL2*LL3*mL3*sin(2*q2 + q3))/2 + LL2*LL3*mL4*sin(2*q2 + q3) + (LL2*LL4*mL4*sin(2*q2 + q3 + q4))/2), (D1q1*(LL2^2*mL2*sin(2*q2) + 4*LL2^2*mL3*sin(2*q2) + 4*LL2^2*mL4*sin(2*q2) + LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4) + LL3^2*mL3*sin(2*q2 + 2*q3) + 4*LL3^2*mL4*sin(2*q2 + 2*q3) + 4*LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4) + 4*LL2*LL3*mL3*sin(2*q2 + q3) + 8*LL2*LL3*mL4*sin(2*q2 + q3) + 4*LL2*LL4*mL4*sin(2*q2 + q3 + q4)))/8, (D1q1*(LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4) + LL3^2*mL3*sin(2*q2 + 2*q3) + 4*LL3^2*mL4*sin(2*q2 + 2*q3) - 2*LL2*LL4*mL4*sin(q3 + q4) - 2*LL2*LL3*mL3*sin(q3) - 4*LL2*LL3*mL4*sin(q3) + 4*LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4) + 2*LL2*LL3*mL3*sin(2*q2 + q3) + 4*LL2*LL3*mL4*sin(2*q2 + q3) + 2*LL2*LL4*mL4*sin(2*q2 + q3 + q4)))/8, (D1q1*LL4*mL4*(2*LL2*sin(2*q2 + q3 + q4) - 2*LL2*sin(q3 + q4) - 2*LL3*sin(q4) + 2*LL3*sin(2*q2 + 2*q3 + q4) + LL4*sin(2*q2 + 2*q3 + 2*q4)))/8;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            -(D1q1*(LL2^2*mL2*sin(2*q2) + 4*LL2^2*mL3*sin(2*q2) + 4*LL2^2*mL4*sin(2*q2) + LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4) + LL3^2*mL3*sin(2*q2 + 2*q3) + 4*LL3^2*mL4*sin(2*q2 + 2*q3) + 4*LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4) + 4*LL2*LL3*mL3*sin(2*q2 + q3) + 8*LL2*LL3*mL4*sin(2*q2 + q3) + 4*LL2*LL4*mL4*sin(2*q2 + q3 + q4)))/8,                                                                                                                                                                  - D1q4*((LL2*LL4*mL4*sin(q3 + q4))/2 + (LL3*LL4*mL4*sin(q4))/2) - D1q3*((LL2*LL4*mL4*sin(q3 + q4))/2 + (LL2*LL3*mL3*sin(q3))/2 + LL2*LL3*mL4*sin(q3)),                                                                                                                   - (D1q2*LL2*(LL4*mL4*sin(q3 + q4) + LL3*mL3*sin(q3) + 2*LL3*mL4*sin(q3)))/2 - (D1q3*LL2*(LL4*mL4*sin(q3 + q4) + LL3*mL3*sin(q3) + 2*LL3*mL4*sin(q3)))/2 - (D1q4*LL4*mL4*(LL2*sin(q3 + q4) + LL3*sin(q4)))/2,                                                                            -(LL4*mL4*(LL2*sin(q3 + q4) + LL3*sin(q4))*(D1q2 + D1q3 + D1q4))/2;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     -(D1q1*(LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4) + LL3^2*mL3*sin(2*q2 + 2*q3) + 4*LL3^2*mL4*sin(2*q2 + 2*q3) - 2*LL2*LL4*mL4*sin(q3 + q4) - 2*LL2*LL3*mL3*sin(q3) - 4*LL2*LL3*mL4*sin(q3) + 4*LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4) + 2*LL2*LL3*mL3*sin(2*q2 + q3) + 4*LL2*LL3*mL4*sin(2*q2 + q3) + 2*LL2*LL4*mL4*sin(2*q2 + q3 + q4)))/8,                                                                                                                                                                                                     D1q2*((LL2*LL4*mL4*sin(q3 + q4))/2 + (LL2*LL3*mL3*sin(q3))/2 + LL2*LL3*mL4*sin(q3)) - (D1q4*LL3*LL4*mL4*sin(q4))/2,                                                                                                                                                                                                                                                                                                 -(D1q4*LL3*LL4*mL4*sin(q4))/2,                                                                                                 -(LL3*LL4*mL4*sin(q4)*(D1q2 + D1q3 + D1q4))/2;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     -(D1q1*LL4*mL4*(2*LL2*sin(2*q2 + q3 + q4) - 2*LL2*sin(q3 + q4) - 2*LL3*sin(q4) + 2*LL3*sin(2*q2 + 2*q3 + q4) + LL4*sin(2*q2 + 2*q3 + 2*q4)))/8,                                                                                                                                                                                                                           D1q2*((LL2*LL4*mL4*sin(q3 + q4))/2 + (LL3*LL4*mL4*sin(q4))/2) + (D1q3*LL3*LL4*mL4*sin(q4))/2,                                                                                                                                                                                                                                                                                         (LL3*LL4*mL4*sin(q4)*(D1q2 + D1q3))/2,                                                                                                                                             0         ];
 

    GG= [                                                                                                                                         0;
          - g*mL4*(LL3*sin(q2 + q3) + LL2*sin(q2) + (LL4*sin(q2 + q3 + q4))/2) - g*mL3*((LL3*sin(q2 + q3))/2 + LL2*sin(q2)) - (LL2*g*mL2*sin(q2))/2;
                                                                - g*mL4*(LL3*sin(q2 + q3) + (LL4*sin(q2 + q3 + q4))/2) - (LL3*g*mL3*sin(q2 + q3))/2;
                                                                                                                   -(LL4*g*mL4*sin(q2 + q3 + q4))/2 ];

    Torque(:,i) = MM*D2q(:,i) + CC*Dq(:,i) + GG;
    
end