function Cost=CF3b_TorqueCost_4R_3D(Alpha,Alpha_Q1,Time,Degree,Tres,Weight,Landa,SampleRate,g,mL1,mL2,mL3,mL4,LL1,LL2,LL3,LL4)
   
% Torque=zeros(3,length(Time));
nn=Degree(1);
rQ=Degree(2);
rU=Degree(3);
rB=Degree(4);

% Alpha=[Alpha_Q1 Alpha];

Alpha_Q1=Alpha(0*(rQ+1)+1:1*(rQ+1));
Alpha_Q2=Alpha(1*(rQ+1)+1:2*(rQ+1));
Alpha_Q3=Alpha(2*(rQ+1)+1:3*(rQ+1));
Alpha_Q4=Alpha(3*(rQ+1)+1:4*(rQ+1));

Alpha_D1Q1=Alpha_Q1(1:end-1).*(rQ:-1:1);
Alpha_D1Q2=Alpha_Q2(1:end-1).*(rQ:-1:1);
Alpha_D1Q3=Alpha_Q3(1:end-1).*(rQ:-1:1);
Alpha_D1Q4=Alpha_Q4(1:end-1).*(rQ:-1:1);

Alpha_D2Q1=Alpha_D1Q1(1:end-1).*(rQ-1:-1:1);
Alpha_D2Q2=Alpha_D1Q2(1:end-1).*(rQ-1:-1:1);
Alpha_D2Q3=Alpha_D1Q3(1:end-1).*(rQ-1:-1:1);
Alpha_D2Q4=Alpha_D1Q4(1:end-1).*(rQ-1:-1:1);

Q1=polyval(Alpha_Q1,Time);
Q2=polyval(Alpha_Q2,Time);
Q3=polyval(Alpha_Q3,Time);
Q4=polyval(Alpha_Q4,Time);
QJ=[Q1;Q2;Q3;Q4];

Qhat1=zeros(size(Q1));
Qhat2=Q2+Q3;
Qhat3=Q3+Q4;
QhatJ=[Qhat1;Qhat2;Qhat3];


D1Q1=polyval(Alpha_D1Q1,Time);
D1Q2=polyval(Alpha_D1Q2,Time);
D1Q3=polyval(Alpha_D1Q3,Time);
D1Q4=polyval(Alpha_D1Q4,Time);

D2Q1=polyval(Alpha_D2Q1,Time);
D2Q2=polyval(Alpha_D2Q2,Time);
D2Q3=polyval(Alpha_D2Q3,Time);
D2Q4=polyval(Alpha_D2Q4,Time);

for i=1:length(Time)
%     q1=Q1(i);
    q2=Q2(i);
    q3=Q3(i);
    q4=Q4(i);
    D1q1=D1Q1(i);
    D1q2=D1Q2(i);
    D1q3=D1Q3(i);
    D1q4=D1Q4(i);
    D2q1=D2Q1(i);
    D2q2=D2Q2(i);
    D2q3=D2Q3(i);
    D2q4=D2Q4(i);

    
    MM= [ (63*mL1)/8000 + (LL3^2*mL3*sin(q2 + q3)^2)/4 + LL3^2*mL4*sin(q2 + q3)^2 + (LL2^2*mL2*sin(q2)^2)/4 + LL2^2*mL3*sin(q2)^2 + LL2^2*mL4*sin(q2)^2 + (LL4^2*mL4*sin(q2 + q3 + q4)^2)/4 - LL2*LL3*mL3*sin(q3/2)^2 - 2*LL2*LL3*mL4*sin(q3/2)^2 - LL3*LL4*mL4*sin(q4/2)^2 - LL2*LL4*mL4*sin(q3/2 + q4/2)^2 + LL3*LL4*mL4*sin(q2 + q3 + q4/2)^2 + LL2*LL4*mL4*sin(q2 + q3/2 + q4/2)^2 + LL2*LL3*mL3*sin(q2 + q3/2)^2 + 2*LL2*LL3*mL4*sin(q2 + q3/2)^2,    0,                        0,              0;
           0, (73*mL1)/20000 + (LL1^2*mL1)/4 + (LL2^2*mL2)/4 + LL2^2*mL3 + LL2^2*mL4 + (LL3^2*mL3)/4 + LL3^2*mL4 + (LL4^2*mL4)/4 + LL2*LL4*mL4*cos(q3 + q4) + LL2*LL3*mL3*cos(q3) + 2*LL2*LL3*mL4*cos(q3) + LL3*LL4*mL4*cos(q4), (13*mL1)/8000 + (LL1^2*mL1)/6 + (LL3^2*mL3)/4 + LL3^2*mL4 + (LL4^2*mL4)/4 + (LL2*LL4*mL4*cos(q3 + q4))/2 + (LL2*LL3*mL3*cos(q3))/2 + LL2*LL3*mL4*cos(q3) + LL3*LL4*mL4*cos(q4), mL1/2500 + (LL1^2*mL1)/12 + (LL4^2*mL4)/4 + (LL2*LL4*mL4*cos(q3 + q4))/2 + (LL3*LL4*mL4*cos(q4))/2;
           0,                                    (13*mL1)/8000 + (LL1^2*mL1)/6 + (LL3^2*mL3)/4 + LL3^2*mL4 + (LL4^2*mL4)/4 + (LL2*LL4*mL4*cos(q3 + q4))/2 + (LL2*LL3*mL3*cos(q3))/2 + LL2*LL3*mL4*cos(q3) + LL3*LL4*mL4*cos(q4),                                                                                (13*mL1)/8000 + (LL1^2*mL1)/6 + (LL3^2*mL3)/4 + LL3^2*mL4 + (LL4^2*mL4)/4 + LL3*LL4*mL4*cos(q4),                                (mL1*LL1^2)/12 + (mL4*LL4^2)/4 + (LL3*mL4*cos(q4)*LL4)/2 + mL1/2500;
           0,                                                                                                                mL1/2500 + (LL1^2*mL1)/12 + (LL4^2*mL4)/4 + (LL2*LL4*mL4*cos(q3 + q4))/2 + (LL3*LL4*mL4*cos(q4))/2,                                                                                                            (mL1*LL1^2)/12 + (mL4*LL4^2)/4 + (LL3*mL4*cos(q4)*LL4)/2 + mL1/2500,                                                          (mL1*LL1^2)/12 + (mL4*LL4^2)/4 + mL1/2500 ];
 

    CC= [ D1q4*((LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4))/8 - (LL2*LL4*mL4*sin(q3 + q4))/4 - (LL3*LL4*mL4*sin(q4))/4 + (LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4))/4 + (LL2*LL4*mL4*sin(2*q2 + q3 + q4))/4) + D1q3*((LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4))/8 + (LL3^2*mL3*sin(2*q2 + 2*q3))/8 + (LL3^2*mL4*sin(2*q2 + 2*q3))/2 - (LL2*LL4*mL4*sin(q3 + q4))/4 - (LL2*LL3*mL3*sin(q3))/4 - (LL2*LL3*mL4*sin(q3))/2 + (LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4))/2 + (LL2*LL3*mL3*sin(2*q2 + q3))/4 + (LL2*LL3*mL4*sin(2*q2 + q3))/2 + (LL2*LL4*mL4*sin(2*q2 + q3 + q4))/4) + D1q2*((LL2^2*mL2*sin(2*q2))/8 + (LL2^2*mL3*sin(2*q2))/2 + (LL2^2*mL4*sin(2*q2))/2 + (LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4))/8 + (LL3^2*mL3*sin(2*q2 + 2*q3))/8 + (LL3^2*mL4*sin(2*q2 + 2*q3))/2 + (LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4))/2 + (LL2*LL3*mL3*sin(2*q2 + q3))/2 + LL2*LL3*mL4*sin(2*q2 + q3) + (LL2*LL4*mL4*sin(2*q2 + q3 + q4))/2), (D1q1*(LL2^2*mL2*sin(2*q2) + 4*LL2^2*mL3*sin(2*q2) + 4*LL2^2*mL4*sin(2*q2) + LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4) + LL3^2*mL3*sin(2*q2 + 2*q3) + 4*LL3^2*mL4*sin(2*q2 + 2*q3) + 4*LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4) + 4*LL2*LL3*mL3*sin(2*q2 + q3) + 8*LL2*LL3*mL4*sin(2*q2 + q3) + 4*LL2*LL4*mL4*sin(2*q2 + q3 + q4)))/8, (D1q1*(LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4) + LL3^2*mL3*sin(2*q2 + 2*q3) + 4*LL3^2*mL4*sin(2*q2 + 2*q3) - 2*LL2*LL4*mL4*sin(q3 + q4) - 2*LL2*LL3*mL3*sin(q3) - 4*LL2*LL3*mL4*sin(q3) + 4*LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4) + 2*LL2*LL3*mL3*sin(2*q2 + q3) + 4*LL2*LL3*mL4*sin(2*q2 + q3) + 2*LL2*LL4*mL4*sin(2*q2 + q3 + q4)))/8, (D1q1*LL4*mL4*(2*LL2*sin(2*q2 + q3 + q4) - 2*LL2*sin(q3 + q4) - 2*LL3*sin(q4) + 2*LL3*sin(2*q2 + 2*q3 + q4) + LL4*sin(2*q2 + 2*q3 + 2*q4)))/8;
         -(D1q1*(LL2^2*mL2*sin(2*q2) + 4*LL2^2*mL3*sin(2*q2) + 4*LL2^2*mL4*sin(2*q2) + LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4) + LL3^2*mL3*sin(2*q2 + 2*q3) + 4*LL3^2*mL4*sin(2*q2 + 2*q3) + 4*LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4) + 4*LL2*LL3*mL3*sin(2*q2 + q3) + 8*LL2*LL3*mL4*sin(2*q2 + q3) + 4*LL2*LL4*mL4*sin(2*q2 + q3 + q4)))/8,                                                                                                                                                                  - D1q4*((LL2*LL4*mL4*sin(q3 + q4))/2 + (LL3*LL4*mL4*sin(q4))/2) - D1q3*((LL2*LL4*mL4*sin(q3 + q4))/2 + (LL2*LL3*mL3*sin(q3))/2 + LL2*LL3*mL4*sin(q3)),                                                                                                                   - (D1q2*LL2*(LL4*mL4*sin(q3 + q4) + LL3*mL3*sin(q3) + 2*LL3*mL4*sin(q3)))/2 - (D1q3*LL2*(LL4*mL4*sin(q3 + q4) + LL3*mL3*sin(q3) + 2*LL3*mL4*sin(q3)))/2 - (D1q4*LL4*mL4*(LL2*sin(q3 + q4) + LL3*sin(q4)))/2,                                                                            -(LL4*mL4*(LL2*sin(q3 + q4) + LL3*sin(q4))*(D1q2 + D1q3 + D1q4))/2;
         -(D1q1*(LL4^2*mL4*sin(2*q2 + 2*q3 + 2*q4) + LL3^2*mL3*sin(2*q2 + 2*q3) + 4*LL3^2*mL4*sin(2*q2 + 2*q3) - 2*LL2*LL4*mL4*sin(q3 + q4) - 2*LL2*LL3*mL3*sin(q3) - 4*LL2*LL3*mL4*sin(q3) + 4*LL3*LL4*mL4*sin(2*q2 + 2*q3 + q4) + 2*LL2*LL3*mL3*sin(2*q2 + q3) + 4*LL2*LL3*mL4*sin(2*q2 + q3) + 2*LL2*LL4*mL4*sin(2*q2 + q3 + q4)))/8,                                                                                                                                                                                                     D1q2*((LL2*LL4*mL4*sin(q3 + q4))/2 + (LL2*LL3*mL3*sin(q3))/2 + LL2*LL3*mL4*sin(q3)) - (D1q4*LL3*LL4*mL4*sin(q4))/2,                                                                                                                                                                                                                                                                                                 -(D1q4*LL3*LL4*mL4*sin(q4))/2,                                                                                                 -(LL3*LL4*mL4*sin(q4)*(D1q2 + D1q3 + D1q4))/2;
         -(D1q1*LL4*mL4*(2*LL2*sin(2*q2 + q3 + q4) - 2*LL2*sin(q3 + q4) - 2*LL3*sin(q4) + 2*LL3*sin(2*q2 + 2*q3 + q4) + LL4*sin(2*q2 + 2*q3 + 2*q4)))/8,                                                                                                                                                                                                                           D1q2*((LL2*LL4*mL4*sin(q3 + q4))/2 + (LL3*LL4*mL4*sin(q4))/2) + (D1q3*LL3*LL4*mL4*sin(q4))/2,                                                                                                                                                                                                                                                                                         (LL3*LL4*mL4*sin(q4)*(D1q2 + D1q3))/2,                                                                                                                                             0         ];
 

    GG= [ 	0;
          - g*mL4*(LL3*sin(q2 + q3) + LL2*sin(q2) + (LL4*sin(q2 + q3 + q4))/2) - g*mL3*((LL3*sin(q2 + q3))/2 + LL2*sin(q2)) - (LL2*g*mL2*sin(q2))/2;
          - g*mL4*(LL3*sin(q2 + q3) + (LL4*sin(q2 + q3 + q4))/2) - (LL3*g*mL3*sin(q2 + q3))/2;
          -(LL4*g*mL4*sin(q2 + q3 + q4))/2 ];



     TorqueDesire(:,i) = MM*[D2q1;D2q2;D2q3;D2q4] + CC*[D1q1;D1q2;D1q3;D1q4] + GG;
    
end

    
% Integral Matrix
Cost=OptimalParam_SV_4R_3D(QJ,QhatJ,TorqueDesire,nn,rU,rB,Landa,Weight,SampleRate)*Tres;

end