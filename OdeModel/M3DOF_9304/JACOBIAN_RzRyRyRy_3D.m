function J = JACOBIAN_RzRyRyRy_3D(q1,q2,q3,q4,d1,a2,a3,a4)

LL1=d1;
LL2=a2;
LL3=a3;
LL4=a4;

J=[ - LL3*(cos(q2)*sin(q1)*sin(q3) + cos(q3)*sin(q1)*sin(q2)) - LL4*(cos(q4)*(cos(q2)*sin(q1)*sin(q3) + cos(q3)*sin(q1)*sin(q2)) - sin(q4)*(sin(q1)*sin(q2)*sin(q3) - cos(q2)*cos(q3)*sin(q1))) - LL2*sin(q1)*sin(q2), LL3*(cos(q1)*cos(q2)*cos(q3) - cos(q1)*sin(q2)*sin(q3)) + LL4*(cos(q4)*(cos(q1)*cos(q2)*cos(q3) - cos(q1)*sin(q2)*sin(q3)) - sin(q4)*(cos(q1)*cos(q2)*sin(q3) + cos(q1)*cos(q3)*sin(q2))) + LL2*cos(q1)*cos(q2),   LL3*(cos(q1)*cos(q2)*cos(q3) - cos(q1)*sin(q2)*sin(q3)) + LL4*(cos(q4)*(cos(q1)*cos(q2)*cos(q3) - cos(q1)*sin(q2)*sin(q3)) - sin(q4)*(cos(q1)*cos(q2)*sin(q3) + cos(q1)*cos(q3)*sin(q2))),  LL4*(cos(q4)*(cos(q1)*cos(q2)*cos(q3) - cos(q1)*sin(q2)*sin(q3)) - sin(q4)*(cos(q1)*cos(q2)*sin(q3) + cos(q1)*cos(q3)*sin(q2)));
      LL3*(cos(q1)*cos(q2)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + LL4*(cos(q4)*(cos(q1)*cos(q2)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + sin(q4)*(cos(q1)*cos(q2)*cos(q3) - cos(q1)*sin(q2)*sin(q3))) + LL2*cos(q1)*sin(q2), LL2*cos(q2)*sin(q1) - LL4*(cos(q4)*(sin(q1)*sin(q2)*sin(q3) - cos(q2)*cos(q3)*sin(q1)) + sin(q4)*(cos(q2)*sin(q1)*sin(q3) + cos(q3)*sin(q1)*sin(q2))) - LL3*(sin(q1)*sin(q2)*sin(q3) - cos(q2)*cos(q3)*sin(q1)), - LL3*(sin(q1)*sin(q2)*sin(q3) - cos(q2)*cos(q3)*sin(q1)) - LL4*(cos(q4)*(sin(q1)*sin(q2)*sin(q3) - cos(q2)*cos(q3)*sin(q1)) + sin(q4)*(cos(q2)*sin(q1)*sin(q3) + cos(q3)*sin(q1)*sin(q2))), -LL4*(cos(q4)*(sin(q1)*sin(q2)*sin(q3) - cos(q2)*cos(q3)*sin(q1)) + sin(q4)*(cos(q2)*sin(q1)*sin(q3) + cos(q3)*sin(q1)*sin(q2)));
                                                                                                                                                                                                                   0,                                                       - LL3*(cos(q2)*sin(q3) + cos(q3)*sin(q2)) - LL2*sin(q2) - LL4*(cos(q4)*(cos(q2)*sin(q3) + cos(q3)*sin(q2)) + sin(q4)*(cos(q2)*cos(q3) - sin(q2)*sin(q3))),                                                 - LL3*(cos(q2)*sin(q3) + cos(q3)*sin(q2)) - LL4*(cos(q4)*(cos(q2)*sin(q3) + cos(q3)*sin(q2)) + sin(q4)*(cos(q2)*cos(q3) - sin(q2)*sin(q3))),                                 -LL4*(cos(q4)*(cos(q2)*sin(q3) + cos(q3)*sin(q2)) + sin(q4)*(cos(q2)*cos(q3) - sin(q2)*sin(q3)))  ];
 

% A1=DHP(q1,d1,0,pi/2);
% A2=DHP(q2+pi/2,0 ,a2,0);
% A3=DHP(q3,0 ,a3,0);
% A4=DHP(q4,0 ,a4,0);
% 
% H_1_0=A1;
% H_2_0=A1*A2;
% H_3_0=A1*A2*A3;
% H_4_0=A1*A2*A3*A4;
% 
% z0=[0;0;1];
% Vv=[0;0;0;1];
% 
% L0=H_4_0*Vv;
% L1=H_4_0*Vv-H_1_0*Vv;
% L2=H_4_0*Vv-H_2_0*Vv;
% L3=H_4_0*Vv-H_3_0*Vv;
% 
% J=[cross(z0,L0(1:3)), cross(H_1_0(1:3,1:3)*z0,L1(1:3))  , cross(H_2_0(1:3,1:3)*z0,L2(1:3)), cross(H_3_0(1:3,1:3)*z0,L3(1:3)) ];
end

function A=DHP(teta,d,a,alfa)
    A=[cos(teta) -sin(teta)*cos(alfa) sin(teta)*sin(alfa) a*cos(teta);...
      sin(teta) cos(teta)*cos(alfa) -cos(teta)*sin(alfa) a*sin(teta);...
      0 sin(alfa) cos(alfa) d;...
      0 0 0 1];
end
