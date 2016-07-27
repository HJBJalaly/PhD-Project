function Pos=FK_RzRyRyRy_3D(q1,q2,q3,q4,d1,a2,a3,a4)

LL1=d1;
LL2=a2;
LL3=a3;
LL4=a4;

Pos=[];
for i=1:length(q1)
    Pos(:,i)=[  LL3*(cos(q1(i))*cos(q2(i))*sin(q3(i)) + cos(q1(i))*cos(q3(i))*sin(q2(i))) + LL4*(cos(q4(i))*(cos(q1(i))*cos(q2(i))*sin(q3(i)) + cos(q1(i))*cos(q3(i))*sin(q2(i))) + sin(q4(i))*(cos(q1(i))*cos(q2(i))*cos(q3(i)) - cos(q1(i))*sin(q2(i))*sin(q3(i)))) + LL2*cos(q1(i))*sin(q2(i));
                LL3*(cos(q2(i))*sin(q1(i))*sin(q3(i)) + cos(q3(i))*sin(q1(i))*sin(q2(i))) + LL4*(cos(q4(i))*(cos(q2(i))*sin(q1(i))*sin(q3(i)) + cos(q3(i))*sin(q1(i))*sin(q2(i))) - sin(q4(i))*(sin(q1(i))*sin(q2(i))*sin(q3(i)) - cos(q2(i))*cos(q3(i))*sin(q1(i)))) + LL2*sin(q1(i))*sin(q2(i));
                LL1 + LL3*(cos(q2(i))*cos(q3(i)) - sin(q2(i))*sin(q3(i))) + LL2*cos(q2(i)) + LL4*(cos(q4(i))*(cos(q2(i))*cos(q3(i)) - sin(q2(i))*sin(q3(i))) - sin(q4(i))*(cos(q2(i))*sin(q3(i)) + cos(q3(i))*sin(q2(i))))];
end
% A1=DHP(q1(i),d1,0,pi/2);
% A2=DHP(q2(i)+pi/2,0 ,a2,0);
% A3=DHP(q3(i),0 ,a3,0);
% A4=DHP(q4(i),0 ,a4,0);
% 
% H_4_0=A1*A2*A3*A4;
% 
% v1=[0;0;0;1];
% pos=H_4_0*v1;
% Pos=pos(1:3);
end

function A=DHP(teta,d,a,alfa)
    A=[cos(teta) -sin(teta)*cos(alfa) sin(teta)*sin(alfa) a*cos(teta);...
      sin(teta) cos(teta)*cos(alfa) -cos(teta)*sin(alfa) a*sin(teta);...
      0 sin(alfa) cos(alfa) d;...
      0 0 0 1];
end
