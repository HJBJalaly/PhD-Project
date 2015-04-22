function Torque=TorqueCalculator(t,Q,g,mL1,mL2,mL3,LL1,LL2,LL3,Pos,DPos,D2Pos,Time,Kp,Kd)
    

% Torque=zeros(3,length(Q));


for i=1:size(Q,1)
    q1=Q(i,1);
    q2=Q(i,2);
    q3=Q(i,3);

    D1q1=Q(i,4);
    D1q2=Q(i,5);
    D1q3=Q(i,6);
   
    D1q=[D1q1 D1q2 D1q3]';


    % Dynamic
    MM=[mL3*LL1*LL3*cos(q2+q3)+LL1*LL2*cos(q2)*mL2+2*LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL3*LL1^2+LL1^2*mL2+mL1*LL1^2/3+mL2*LL2^2/3,mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6;
        mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6;
        mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6,mL3*LL3^2/3];


    CC= [-mL3*LL1*LL3*sin(q2+q3)*(D1q2+D1q3)/2-LL2*(LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL1*LL3*(D1q1+D1q2+D1q3)*sin(q2+q3)/2-LL2*(sin(q2)*D1q1*LL1*(mL2+2*mL3)+LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL3*(LL1*sin(q2+q3)+LL2*sin(q3))*(D1q1+D1q2+D1q3)/2;
          mL3*LL1*D1q1*LL3*sin(q2+q3)/2+(sin(q2)*D1q1*LL1*(mL2+2*mL3)-mL3*LL3*sin(q3)*D1q3)*LL2/2,-mL3*LL2*LL3*sin(q3)*D1q3/2,-mL3*LL3*sin(q3)*(D1q1+D1q2+D1q3)*LL2/2;
          mL3*(LL1*D1q1*sin(q2+q3)+LL2*sin(q3)*(D1q1+D1q2))*LL3/2,LL3*LL2*sin(q3)*(D1q1+D1q2)*mL3/2,0];


    GG= [g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2+2*cos(q1)*LL1*mL3+2*cos(q1)*LL1*mL2+cos(q1)*LL1*mL1)/2;
         g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2)/2;
         g*mL3*LL3*cos(q1+q2+q3)/2];

    JJ=[-LL1*sin(q1)-LL2*sin(q1+q2)-LL3*sin(q1+q2+q3), -LL2*sin(q1+q2)-LL3*sin(q1+q2+q3), -LL3*sin(q1+q2+q3);
          LL1*cos(q1)+LL2*cos(q1+q2)+LL3*cos(q1+q2+q3),  LL2*cos(q1+q2)+LL3*cos(q1+q2+q3),  LL3*cos(q1+q2+q3)]; 

    X=[LL1*cos(q1)+LL2*cos(q1+q2)+LL3*cos(q1+q2+q3);
       LL1*sin(q1)+LL2*sin(q1+q2)+LL3*sin(q1+q2+q3)];


    Xd=interp1(Time,Pos',t(i),'spline')';
    D1Xd=interp1(Time,DPos',t(i),'spline')';
    D2Xd=interp1(Time,D2Pos',t(i),'spline')';

    D2xr=D2Xd+Kd*(D1Xd-JJ*D1q)+Kp*(Xd-X);

    Jsharp=JJ'*(JJ*JJ')^-1;
    DJ=[-LL1*cos(q1)*D1q1-LL2*cos(q1+q2)*(D1q1+D1q2)-LL3*cos(q1+q2+q3)*(D1q1+D1q2+D1q3)  , -LL2*cos(q1+q2)*(D1q1+D1q2)-LL3*cos(q1+q2+q3)*(D1q1+D1q2+D1q3) , -LL3*cos(q1+q2+q3)*(D1q1+D1q2+D1q3);
        -LL1*sin(q1)*D1q1-LL2*sin(q1+q2)*(D1q1+D1q2)-LL3*sin(q1+q2+q3)*(D1q1+D1q2+D1q3)  , -LL2*sin(q1+q2)*(D1q1+D1q2)-LL3*sin(q1+q2+q3)*(D1q1+D1q2+D1q3) , -LL3*sin(q1+q2+q3)*(D1q1+D1q2+D1q3)]; 
   
    D2qr=Jsharp*(D2xr-DJ*D1q);
    Torque(:,i)=MM*D2qr+CC*D1q+GG;
    
end
1;