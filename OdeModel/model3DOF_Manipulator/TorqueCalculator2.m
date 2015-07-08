function Torque=TorqueCalculator2(t,Q,g,mL1,mL2,LL1,LL2,Pos,DPos,D2Pos,Time,Kp,Kd)
    

% Torque=zeros(3,length(Q));


for i=1:size(Q,1)
    q1=Q(i,1);
    q2=Q(i,2);

    D1q1=Q(i,3);
    D1q2=Q(i,4);
   
    D1q=[D1q1 D1q2]';


    % Dynamic
    MM=[mL2*LL1*LL2*cos(q2)+mL2*LL1^2+mL2*LL2^2/3+mL1*LL1^2/3, mL2*LL2*(3*LL1*cos(q2)+2*LL2)/6;
        mL2*LL2*(3*LL1*cos(q2)+2*LL2)/6, mL2*LL2^2/3];


    CC= [-mL2*LL2*sin(q2)*D1q2*LL1/2, -mL2*LL2*sin(q2)*(D1q1+D1q2)*LL1/2; 
          mL2*LL2*LL1*D1q1*sin(q2)/2, 0];


    GG= [g*(mL2*LL2*cos(q1+q2)+LL1*cos(q1)*mL1+2*LL1*cos(q1)*mL2)/2;
         g*mL2*LL2*cos(q1+q2)/2];

%  Control Law
    JJ=[-LL1*sin(q1)-LL2*sin(q1+q2), -LL2*sin(q1+q2);
         LL1*cos(q1)+LL2*cos(q1+q2),  LL2*cos(q1+q2)]; 

    X=[LL1*cos(q1)+LL2*cos(q1+q2);
       LL1*sin(q1)+LL2*sin(q1+q2)];
        
    Xd=interp1(Time,Pos',t(i),'spline')';
    D1Xd=interp1(Time,DPos',t(i),'spline')';
    D2Xd=interp1(Time,D2Pos',t(i),'spline')';

    D2xr=D2Xd+Kd*(D1Xd-JJ*D1q)+Kp*(Xd-X);

    Jsharp=JJ'*(JJ*JJ')^-1;
    
    DJ=[-LL1*cos(q1)*D1q1-LL2*cos(q1+q2)*(D1q1+D1q2) , -LL2*cos(q1+q2)*(D1q1+D1q2);
        -LL1*sin(q1)*D1q1-LL2*sin(q1+q2)*(D1q1+D1q2), -LL2*sin(q1+q2)*(D1q1+D1q2)]; 


   
    D2qr=Jsharp*(D2xr-DJ*D1q);
    
    Torque(:,i)=MM*D2qr+CC*D1q+GG;
    
end
1;