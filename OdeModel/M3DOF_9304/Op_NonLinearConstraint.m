function [Cneq,Ceq]=NonLinearConstraintFun1( Coef,time,Tres,Degree,L,Xef,Yef)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Pos=[Xef;Yef];


CoefP_q1=Coef(1:(Degree+1));
CoefP_q2=Coef((Degree+1)+1:2*(Degree+1));
CoefP_q3=Coef(2*(Degree+1)+1:3*(Degree+1));
q1=polyval(CoefP_q1,time);
q2=polyval(CoefP_q2,time);
q3=polyval(CoefP_q3,time);

D1CoefP_q1=CoefP_q1(1:end-1).*(Degree:-1:1);
D1CoefP_q2=CoefP_q2(1:end-1).*(Degree:-1:1);
D1CoefP_q3=CoefP_q3(1:end-1).*(Degree:-1:1);

D2CoefP_q1=D1CoefP_q1(1:end-1).*(Degree-1:-1:1);
D2CoefP_q2=D1CoefP_q2(1:end-1).*(Degree-1:-1:1);
D2CoefP_q3=D1CoefP_q3(1:end-1).*(Degree-1:-1:1);

D1Q1=polyval(D1CoefP_q1,time([1,end]));
D1Q2=polyval(D1CoefP_q2,time([1,end]));
D1Q3=polyval(D1CoefP_q3,time([1,end]));

D2Q1=polyval(D2CoefP_q1,time([1,end]));
D2Q2=polyval(D2CoefP_q2,time([1,end]));
D2Q3=polyval(D2CoefP_q3,time([1,end]));


RPos=L*[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
        sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];

    
Cneq=sum(sum((RPos-Pos).^2))*Tres-.0002;
Ceq=[RPos(:,1)-Pos(:,1);
%      RPos(:,end)-Pos(:,end);
     q1(1)-q1(end);
     q2(1)-q2(end);
     q3(1)-q3(end);
     D1Q1(1)-D1Q1(end);
     D1Q2(1)-D1Q2(end);
     D1Q3(1)-D1Q3(end);
     D2Q1(1)-D2Q1(end);
     D2Q2(1)-D2Q2(end);
     D2Q3(1)-D2Q3(end)];    
end

