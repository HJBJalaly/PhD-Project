function [Cneq,Ceq]=CF2_NonLinearConstraint(Alpha,Time,Tres,Degree,L,Xef,Yef)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Pos=[Xef;Yef];

rQ=Degree(2);

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
% QVal=[Q1;Q2;Q3];

D1Q1=polyval(Alpha_D1Q1,Time);
D1Q2=polyval(Alpha_D1Q2,Time);
D1Q3=polyval(Alpha_D1Q3,Time);

D2Q1=polyval(Alpha_D2Q1,Time);
D2Q2=polyval(Alpha_D2Q2,Time);
D2Q3=polyval(Alpha_D2Q3,Time);



RPos=L*[cos(Q1)+cos(Q1+Q2)+cos(Q1+Q2+Q3);
        sin(Q1)+sin(Q1+Q2)+sin(Q1+Q2+Q3)];

    
Cneq=sum(sum((RPos-Pos).^2))*Tres-.0001;
Ceq=[RPos(:,1)-Pos(:,1);
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