function [Cneq,Ceq]=CF3b_NonLinearConstraint(Alpha,Time,Tres,Degree,L,Xef,Yef,g,mL1,mL2,mL3,LL1,LL2,LL3)
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
 
TorqueDesire=[];    
for i=[1,2,length(Time)-1,length(Time)]
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


     TorqueDesire(:,end+1) = MM*[D2q1;D2q2;D2q3] + CC*[D1q1;D1q2;D1q3] + GG;
    
end
    
    
    
Cneq=sum(sum((RPos-Pos).^2))*Tres-.000005;
Middle1=ceil(1*length(Time)/8);
Middle2=ceil(2*length(Time)/8);
Middle3=ceil(3*length(Time)/8);
Middle4=ceil(4*length(Time)/8);
Middle5=ceil(5*length(Time)/8);
Middle6=ceil(6*length(Time)/8);
Middle7=ceil(7*length(Time)/8);
Ceq=[RPos(:,1)-Pos( :,1);
    (RPos(:,Middle1)-Pos( :,Middle1))*1;
    RPos(:,Middle2)-Pos( :,Middle2);
    (RPos(:,Middle3)-Pos( :,Middle3))*.1;
    RPos(:,Middle4)-Pos( :,Middle4);
    RPos(:,Middle5)-Pos( :,Middle5);
    RPos(:,Middle6)-Pos( :,Middle6);
    RPos(:,Middle7)-Pos( :,Middle7);
%      RPos(:,end)-Pos(:,end);
     Q1(1)-Q1(end);
     Q2(1)-Q2(end);
     Q3(1)-Q3(end);
     D1Q1(1)-D1Q1(end);
     D1Q2(1)-D1Q2(end);
     D1Q3(1)-D1Q3(end);
     D2Q1(1)-D2Q1(end);
     D2Q2(1)-D2Q2(end);
     D2Q3(1)-D2Q3(end);
     (diff(TorqueDesire([1],1:2)')'-diff(TorqueDesire([1],end-1:end)')')*(1e-100)];
end