function [Cneq,Ceq]=CF3b_NonLinearConstraint_WithJointLimit_4R_3D(Alpha,Time,Tres,Degree,Xef,Yef,Zef,g,mL1,mL2,mL3,mL4,LL1,LL2,LL3,LL4,Q_limit,DQ_limit,SampleRate)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Pos=[Xef;Yef;Zef];

rQ=Degree(2);

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
Qq=[Q1;Q2;Q3;Q4];

D1Q1=polyval(Alpha_D1Q1,Time);
D1Q2=polyval(Alpha_D1Q2,Time);
D1Q3=polyval(Alpha_D1Q3,Time);
D1Q4=polyval(Alpha_D1Q4,Time);
D1Qq=[D1Q1;D1Q2;D1Q3;D1Q4];

D2Q1=polyval(Alpha_D2Q1,Time);
D2Q2=polyval(Alpha_D2Q2,Time);
D2Q3=polyval(Alpha_D2Q3,Time);
D2Q4=polyval(Alpha_D2Q4,Time);
D2Qq=[D2Q1;D2Q2;D2Q3;D2Q4];


RPos=FK_RzRyRyRy_3D(Q1,Q2,Q3,Q4,LL1,LL2,LL3,LL4);
 
TorqueDesire=TorqueCalculator_4R_3D(D2Qq,D1Qq,Qq,g,mL1,mL2,mL3,mL4,LL1,LL2,LL3,LL4);    
    
    
Cneq=[sum(sum((RPos-Pos).^2))*Tres-.00000005;
      InRangeShifter(Q1(1:SampleRate:end))'-Q_limit(1,2);
     -InRangeShifter(Q1(1:SampleRate:end))'+Q_limit(1,1);
      InRangeShifter(Q2(1:SampleRate:end))'-Q_limit(2,2);
     -InRangeShifter(Q2(1:SampleRate:end))'+Q_limit(2,1);
      InRangeShifter(Q3(1:SampleRate:end))'-Q_limit(3,2);
     -InRangeShifter(Q3(1:SampleRate:end))'+Q_limit(3,1);
      InRangeShifter(Q4(1:SampleRate:end))'-Q_limit(4,2);
     -InRangeShifter(Q4(1:SampleRate:end))'+Q_limit(4,1);
      D1Q1(1:SampleRate:end)'-DQ_limit(1,2);
     -D1Q1(1:SampleRate:end)'+DQ_limit(1,1);
      D1Q2(1:SampleRate:end)'-DQ_limit(2,2);
     -D1Q2(1:SampleRate:end)'+DQ_limit(2,1);
      D1Q3(1:SampleRate:end)'-DQ_limit(3,2);
     -D1Q3(1:SampleRate:end)'+DQ_limit(3,1);
      D1Q4(1:SampleRate:end)'-DQ_limit(4,2);
     -D1Q4(1:SampleRate:end)'+DQ_limit(4,1)];
     

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
     Q4(1)-Q4(end);
     D1Q1(1)-D1Q1(end);
     D1Q2(1)-D1Q2(end);
     D1Q3(1)-D1Q3(end);
     D1Q4(1)-D1Q4(end);
     D2Q1(1)-D2Q1(end);
     D2Q2(1)-D2Q2(end);
     D2Q3(1)-D2Q3(end);
     D2Q4(1)-D2Q4(end);
    (diff(TorqueDesire([1],1:2)')'-diff(TorqueDesire([1],end-1:end)')')*(1e-100)];
end