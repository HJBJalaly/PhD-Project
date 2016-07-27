

%% Ode in M, D and C Matrix for my robot for 4 Dof Kuka (R1100 fivve):

home 
clear
syms mL1 mL2 mL3 mL4 LL1 LL2 LL3 LL4 g real positive
syms q1 q2 q3 q4 real
syms D1q1 D1q2 D1q3 D1q4 real
syms D2q1 D2q2 D2q3 D2q4 real

II1=1/12*mL1*diag([3*(.13)^2+LL1^2, 3*(.13)^2+LL1^2, 3*(.13)^2]);
II2=1/12*mL1*diag([3*(.09)^2+LL1^2, 3*(.09)^2+LL1^2, 3*(.09)^2]);
II3=1/12*mL1*diag([3*(.07)^2+LL1^2, 3*(.07)^2+LL1^2, 3*(.07)^2]);
II4=1/12*mL1*diag([3*(.04)^2+LL1^2, 3*(.04)^2+LL1^2, 3*(.04)^2]);

QQ=[q1 q2 q3 q4]';
DQ=[D1q1 D1q2 D1q3 D1q4]';
D2Q=[D2q1 D2q2 D2q3 D2q4]';
%%

i0=[1 0 0];             
j0=[0 1 0];
k0=[0 0 1];

R1=[ cos(q1) sin(q1) 0;
    -sin(q1) cos(q1) 0;
     0       0       1];
Unit1=R1*[i0;j0;k0];       % unit vector of coordinate 1
i1=Unit1(1,:);
j1=Unit1(2,:);
k1=Unit1(3,:);
 
R2=[ cos(q2) 0 -sin(q2);
      0       1  0;
      sin(q2) 0  cos(q2)];
Unit2=R2*[i1;j1;k1];       % unit vector of coordinate 2
i2=Unit2(1,:);
j2=Unit2(2,:);
k2=Unit2(3,:);
 
R3=[ cos(q3) 0 -sin(q3);
     0       1  0;
     sin(q3) 0  cos(q3)];
Unit3=R3*[i2;j2;k2];       % unit vector of coordinate 3
i3=Unit3(1,:);
j3=Unit3(2,:);
k3=Unit3(3,:);

R4=[ cos(q4) 0 -sin(q4);
     0       1  0;
     sin(q4) 0  cos(q4)];
Unit4=R4*[i3;j3;k3];       % unit vector of coordinate 4
i4=Unit4(1,:);
j4=Unit4(2,:);
k4=Unit4(3,:);

Pc1=0;              % orgin of coordinate of Link 1
Pg1=Pc1+LL1/2*k1;   % position of center of mass of Link 1   
w1=D1q1*k1;         % angular velocity of Link 1
Vc1=0;                      % velocity of orgin of coordinate of Link 1
Vg1=Vc1+cross(w1,LL1/2*k1); % velocity of center of mass of Link 1   

Pc2=Pc1+LL1*k1;     % orgin of coordinate of Link 2
Pg2=Pc2+LL2/2*k2;   % position of center of mass of Link 2
w2=w1+D1q2*j2;         % angular velocity of Link 2
Vc2=Vc1+cross(w1,LL1*k1); % velocity of orgin of coordinate of Link 2
Vg2=Vc2+cross(w2,LL2/2*k2); % velocity of center of mass of Link 2

Pc3=Pc2+LL2*k2;     % orgin of coordinate of Link 3
Pg3=Pc3+LL3/2*k3;   % position of center of mass of Link 3
w3=w2+D1q3*j3;         % angular velocity of Link 3
Vc3=Vc2+cross(w2,LL2*k2); % velocity of orgin of coordinate of Link 3
Vg3=Vc3+cross(w3,LL3/2*k3); % velocity of center of mass of Link 3

Pc4=Pc3+LL3*k3;     % orgin of coordinate of Link 4
Pg4=Pc4+LL4/2*k4;   % position of center of mass of Link 4
w4=w3+D1q4*j4;         % angular velocity of Link 4
Vc4=Vc3+cross(w3,LL3*k3); % velocity of orgin of coordinate of Link 4
Vg4=Vc4+cross(w4,LL4/2*k4); % velocity of center of mass of Link 4

PEE=Pc4+LL4*k4;     % EE

%% Lagrange

PotEnergy=g*(mL1*Pg1(3)+mL2*Pg2(3)+mL3*Pg3(3)+mL4*Pg4(3));
PotEnergy = simple(PotEnergy);

KinEnergy=1/2*mL1*Vg1*Vg1'+1/2*mL2*Vg2*Vg2'+1/2*mL3*Vg3*Vg3'+1/2*mL4*Vg4*Vg4'+...
            1/2*w1*II1*w1'+1/2*w2*II2*w2'+1/2*w3*II3*w3'+1/2*w4*II4*w4';
KinEnergy=simple(KinEnergy);
%%

GG = jacobian(PotEnergy,QQ).';
GG = simple(GG);

% mass-inertial matrix
DD = simple(jacobian(KinEnergy,DQ).');
DD = simple(jacobian(DD,DQ));

% Coriolis and centrifugal matrix
syms CC real
n=max(size(QQ));
for k=1:n
	for j=1:n
		CC(k,j)=0*g;
		for i=1:n
			CC(k,j)=CC(k,j)+1/2*(diff(DD(k,j),QQ(i)) + ...
				diff(DD(k,i),QQ(j)) - ...
				diff(DD(i,j),QQ(k)))*DQ(i);
		end
	end
end
CC=simple(CC);

'fin'
%%

Jacobian=jacobian(PEE,QQ)
