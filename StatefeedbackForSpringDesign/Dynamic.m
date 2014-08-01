syms m1 m2 B1 B2
syms k11 k12 k21 k22

m1=1;
m2=2;
B1=0.5;
B2=0.75;

% M.D2r + C.D1r + G = F
% F=[ u1;u2]
M=[m1+m2 , m2 ; 0,  m2];
C=[B1 0; 0 B2];
G=[0 0;0 0];

% x=Tr
% Dx= Ax+Bu
A = [zeros(2)   eye(2)  ; zeros(2)  -M^-1 * C];
B=[zeros(2);eye(2)];

Poles=eig(A)

K=[k11 k12 0 0 ; k21 k22 0  0]

Ahat=A-B*K

% m1=1;
% m2=1;
% B1=0.5;
% B2=0.5;
Ahats=subs(Ahat,[m1 m2 B1 B2],[1,1,0.5,0.5])

A