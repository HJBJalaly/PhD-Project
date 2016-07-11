% A MATLAB script to generate the equations of motion for a 5-link,
% planar, kneed, biped walker using the method of Lagrange.
%
% This routine is for a robot with an upright torso, two legs and
% knees.  The model is for five degrees of freedom of the robot,
% with the stance foot (Foot1) pinned to the ground.  The notation
% used for the equation of motions are as in Robot Dynamics and
% Control by Spong and Vidyasagar (1989), page 142, Eq. (6.3.12)
%
%    D(q)ddq + C(q,dq)*dq + G(q) = B*tau + E2*F_external
%
% Note the following convention
%
%  fem = femur, tib = tibia, 1 = stance leg, 2 = swing leg
%
% For simplicity, the model is derived using absolute coorinates
% measured in the trigonometric sense from the vertical.  The model
% is trasnformed in to relative coordates plus one absolute
% coordinates via a coordinate transformation.
%
% This file is associated with the book Feedback Control of Dynamic 
% Bipedal Robot Locomotion by Eric R. Westervelt, Jessy W. Grizzle, 
% Christine Chevallereau, Jun-Ho Choi, and Benjamin Morris published 
% by Taylor & Francis/CRC Press in 2007.
% 
% Copyright (c) 2007 by Eric R. Westervelt, Jessy W. Grizzle, Christine
% Chevallereau, Jun-Ho Choi, and Benjamin Morris.  This code may be
% freely used for noncommercial ends.  If use of this code in part or in
% whole results in publication, proper citation must be included in that
% publication.  This code comes with no guarantees or support.
% 
% Eric Westervelt
% 20 February 2007

clear all

% -----------------------------------------------------------------
%
%  Model variables
%
% -----------------------------------------------------------------

% absolute joint angles and velocities
syms q_torso dq_torso real
syms q_fem1 dq_fem1 real
syms q_fem2 dq_fem2 real
syms q_tib1 dq_tib1 real
syms q_tib2 dq_tib2 real

% gravity
syms g real

% link lengths
syms L_torso real
syms L_fem real
syms L_tib real

% link masses
syms M_torso real
syms M_fem real
syms M_tib real

% center of mass (distance form top)
syms Lc_torso real
syms Lc_fem real
syms Lc_tib real

% link inertias
syms XX_torso real
syms XX_fem real
syms XX_tib real
1
%%
% relative joint angles
syms q1 q2 q3 q4 q5 real
qs  = [q1; q2; q3; q4; q5];
syms x1 y1 real % position of stance leg
qe=[qs;x1;y1];

% relative joint velocities
syms dq1 dq2 dq3 dq4 dq5 real
dqs = [dq1; dq2; dq3; dq4; dq5];
syms dx1 dy1 real % velocity of stance leg
dqe=[dqs;dx1;dy1];

% relative joint accelaration
syms d2q1 d2q2 d2q3 d2q4 d2q5 real
d2qs = [d2q1; d2q2; d2q3; d2q4; d2q5];
syms d2x1 d2y1 real % accelaration of stance leg
d2qe=[d2qs;d2x1;d2y1];

% -----------------------------------------------------------------
%
%  Change coordinates to relative coordinates
%
% -----------------------------------------------------------------

T = [ 1   0   0   0  -1;
      0   1   0   0  -1;
     -1   0   1   0   0;
	  0  -1   0   1   0;
	  0   0   0   0   1];% q_abs= inv(T)* q_rel

% old absolute coordinates in terms of relative coords
q_new  = inv(T)*qs;
dq_new = inv(T)*dqs;

q_fem1  = q_new(1);
q_fem2  = q_new(2);
q_tib1  = q_new(3);
q_tib2  = q_new(4);
q_torso = q_new(5);
dq_fem1  = dq_new(1);
dq_fem2  = dq_new(2);
dq_tib1  = dq_new(3);
dq_tib2  = dq_new(4);
dq_torso = dq_new(5);
2
%%
% -----------------------------------------------------------------
%
%  Position and velocity
%
% -----------------------------------------------------------------

% hip, knee and center of mass positions
p_tib1     =[x1;y1];
p_cen_tib1 =p_tib1 + (L_tib-Lc_tib    )*[ sin(pi-q_tib1); cos(pi-q_tib1)];
p_knee1    =p_tib1 + (L_tib           )*[ sin(pi-q_tib1); cos(pi-q_tib1)];
p_cen_fem1 =p_knee1+ (L_fem-Lc_fem    )*[ sin(pi-q_fem1); cos(pi-q_fem1)];
p_hip      =p_knee1+ (L_fem           )*[ sin(pi-q_fem1); cos(pi-q_fem1)];
p_cen_torso=p_hip  + (L_torso-Lc_torso)*[-sin(q_torso)  ; cos(q_torso)  ];
p_cen_fem2 =p_hip  + (Lc_fem          )*[ sin(q_fem2-pi);-cos(q_fem2-pi)];
p_knee2    =p_hip  + (L_fem           )*[ sin(q_fem2-pi);-cos(q_fem2-pi)];
p_cen_tib2 =p_knee2+ (Lc_tib          )*[ sin(q_tib2-pi);-cos(pi-q_tib2)];
p_tib2     =p_knee2+ (L_tib           )*[ sin(q_tib2-pi);-cos(pi-q_tib2)];
p_com      = (M_tib*p_cen_tib1 + M_fem*p_cen_fem1 + M_torso*p_cen_torso + M_fem*p_cen_fem2+M_tib*p_cen_tib2)/(2*M_fem+2*M_tib+M_torso);


% hip and knee velocities
v_tib1     =jacobian(p_tib1     ,qe)*dqe;
v_cen_tib1 =jacobian(p_cen_tib1 ,qe)*dqe;
v_knee1    =jacobian(p_knee1    ,qe)*dqe;
v_cen_fem1 =jacobian(p_cen_fem1 ,qe)*dqe;
v_hip      =jacobian(p_hip      ,qe)*dqe;
v_cen_torso=jacobian(p_cen_torso,qe)*dqe;
v_cen_fem2 =jacobian(p_cen_fem2 ,qe)*dqe;
v_knee2    =jacobian(p_knee2    ,qe)*dqe;
v_cen_tib2 =jacobian(p_cen_tib2 ,qe)*dqe;
v_tib2     =jacobian(p_tib2     ,qe)*dqe;
v_com     =jacobian(p_com     ,qe)*dqe;

a_com(1)  =dqe'*jacobian( jacobian(p_com(1)     ,qe),qe)*dqe + jacobian(p_com(1)     ,qe)*d2qe;
a_com(2)  =dqe'*jacobian( jacobian(p_com(2)     ,qe),qe)*dqe + jacobian(p_com(2)     ,qe)*d2qe;
% a_com(1)*m = GRF_x
% a_com(1)*m +m*g = GRF_y


%angular velocties
w_torso = dq_torso;
w_fem1  = dq_fem1;
w_fem2  = dq_fem2;
w_tib1  = dq_tib1;
w_tib2  = dq_tib2;
3
%% Construct Lagrange eqaution
% -----------------------------------------------------------------
%
%  Calculate kinetic energy
%
% -----------------------------------------------------------------

% kinetic energy of links
KE_torso = 1/2*M_torso*(v_cen_torso'*v_cen_torso) + 1/2*XX_torso*(w_torso)^2;
KE_torso =simplify(KE_torso);
KE_fem1  = 1/2*M_fem*(v_cen_fem1'*v_cen_fem1)     + 1/2*XX_fem*(w_fem1)^2;
KE_fem1  =simplify(KE_fem1);
KE_fem2  = 1/2*M_fem*(v_cen_fem2'*v_cen_fem2)     + 1/2*XX_fem*(w_fem2)^2;
KE_fem2  =simplify(KE_fem2);
KE_tib1  = 1/2*M_tib*(v_cen_tib1'*v_cen_tib1)     + 1/2*XX_tib*(w_tib1)^2;
KE_tib1  =simplify(KE_tib1);
KE_tib2 = 1/2*M_tib*(v_cen_tib2'*v_cen_tib2)      + 1/2*XX_tib*(w_tib2)^2;
KE_tib2 = simplify(KE_tib2);

% total kinetic energy
KE = KE_torso + KE_fem1 + KE_fem2 + KE_tib1 + KE_tib2;
KE = simple(KE);

% -----------------------------------------------------------------
%
%  Calculate potential energy
%
% -----------------------------------------------------------------

% total potential energy
PE = g*(M_torso*(p_cen_torso(2) +p_tib1(2))+...
        M_fem*(p_cen_fem1(2)+p_tib1(2)) + M_fem*(p_cen_fem2(2) +p_tib1(2))+ ...
	    M_tib*(p_cen_tib1(2)+p_tib1(2)) + M_tib*(p_cen_tib2(2) +p_tib1(2)) );
PE = simple(PE);
4
%%
% -----------------------------------------------------------------
%
%  Calculate model matrices
%
% -----------------------------------------------------------------

% gravity vector
G_vect = jacobian(PE,qe).';
G_vect = simple(G_vect);
GG=G_vect(1:end-2);

% mass-inertial matrix
D_mtx = simple(jacobian(KE,dqe).');
D_mtx = simple(jacobian(D_mtx,dqe));
DD=D_mtx(1:end-2,1:end-2);

% Coriolis and centrifugal matrix
syms C_mtx real
n=max(size(qe));
for k=1:n
	for j=1:n
		C_mtx(k,j)=0*g;
		for i=1:n
			C_mtx(k,j)=C_mtx(k,j)+1/2*(diff(D_mtx(k,j),qe(i)) + ...
				diff(D_mtx(k,i),qe(j)) - ...
				diff(D_mtx(i,j),qe(k)))*dqe(i);
		end
	end
end
C_mtx=simple(C_mtx);
CC=C_mtx(1:end-2,1:end-2);

5
%% ground Reaction Force
Ft=D_mtx(6,1:end-2)*d2qs+C_mtx(6,1:end-2)*dqs+G_vect(6); 
Fn=D_mtx(7,1:end-2)*d2qs+C_mtx(7,1:end-2)*dqs+G_vect(7);

6
%% Impact Model

Gamma_Pe=[0;0];
Gamma_P2=p_tib2;

De=D_mtx;
E2=jacobian(p_tib2,qe)

% A X = B dq
A=[ De -E2'; 
    E2 zeros(2)];
   
B=[De*[eye(5);jacobian(Gamma_Pe,qs)] ; zeros(2,5)];

R=[ 0 1 0 0 0;
    1 0 0 0 0;
    0 0 0 1 0;
    0 0 1 0 0;
    0 0 0 0 1];
    
7

save model_5_link