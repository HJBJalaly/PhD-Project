%
% file name: symb_model_abs_red is short for symbolic model in
% absolute coordinates, reduced form [one foot is pinned to the
% ground]
%
% This file is known to work with MATLAB Version 6.1.
% 
disp('This file generates the model used in the paper: ') 
disp(' ')
disp(' Title: STABLE WALKING OF A 7-DOF BIPED ROBOT')
disp(' Authors: Plestan, Grizzle, Westervelt and Abba')
disp(' Journal: IEEE Tran. on Robotics and Automation')
disp(' ')
disp('For a copy of a preprint, see');
disp(' ');
disp('   http://www.eecs.umich.edu/~grizzle/papers/robotics.html')
disp(' ')
disp('See the paper for the nomenclature  ')
disp(' ')
disp('Computing.....')
disp(' ')

clear
% This is for a robot with an upright trunk, two legs and knees. The
% model is for five degrees of freedom of the robot, with Foot1
% touching the ground
%
syms q1 q31 q32 q41 q42 dq1 dq31 dq32 dq41 dq42 real
syms g L1 L3 L4 M1 M3 M4 real
syms MY1 MZ1 MZ3 MZ4 XX1 XX3 XX4 real 
%
% reference = absolute angles for everything; trigonometric sense for
%             all angles, measured from the vertical,
%
q=[q31 q32 q41 q42 q1].';
dq=[dq31 dq32 dq41 dq42 dq1].';
%
% 
pH=-([-L3*sin(q31);L3*cos(q31)] + [-L4*sin(q41);L4*cos(q41)]);
pG1 = pH + L3*[-sin(q31);cos(q31)];
pG2 = pH + L3*[-sin(q32);cos(q32)];
%
%
vH =  jacobian(pH,q)*dq;
vG1 = jacobian(pG1,q)*dq;
vG2 = jacobian(pG2,q)*dq;
%
%
R_T=[cos(q1) -sin(q1); sin(q1) cos(q1)];
vH_RT = R_T.'*vH;
R_Fem1=[cos(q31) -sin(q31); sin(q31) cos(q31)];
vH_Fem1 = R_Fem1.'*vH;
R_Fem2=[cos(q32) -sin(q32); sin(q32) cos(q32)];
vH_Fem2 = R_Fem2.'*vH;
R_Tib1=[cos(q41) -sin(q41); sin(q41) cos(q41)];
vG1_Tib1 = R_Tib1.'*vG1;
R_Tib2=[cos(q42) -sin(q42); sin(q42) cos(q42)];
vG2_Tib2 = R_Tib2.'*vG2;
%
%
KET    = simplify(1/2*M1*vH.'*vH + vH_RT.'*[-MZ1;MY1]*dq1 + 1/2*XX1*dq1^2);
KEFem1 = simplify(1/2*M3*vH.'*vH + vH_Fem1.'*[-MZ3;0]*(dq31) ...
		  +1/2*XX3*(dq31)^2);
KEFem2 = simplify(1/2*M3*vH.'*vH + vH_Fem2.'*[-MZ3;0]*(dq32) ...
		  +1/2*XX3*(dq32)^2);
KETib1 = simplify(1/2*M4*vG1.'*vG1 + vG1_Tib1.'*[-MZ4;0]*(dq41) ...
		  +1/2*XX4*(dq41)^2);
KETib2 = simplify(1/2*M4*vG2.'*vG2 + vG2_Tib2.'*[-MZ4;0]*(dq42) ...
		  +1/2*XX4*(dq42)^2);
%
%
KE = simplify(KET + KEFem1 + KEFem2 + KETib1 + KETib2);
KE=simple(KE);
%
%
% centers of gravity and positions of the various members 
%
pT    = pH + (1/M1)*[-(sin(q1)*MZ1-cos(q1)*MY1);(cos(q1)*MZ1+sin(q1)*MY1)];
pFem1 = pH + (MZ3/M3)*[-sin(q31);cos(q31)];
pFem2 = pH + (MZ3/M3)*[-sin(q32);cos(q32)];
pG1   = pH +[-L3*sin(q31);L3*cos(q31)];
pG2   = pH +[-L3*sin(q32);L3*cos(q32)];
pTib1 = pG1 + [-(MZ4/M4)*sin(q41);(MZ4/M4)*cos(q41)];
pTib2 = pG2 + [-(MZ4/M4)*sin(q42);(MZ4/M4)*cos(q42)];
pFoot1= pG1 + [-L4*sin(q41);L4*cos(q41)];
pFoot2= pG2 + [-L4*sin(q42);L4*cos(q42)];
%
%
vFoot1=jacobian(pFoot1,q)*dq;
vFoot2=jacobian(pFoot2,q)*dq;
vT=jacobian(pT,q)*dq;
%
%
PE = g*(M1*pT(2) + M3*pFem1(2) + M3*pFem2(2) + M4*pTib1(2) + M4*pTib2(2));
PE = simple(PE);
%
%
% Model NOTATION: Spong and Vidyasagar, page 142, Eq. (6.3.12)
%                 D(q)ddq + C(q,dq)*dq + G(q) = B*tau + E2*F_external
%
L = KE-PE;
save work_symb_model_abs_red
G = jacobian(PE,q).';
G = simple(G);
D = simple(jacobian(KE,dq).');
D = simple(jacobian(D,dq));
save work_symb_model_abs_red
syms C real
n = max(size(q));
for k = 1:n
  for j = 1:n
    C(k,j) = 0*g;
    for i = 1:n
      C(k,j) = C(k,j)+1/2*(diff(D(k,j),q(i))+diff(D(k,i), ...
						  q(j))- ...
			   diff(D(i,j),q(k)))*dq(i);
    end
  end
end
C = simple(C);
%
% Compute the matrix for the input torques and then the contact
% points.
%
Phi_0 = [q31-q1;q32-q1;q41-q31;q42-q32;q1];
B = jacobian(Phi_0,q);
B = B.'*[eye(4,4);zeros(1,4)];
%
%
% F_ext = [F_T^1;F_N^1;F_T^2;F_N^2];
%
Phi_1 = [pFoot2];
E = jacobian(Phi_1,q);
E = E.';
save work_symb_model_abs_red

%print out results
%
disp(' ')
disp('D(q) d^2q/dt^2 + C(q,dq/dt)dq/dt + G(q) = B u')
disp(' ')
D
C
G
B

% Parameters for the IEEE TRA paper
disp(' ')
disp(['Parameters used for the paper are given in the file. The' ...
      ' values account'])
disp(['for the rotor inertias of the actuators and the gear ratio. '...
      'For the'])
disp(['actual robot, the center of mass of the torso is' ...
      ' off-axis, so MY1=0.2'])
g = 9.81;
L1 = 0.625;
L3 = 0.4;
L4 = 0.4;
M1 = 20;
M3 = 6.8;
M4 = 3.2;
MY1 = 0;
MZ1 = 4.;
MZ3 = 1.11;
MZ4 = .41;
XX1 = 2.22;
XX3 = .25+0.83;
XX4 = .10+0.83;


