
%% when m3 is zero

clear
home
syms g d0 L mL m2  m4 K1 K2 B1 B2 B3pin real positvie
syms r1 r2  theta1 theta2 real
syms D1r1 D2r1 D1r2 D2r2 D1theta1 D2theta1 D1theta2 D2theta2 real
syms GAMMA  real


DD2 = [[2*mL+m4,-(mL*L*sin(theta2))/2-m4*(L+r2+d0)*sin(theta2),-m4*cos(theta2),(mL*L*sin(theta2))/2+m4*(L+r2+d0)*sin(theta2)];
       [-(mL*L*sin(theta2))/2-m4*(L+r2+d0)*sin(theta2),mL*(r1+L/2+d0)^2+(mL*L^2)/6+(mL*((L^2*(sin(theta2))^2)/2+2*(r1+d0+L-(L*cos(theta2))/2)^2))/2+(m4*(2*(L+r2+d0)^2*(sin(theta2))^2+2*(r1+d0+L-(L+r2+d0)*cos(theta2))^2))/2,(m4*(2*cos(theta2)*(L+r2+d0)*sin(theta2)+2*sin(theta2)*(r1+d0+L-(L+r2+d0)*cos(theta2))))/2,(mL*(-(L^2*(sin(theta2))^2)/2+L*cos(theta2)*(r1+d0+L-(L*cos(theta2))/2)))/2-(mL*L^2)/12+(m4*(-2*(L+r2+d0)^2*(sin(theta2))^2+2*(L+r2+d0)*cos(theta2)*(r1+d0+L-(L+r2+d0)*cos(theta2))))/2];
       [-m4*cos(theta2),(m4*(2*cos(theta2)*(L+r2+d0)*sin(theta2)+2*sin(theta2)*(r1+d0+L-(L+r2+d0)*cos(theta2))))/2,(m4*(2*(cos(theta2))^2+2*(sin(theta2))^2))/2,0];
       [(mL*L*sin(theta2))/2+m4*L*sin(theta2)+m4*sin(theta2)*r2+m4*sin(theta2)*d0,-m4*(r2)^2-m4*L^2-m4*d0^2+m4*L*cos(theta2)*r2+m4*cos(theta2)*r2*d0+m4*cos(theta2)*r1*r2+m4*cos(theta2)*r1*d0+(mL*L*cos(theta2)*r1)/2+(mL*L*cos(theta2)*d0)/2+m4*L*cos(theta2)*r1+2*m4*L*cos(theta2)*d0-(mL*L^2)/3-2*m4*r2*d0-2*m4*L*d0-2*m4*L*r2+(mL*L^2*cos(theta2))/2+m4*L^2*cos(theta2)+m4*cos(theta2)*d0^2,0,m4*(r2)^2+m4*L^2+m4*d0^2+(mL*L^2)/3+2*m4*r2*d0+2*m4*L*d0+2*m4*L*r2]];


Dyn1t= mL * D2r1 + mL * ((2 * D2r1) - L * (D2theta1 - D2theta2) * sin(theta2) - L * (D1theta1 - D1theta2) * cos(theta2) * D1theta2) / 0.2e1 + m4 * ((2 * D2r1) - 0.2e1 * D1r2 * (D1theta1 - D1theta2) * sin(theta2) - (2 * L + 2 * r2 + 2 * d0) * (D2theta1 - D2theta2) * sin(theta2) - (2 * L + 2 * r2 + 2 * d0) * (D1theta1 - D1theta2) * cos(theta2) * D1theta2 - 0.2e1 * D2r2 * cos(theta2) + 0.2e1 * D1r2 * sin(theta2) * D1theta2) / 0.2e1 - mL * (r1 + L / 0.2e1 + d0) * D1theta1 ^ 2 - mL * ((r1 + d0 + L) * D1theta1 - L * cos(theta2) * (D1theta1 - D1theta2) / 0.2e1) * D1theta1 - m4 * ((r1 + d0 + L) * D1theta1 - (L + r2 + d0) * cos(theta2) * (D1theta1 - D1theta2) + D1r2 * sin(theta2)) * D1theta1 - g * (0.2e1 * mL * cos(theta1) + m4 * cos(theta1)) + K1 * r1 + B1 * D1r1;
Dyn2t= 0.2e1 * mL * (r1 + L / 0.2e1 + d0) * D1theta1 * D1r1 + mL * (r1 + L / 0.2e1 + d0) ^ 2 * D2theta1 + D2theta1 * mL * (L ^ 2) / 0.12e2 + mL * (-(D2r1 - L * (D2theta1 - D2theta2) * sin(theta2) / 0.2e1 - L * (D1theta1 - D1theta2) * cos(theta2) * D1theta2 / 0.2e1) * L * sin(theta2) - (D1r1 - L * (D1theta1 - D1theta2) * sin(theta2) / 0.2e1) * L * cos(theta2) * D1theta2 + (0.2e1 * D1r1 * D1theta1 + 0.2e1 * (r1 + d0 + L) * D2theta1 + L * sin(theta2) * D1theta2 * (D1theta1 - D1theta2) - L * cos(theta2) * (D2theta1 - D2theta2)) * (r1 + d0 + L - L * cos(theta2) / 0.2e1) + (0.2e1 * (r1 + d0 + L) * D1theta1 - L * cos(theta2) * (D1theta1 - D1theta2)) * (D1r1 + L * sin(theta2) * D1theta2 / 0.2e1)) / 0.2e1 + (D2theta1 / 0.12e2 - D2theta2 / 0.12e2) * mL * (L ^ 2) + m4 * (-(0.2e1 * D2r1 - 0.2e1 * D1r2 * (D1theta1 - D1theta2) * sin(theta2) - 0.2e1 * (L + r2 + d0) * (D2theta1 - D2theta2) * sin(theta2) - 0.2e1 * (L + r2 + d0) * (D1theta1 - D1theta2) * cos(theta2) * D1theta2 - 0.2e1 * D2r2 * cos(theta2) + 0.2e1 * D1r2 * sin(theta2) * D1theta2) * (L + r2 + d0) * sin(theta2) - (0.2e1 * D1r1 - 0.2e1 * (L + r2 + d0) * (D1theta1 - D1theta2) * sin(theta2) - 0.2e1 * D1r2 * cos(theta2)) * D1r2 * sin(theta2) - (0.2e1 * D1r1 - 0.2e1 * (L + r2 + d0) * (D1theta1 - D1theta2) * sin(theta2) - 0.2e1 * D1r2 * cos(theta2)) * (L + r2 + d0) * cos(theta2) * D1theta2 + (0.2e1 * D1r1 * D1theta1 + 0.2e1 * (r1 + d0 + L) * D2theta1 - 0.2e1 * D1r2 * cos(theta2) * (D1theta1 - D1theta2) + 0.2e1 * (L + r2 + d0) * sin(theta2) * D1theta2 * (D1theta1 - D1theta2) - 0.2e1 * (L + r2 + d0) * cos(theta2) * (D2theta1 - D2theta2) + 0.2e1 * D2r2 * sin(theta2) + 0.2e1 * D1r2 * cos(theta2) * D1theta2) * (r1 + d0 + L - (L + r2 + d0) * cos(theta2)) + (0.2e1 * (r1 + d0 + L) * D1theta1 - 0.2e1 * (L + r2 + d0) * cos(theta2) * (D1theta1 - D1theta2) + 0.2e1 * D1r2 * sin(theta2)) * (D1r1 - D1r2 * cos(theta2) + (L + r2 + d0) * sin(theta2) * D1theta2)) / 0.2e1 - g * (-mL * (r1 + L / 0.2e1 + d0) * sin(theta1) + mL * (-(r1 + d0 + L) * sin(theta1) + L * sin(theta1 - theta2) / 0.2e1) + m4 * (-(r1 + d0 + L) * sin(theta1) + (L + r2 + d0) * sin(theta1 - theta2)));
Dyn3t= m4 * (-((2 * D2r1) - 0.2e1 * D1r2 * (D1theta1 - D1theta2) * sin(theta2) - 0.2e1 * (L + r2 + d0) * (D2theta1 - D2theta2) * sin(theta2) - 0.2e1 * (L + r2 + d0) * (D1theta1 - D1theta2) * cos(theta2) * D1theta2 - 0.2e1 * D2r2 * cos(theta2) + 0.2e1 * D1r2 * sin(theta2) * D1theta2) * cos(theta2) + ((2 * D1r1) - 0.2e1 * (L + r2 + d0) * (D1theta1 - D1theta2) * sin(theta2) - 0.2e1 * D1r2 * cos(theta2)) * sin(theta2) * D1theta2 + (0.2e1 * D1r1 * D1theta1 + 0.2e1 * (r1 + d0 + L) * D2theta1 - 0.2e1 * D1r2 * cos(theta2) * (D1theta1 - D1theta2) + 0.2e1 * (L + r2 + d0) * sin(theta2) * D1theta2 * (D1theta1 - D1theta2) - 0.2e1 * (L + r2 + d0) * cos(theta2) * (D2theta1 - D2theta2) + 0.2e1 * D2r2 * sin(theta2) + 0.2e1 * D1r2 * cos(theta2) * D1theta2) * sin(theta2) + (0.2e1 * (r1 + d0 + L) * D1theta1 - 0.2e1 * (L + r2 + d0) * cos(theta2) * (D1theta1 - D1theta2) + 0.2e1 * D1r2 * sin(theta2)) * cos(theta2) * D1theta2) / 0.2e1 - m4 * (-((2 * D1r1) - 0.2e1 * (L + r2 + d0) * (D1theta1 - D1theta2) * sin(theta2) - 0.2e1 * D1r2 * cos(theta2)) * (D1theta1 - D1theta2) * sin(theta2) - (0.2e1 * (r1 + d0 + L) * D1theta1 - 0.2e1 * (L + r2 + d0) * cos(theta2) * (D1theta1 - D1theta2) + 0.2e1 * D1r2 * sin(theta2)) * cos(theta2) * (D1theta1 - D1theta2)) / 0.2e1 + g * m4 * cos(theta1 - theta2) + K2 * r2 + B2 * D1r2;
Dyn4t= (2 * m4 * L * D2theta2 * r2) + (2 * m4 * L * D2theta2 * d0) - (2 * m4 * r2 * D2theta1 * d0) + (2 * m4 * r2 * D2theta2 * d0) + (2 * m4 * D1theta2 * D1r2 * L) + (2 * m4 * D1theta2 * D1r2 * r2) + (2 * m4 * D1theta2 * D1r2 * d0) - (2 * m4 * L * D2theta1 * d0) - (2 * m4 * D1r2 * D1theta1 * L) - (2 * m4 * D1r2 * D1theta1 * r2) - (2 * m4 * D1r2 * D1theta1 * d0) - (2 * m4 * L * D2theta1 * r2) - GAMMA - (D2theta1 * mL * L ^ 2) / 0.3e1 + (m4 * L ^ 2 * D2theta2) + (m4 * d0 ^ 2 * D2theta2) - (m4 * L ^ 2 * D2theta1) - (m4 * D2theta1 * d0 ^ 2) + m4 * L * cos(theta2) * D2theta1 * r2 + m4 * cos(theta2) * r2 * D2theta1 * d0 + 0.2e1 * m4 * cos(theta2) * D1r1 * D1theta1 * r2 + 0.2e1 * m4 * cos(theta2) * D1r1 * D1theta1 * d0 + m4 * cos(theta2) * D2theta1 * r1 * r2 + m4 * cos(theta2) * D2theta1 * r1 * d0 - mL * L * sin(theta2) * (D1theta1 ^ 2) * r1 / 0.2e1 - mL * L * sin(theta2) * (D1theta1 ^ 2) * d0 / 0.2e1 + mL * L * cos(theta2) * D1r1 * D1theta1 + mL * L * cos(theta2) * D2theta1 * r1 / 0.2e1 + mL * L * cos(theta2) * D2theta1 * d0 / 0.2e1 + 0.2e1 * m4 * L * cos(theta2) * D1r1 * D1theta1 + m4 * L * cos(theta2) * D2theta1 * r1 + 0.2e1 * m4 * L * cos(theta2) * D2theta1 * d0 - m4 * L * sin(theta2) * (D1theta1 ^ 2) * r1 - 0.2e1 * m4 * L * sin(theta2) * (D1theta1 ^ 2) * d0 - m4 * (D1theta1 ^ 2) * r1 * sin(theta2) * r2 - m4 * (D1theta1 ^ 2) * r1 * sin(theta2) * d0 - m4 * (D1theta1 ^ 2) * d0 * sin(theta2) * r2 - m4 * (D1theta1 ^ 2) * L * sin(theta2) * r2 - (m4 * r2 ^ 2 * D2theta1) + (m4 * r2 ^ 2 * D2theta2) + (mL * L ^ 2 * D2theta2) / 0.3e1 + mL * L * sin(theta2) * D2r1 / 0.2e1 + mL * (L ^ 2) * cos(theta2) * D2theta1 / 0.2e1 + m4 * L * sin(theta2) * D2r1 + m4 * (L ^ 2) * cos(theta2) * D2theta1 + m4 * sin(theta2) * D2r1 * r2 + m4 * sin(theta2) * D2r1 * d0 + m4 * cos(theta2) * D2theta1 * (d0 ^ 2) - mL * (L ^ 2) * sin(theta2) * (D1theta1 ^ 2) / 0.2e1 - m4 * (L ^ 2) * sin(theta2) * (D1theta1 ^ 2) - m4 * (D1theta1 ^ 2) * (d0 ^ 2) * sin(theta2) + g * mL * L * sin(theta1 - theta2) / 0.2e1 + g * m4 * sin(theta1 - theta2) * L + g * m4 * sin(theta1 - theta2) * r2 + g * m4 * sin(theta1 - theta2) * d0 +  B3pin * D1theta2;


% eq :  D.D2q + C.D1q + g -B.u=0   >
% State Form: Dx= D^-1( B.u -C.D1q -g) = D^-1 * Tx
Tx1=simplify( Dyn1t - DD2(1,:)*[D2r1 D2theta1 D2r2  D2theta2]');
Tx2=simplify( Dyn2t - DD2(2,:)*[D2r1 D2theta1 D2r2  D2theta2]');
Tx3=simplify( Dyn3t - DD2(3,:)*[D2r1 D2theta1 D2r2  D2theta2]');
Tx4=simplify( Dyn4t - DD2(4,:)*[D2r1 D2theta1 D2r2  D2theta2]');

q   = [r1  theta1  r2  theta2]';
D1q = [D1r1 D1theta1 D1r2 D1theta2]';
D2q = -DD2^-1*[Tx1, Tx2, Tx3, Tx4]';


DX = [D1q ; D2q];

%% when m3 is non-zero

clear
home
syms g d0 L mL1 mL2 mp3 mp4 K1 K2 B1 B2 real positvie
syms r1 r2  theta1 theta2 real
syms D1r1 D2r1 D1r2 D2r2 D1theta1 D2theta1 D1theta2 D2theta2 real
syms GAMMA  real

DD2 = [[mL1+mL2+mp3+mp4,-(mL2*L*sin(theta2))/2-mp4*(L+r2+d0)*sin(theta2),-mp4*cos(theta2),(mL2*L*sin(theta2))/2+mp4*(L+r2+d0)*sin(theta2)];
       [-(mL2*L*sin(theta2))/2-mp4*(L+r2+d0)*sin(theta2),mL1*(r1+L/2+d0)^2+(mL1*L^2)/12+(mL2*((L^2*(sin(theta2))^2)/2+2*(r1+d0+L-(L*cos(theta2))/2)^2))/2+(mL2*L^2)/12+mp3*(r1+d0+L)^2+(mp4*(2*(L+r2+d0)^2*(sin(theta2))^2+2*(r1+d0+L-(L+r2+d0)*cos(theta2))^2))/2,(mp4*(2*cos(theta2)*(L+r2+d0)*sin(theta2)+2*sin(theta2)*(r1+d0+L-(L+r2+d0)*cos(theta2))))/2,(mL2*(-(L^2*(sin(theta2))^2)/2+L*cos(theta2)*(r1+d0+L-(L*cos(theta2))/2)))/2-(mL2*L^2)/12+(mp4*(-2*(L+r2+d0)^2*(sin(theta2))^2+2*(L+r2+d0)*cos(theta2)*(r1+d0+L-(L+r2+d0)*cos(theta2))))/2];
       [-mp4*cos(theta2),(mp4*(2*cos(theta2)*(L+r2+d0)*sin(theta2)+2*sin(theta2)*(r1+d0+L-(L+r2+d0)*cos(theta2))))/2,(mp4*(2*(cos(theta2))^2+2*(sin(theta2))^2))/2,0];
       [mp4*sin(theta2)*r2+mp4*sin(theta2)*d0+(mL2*L*sin(theta2))/2+mp4*L*sin(theta2),-(mL2*L^2)/3-mp4*d0^2-mp4*L^2-mp4*(r2)^2+(mL2*L*cos(theta2)*r1)/2+(mL2*L*cos(theta2)*d0)/2+mp4*L*cos(theta2)*r2+mp4*cos(theta2)*r2*d0+mp4*cos(theta2)*r1*r2+mp4*cos(theta2)*r1*d0+mp4*L*cos(theta2)*r1+2*mp4*L*cos(theta2)*d0-2*mp4*r2*d0-2*mp4*L*d0-2*mp4*L*r2+mp4*L^2*cos(theta2)+mp4*cos(theta2)*d0^2+(mL2*L^2*cos(theta2))/2,0,(mL2*L^2)/3+mp4*d0^2+mp4*L^2+mp4*(r2)^2+2*mp4*r2*d0+2*mp4*L*d0+2*mp4*L*r2]];


Dyn1t = mL1*D2r1+mL2*((2*D2r1)-L*(D2theta1-D2theta2)*sin(theta2)-L*(D1theta1-D1theta2)*cos(theta2)*D1theta2)/0.2e1+mp3*D2r1+mp4*((2*D2r1)-0.2e1*D1r2*(D1theta1-D1theta2)*sin(theta2)-(0.2e1*L+(2*r2)+(2*d0))*(D2theta1-D2theta2)*sin(theta2)-(0.2e1*L+(2*r2)+(2*d0))*(D1theta1-D1theta2)*cos(theta2)*D1theta2-0.2e1*D2r2*cos(theta2)+0.2e1*D1r2*sin(theta2)*D1theta2)/0.2e1-mL1*(r1+L/0.2e1+d0)*D1theta1^2-mL2*((r1+d0+L)*D1theta1-L*cos(theta2)*(D1theta1-D1theta2)/0.2e1)*D1theta1-mp3*(r1+d0+L)*D1theta1^2-mp4*((r1+d0+L)*D1theta1-(L+r2+d0)*cos(theta2)*(D1theta1-D1theta2)+D1r2*sin(theta2))*D1theta1-g*(mL1*cos(theta1)+mL2*cos(theta1)+mp3*cos(theta1)+mp4*cos(theta1))+K1*r1+B1*D1r1;
Dyn2t = 0.2e1*mL1*(r1+L/0.2e1+d0)*D1theta1*D1r1+mL1*(r1+L/0.2e1+d0)^2*D2theta1+D2theta1*mL1*L^2/0.12e2+mL2*(-(D2r1-L*(D2theta1-D2theta2)*sin(theta2)/0.2e1-L*(D1theta1-D1theta2)*cos(theta2)*D1theta2/0.2e1)*L*sin(theta2)-(D1r1-L*(D1theta1-D1theta2)*sin(theta2)/0.2e1)*L*cos(theta2)*D1theta2+(0.2e1*D1r1*D1theta1+0.2e1*(r1+d0+L)*D2theta1+L*sin(theta2)*D1theta2*(D1theta1-D1theta2)-L*cos(theta2)*(D2theta1-D2theta2))*(r1+d0+L-L*cos(theta2)/0.2e1)+(0.2e1*(r1+d0+L)*D1theta1-L*cos(theta2)*(D1theta1-D1theta2))*(D1r1+L*sin(theta2)*D1theta2/0.2e1))/0.2e1+(D2theta1/0.12e2-D2theta2/0.12e2)*mL2*L^2+0.2e1*mp3*(r1+d0+L)*D1theta1*D1r1+mp3*(r1+d0+L)^2*D2theta1+mp4*(-(0.2e1*D2r1-0.2e1*D1r2*(D1theta1-D1theta2)*sin(theta2)-0.2e1*(L+r2+d0)*(D2theta1-D2theta2)*sin(theta2)-0.2e1*(L+r2+d0)*(D1theta1-D1theta2)*cos(theta2)*D1theta2-0.2e1*D2r2*cos(theta2)+0.2e1*D1r2*sin(theta2)*D1theta2)*(L+r2+d0)*sin(theta2)-(0.2e1*D1r1-0.2e1*(L+r2+d0)*(D1theta1-D1theta2)*sin(theta2)-0.2e1*D1r2*cos(theta2))*D1r2*sin(theta2)-(0.2e1*D1r1-0.2e1*(L+r2+d0)*(D1theta1-D1theta2)*sin(theta2)-0.2e1*D1r2*cos(theta2))*(L+r2+d0)*cos(theta2)*D1theta2+(0.2e1*D1r1*D1theta1+0.2e1*(r1+d0+L)*D2theta1-0.2e1*D1r2*cos(theta2)*(D1theta1-D1theta2)+0.2e1*(L+r2+d0)*sin(theta2)*D1theta2*(D1theta1-D1theta2)-0.2e1*(L+r2+d0)*cos(theta2)*(D2theta1-D2theta2)+0.2e1*D2r2*sin(theta2)+0.2e1*D1r2*cos(theta2)*D1theta2)*(r1+d0+L-(L+r2+d0)*cos(theta2))+(0.2e1*(r1+d0+L)*D1theta1-0.2e1*(L+r2+d0)*cos(theta2)*(D1theta1-D1theta2)+0.2e1*D1r2*sin(theta2))*(D1r1-D1r2*cos(theta2)+(L+r2+d0)*sin(theta2)*D1theta2))/0.2e1-g*(-mL1*(r1+L/0.2e1+d0)*sin(theta1)+mL2*(-(r1+d0+L)*sin(theta1)+L*sin(theta1-theta2)/0.2e1)-mp3*(r1+d0+L)*sin(theta1)+mp4*(-(r1+d0+L)*sin(theta1)+(L+r2+d0)*sin(theta1-theta2)));
Dyn3t = mp4*(-((2*D2r1)-0.2e1*D1r2*(D1theta1-D1theta2)*sin(theta2)-0.2e1*(L+r2+d0)*(D2theta1-D2theta2)*sin(theta2)-0.2e1*(L+r2+d0)*(D1theta1-D1theta2)*cos(theta2)*D1theta2-0.2e1*D2r2*cos(theta2)+0.2e1*D1r2*sin(theta2)*D1theta2)*cos(theta2)+((2*D1r1)-0.2e1*(L+r2+d0)*(D1theta1-D1theta2)*sin(theta2)-0.2e1*D1r2*cos(theta2))*sin(theta2)*D1theta2+(0.2e1*D1r1*D1theta1+0.2e1*(r1+d0+L)*D2theta1-0.2e1*D1r2*cos(theta2)*(D1theta1-D1theta2)+0.2e1*(L+r2+d0)*sin(theta2)*D1theta2*(D1theta1-D1theta2)-0.2e1*(L+r2+d0)*cos(theta2)*(D2theta1-D2theta2)+0.2e1*D2r2*sin(theta2)+0.2e1*D1r2*cos(theta2)*D1theta2)*sin(theta2)+(0.2e1*(r1+d0+L)*D1theta1-0.2e1*(L+r2+d0)*cos(theta2)*(D1theta1-D1theta2)+0.2e1*D1r2*sin(theta2))*cos(theta2)*D1theta2)/0.2e1-mp4*(-((2*D1r1)-0.2e1*(L+r2+d0)*(D1theta1-D1theta2)*sin(theta2)-0.2e1*D1r2*cos(theta2))*(D1theta1-D1theta2)*sin(theta2)-(0.2e1*(r1+d0+L)*D1theta1-0.2e1*(L+r2+d0)*cos(theta2)*(D1theta1-D1theta2)+0.2e1*D1r2*sin(theta2))*cos(theta2)*(D1theta1-D1theta2))/0.2e1+g*mp4*cos(theta1-theta2)+K2*r2+B2*D1r2;
Dyn4t =-mL2*L^2*D2theta1/0.3e1+mL2*L^2*D2theta2/0.3e1+mp4*L*sin(theta2)*D2r1+mp4*L^2*cos(theta2)*D2theta1+mp4*sin(theta2)*D2r1*r2+mp4*sin(theta2)*D2r1*d0+mp4*cos(theta2)*D2theta1*d0^2-mL2*L^2*sin(theta2)*D1theta1^2/0.2e1-mp4*L^2*sin(theta2)*D1theta1^2-mp4*D1theta1^2*d0^2*sin(theta2)+g*mL2*L*sin(theta1-theta2)/0.2e1+g*mp4*sin(theta1-theta2)*L+g*mp4*sin(theta1-theta2)*r2+g*mp4*sin(theta1-theta2)*d0+mL2*L*sin(theta2)*D2r1/0.2e1+mL2*L^2*cos(theta2)*D2theta1/0.2e1+mL2*L*cos(theta2)*D1r1*D1theta1+mL2*L*cos(theta2)*D2theta1*r1/0.2e1+mL2*L*cos(theta2)*D2theta1*d0/0.2e1+mp4*L*cos(theta2)*D2theta1*r2+mp4*cos(theta2)*r2*D2theta1*d0+0.2e1*mp4*cos(theta2)*D1r1*D1theta1*r2+0.2e1*mp4*cos(theta2)*D1r1*D1theta1*d0+mp4*cos(theta2)*D2theta1*r1*r2+mp4*cos(theta2)*D2theta1*r1*d0+0.2e1*mp4*L*cos(theta2)*D1r1*D1theta1+mp4*L*cos(theta2)*D2theta1*r1+0.2e1*mp4*L*cos(theta2)*D2theta1*d0-mL2*L*sin(theta2)*D1theta1^2*r1/0.2e1-mL2*L*sin(theta2)*D1theta1^2*d0/0.2e1-mp4*L*sin(theta2)*D1theta1^2*r1-0.2e1*mp4*L*sin(theta2)*D1theta1^2*d0-mp4*D1theta1^2*r1*sin(theta2)*r2-mp4*D1theta1^2*r1*sin(theta2)*d0-mp4*D1theta1^2*d0*sin(theta2)*r2-mp4*D1theta1^2*L*sin(theta2)*r2-mp4*L^2*D2theta1+mp4*L^2*D2theta2+mp4*r2^2*D2theta2-mp4*D2theta1*d0^2-0.2e1*mp4*r2*D2theta1*d0+0.2e1*mp4*L*D2theta2*d0-0.2e1*mp4*L*D2theta1*d0+0.2e1*mp4*L*D2theta2*r2+0.2e1*mp4*r2*D2theta2*d0-0.2e1*mp4*D1r2*D1theta1*L-0.2e1*mp4*D1r2*D1theta1*r2-0.2e1*mp4*D1r2*D1theta1*d0+0.2e1*mp4*D1r2*D1theta2*L+0.2e1*mp4*D1r2*D1theta2*r2+0.2e1*mp4*D1r2*D1theta2*d0-0.2e1*mp4*L*D2theta1*r2-GAMMA+mp4*d0^2*D2theta2-mp4*r2^2*D2theta1;



%eq:D.D2q+C.D1q+g-B.u=0>
%StateForm:Dx=D^-1(B.u-C.D1q-g)=D^-1*Tx
Tx1=simplify( Dyn1t - DD2(1,:)*[D2r1 D2theta1 D2r2  D2theta2]');
Tx2=simplify( Dyn2t - DD2(2,:)*[D2r1 D2theta1 D2r2  D2theta2]');
Tx3=simplify( Dyn3t - DD2(3,:)*[D2r1 D2theta1 D2r2  D2theta2]');
Tx4=simplify( Dyn4t - DD2(4,:)*[D2r1 D2theta1 D2r2  D2theta2]');

q   = [r1  theta1  r2  theta2]';
D1q = [D1r1 D1theta1 D1r2 D1theta2]';
D2q = -DD2^-1*[Tx1, Tx2, Tx3, Tx4]';


DX = [D1q ; D2q];

%% Impact De With non-zero m3

clear
home
syms g d0 L mL1 mL2 mp3 mp4 K1 K2 B1 B2 real positvie
syms r1 r2  theta1 theta2 real
syms D1r1 D2r1 D1r2 D2r2 D1theta1 D2theta1 D1theta2 D2theta2 real
syms GAMMA  real

De= [mp4+mL1+mL2+mp3,-sin(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*cos(-theta2+theta1)+cos(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*sin(-theta2+theta1),-mp4*cos(theta1)*cos(-theta2+theta1)-mp4*sin(theta1)*sin(-theta2+theta1),sin(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*cos(-theta2+theta1)-cos(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*sin(-theta2+theta1),sin(theta1)*(mp4+mL1+mL2+mp3),-cos(theta1)*(mp4+mL1+mL2+mp3);
     -sin(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*cos(-theta2+theta1)+cos(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*sin(-theta2+theta1),-2*(r1+d0+L)*sin(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*sin(-theta2+theta1)-2*(r1+d0+L)*cos(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*cos(-theta2+theta1)+1/6*(6*mp3+6*mp4+6*mL1+6*mL2)*(r1)^2+1/6*((12*d0+12*L)*mp4+(12*mp3+6*mL1+12*mL2)*L+12*d0*(mp3+mL1+mL2))*r1+mp4*(r2)^2+2*(L+d0)*mp4*r2+2*(L+d0)^2*mp4+1/6*(2*mL1+8*mL2+6*mp3)*L^2+d0*(2*mp3+mL1+2*mL2)*L+d0^2*(mp3+mL1+mL2),-mp4*cos(theta1)*(r1+d0+L)*sin(-theta2+theta1)+mp4*sin(theta1)*(r1+d0+L)*cos(-theta2+theta1),(r1+d0+L)*sin(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*sin(-theta2+theta1)+(r1+d0+L)*cos(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*cos(-theta2+theta1)-mp4*(r2)^2-2*(L+d0)*mp4*r2-(L+d0)^2*mp4-1/3*mL2*L^2,1/6*(-6*mp4*r2+(-6*L-6*d0)*mp4-3*mL2*L)*cos(-theta2+theta1)+cos(theta1)*((mp4+mL1+mL2+mp3)*r1+(L+d0)*mp4+(mL2+1/2*mL1+mp3)*L+d0*(mp3+mL1+mL2)),1/6*(-6*mp4*r2+(-6*L-6*d0)*mp4-3*mL2*L)*sin(-theta2+theta1)+sin(theta1)*((mp4+mL1+mL2+mp3)*r1+(L+d0)*mp4+(mL2+1/2*mL1+mp3)*L+d0*(mp3+mL1+mL2));
     -mp4*cos(theta1)*cos(-theta2+theta1)-mp4*sin(theta1)*sin(-theta2+theta1),-mp4*cos(theta1)*(r1+d0+L)*sin(-theta2+theta1)+mp4*sin(theta1)*(r1+d0+L)*cos(-theta2+theta1),mp4,0,-mp4*sin(-theta2+theta1),mp4*cos(-theta2+theta1);
     sin(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*cos(-theta2+theta1)-cos(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*sin(-theta2+theta1),(r1+d0+L)*sin(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*sin(-theta2+theta1)+(r1+d0+L)*cos(theta1)*(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*cos(-theta2+theta1)-mp4*(r2)^2-2*(L+d0)*mp4*r2-(L+d0)^2*mp4-1/3*mL2*L^2,0,mp4*(r2)^2+2*(L+d0)*mp4*r2+(L+d0)^2*mp4+1/3*mL2*L^2,(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*cos(-theta2+theta1),(mp4*r2+(L+d0)*mp4+1/2*mL2*L)*sin(-theta2+theta1);
     sin(theta1)*(mp4+mL1+mL2+mp3),1/2*(-2*mp4*r2+(-2*d0-2*L)*mp4-mL2*L)*cos(-theta2+theta1)+1/2*((2*mp4+2*mp3+2*mL2+2*mL1)*r1+(2*d0+2*L)*mp4+(2*mp3+mL1+2*mL2)*L+2*d0*(mp3+mL1+mL2))*cos(theta1),-mp4*sin(-theta2+theta1),1/2*(2*mp4*r2+(2*d0+2*L)*mp4+mL2*L)*cos(-theta2+theta1),mp4+mL1+mL2+mp3,0;
     -cos(theta1)*(mp4+mL1+mL2+mp3),1/2*(-2*mp4*r2+(-2*d0-2*L)*mp4-mL2*L)*sin(-theta2+theta1)+sin(theta1)*((mp4+mL1+mL2+mp3)*r1+(L+d0)*mp4+(mL2+1/2*mL1+mp3)*L+d0*(mp3+mL1+mL2)),mp4*cos(-theta2+theta1),1/2*(2*mp4*r2+(2*d0+2*L)*mp4+mL2*L)*sin(-theta2+theta1),0,mp4+mL1+mL2+mp3];


De=simplify(De);

 E2=[ sin(theta1),(r1+d0+L)*cos(theta1)-(r2+d0+L)*cos(theta2-theta1),-sin(-theta2+theta1), (r2 + d0 + L)*cos(theta2 - theta1),1, 0;
     -cos(theta1), (r1 + d0 + L) * sin(theta1) + (r2 + d0 + L) * sin(theta2 - theta1), -cos(theta2 - theta1), -(r2 + d0 + L) * sin(theta2 - theta1), 0, 1;];

R = [ 0  0  1   0;
      0  1  0  -1;
      1  0  0   0;
      0  0  0 - 1];
 
DF=-simplify( simplify((simplify(E2*De^-1*E2'))^-1) * (E2*[1 0;0 1; 0 0; 0 0]));
Dqd=simplify(simplify(simplify(De^-1)*E2')*DF) + [1 0; 0 1; 0 0; 0 0];
DD=simplify([R ,[0 0;0 0]]*Dqd);


