%% Ode in M, D and C Matrix for my robot for 4 DOF :
%                                       theta1 , theta2 , r1, r2
home 

syms mL1 mL2 mp1 mp3 mp4 d0 L g real positive
syms B1s B2s K1 K2 real positive
syms r1 r2 theta1 theta2 real
syms D1r1 D1r2 D1theta1 D1theta2 real
syms GAMMA real


MM=[ mL1+mL2+mp3+mp4,-sin(theta2)*(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)/2,-mp4*cos(theta2),sin(theta2)*(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)/2;
    -sin(theta2)*(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)/2,mL1*(r1+L/2+d0)^2+(mL1*L^2)/0.12e2+((L^2)*cos(theta2)^2-4*L*(r1+d0+L)*cos(theta2)+(L^2)*sin(theta2)^2+4*(r1+d0+L)^2)*mL2/4+(mL2*L^2)/0.12e2+mp3*(r1+d0+L)^2+mp4*(((L+r2+d0)^2)*cos(theta2)^2-2*(L+r2+d0)*(r1+d0+L)*cos(theta2)+((L+r2+d0)^2)*sin(theta2)^2+(r1+d0+L)^2),mp4*sin(theta2)*(r1+d0+L),-L*(L*cos(theta2)^2+(-2*r1-(2*d0)-(2*L))*cos(theta2)+L*sin(theta2)^2)*mL2/4-(mL2*L^2)/0.12e2-(L+r2+d0)*((L+r2+d0)*cos(theta2)^2+(-r1-d0-L)*cos(theta2)+(L+r2+d0)*sin(theta2)^2)*mp4;
    -mp4*cos(theta2),mp4*sin(theta2)*(r1+d0+L),mp4*(cos(theta2)^2+sin(theta2)^2),0;
     sin(theta2)*(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)/2,(r1+d0+L)*(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*cos(theta2)/2-(mp4*r2^2)-(2*mp4*(L+d0)*r2)-((L+d0)^2*mp4)-(mL2*L^2)/3,0,(mp4*r2^2)+(2*mp4*(L+d0)*r2)+((L+d0)^2*mp4)+(mL2*L^2)/3];


CC=[ 0,((2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*cos(theta2)+((-2*mL2-2*mL1-2*mp4-2*mp3)*r1)+((-2*d0-2*L)*mp4)+((-2*mL2-mL1-2*mp3)*L)-(2*d0*(mL2+mp3+mL1)))*D1theta1/2-(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*cos(theta2)*D1theta2/2-mp4*D1r2*sin(theta2),mp4*sin(theta2)*(-D1theta1+D1theta2),-(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*cos(theta2)*D1theta1/2+(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*cos(theta2)*D1theta2/2+mp4*D1r2*sin(theta2);
    -((2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*cos(theta2)+((-2*mL2-2*mL1-2*mp4-2*mp3)*r1)+((-2*d0-2*L)*mp4)+((-2*mL2-mL1-2*mp3)*L)-(2*d0*(mL2+mp3+mL1)))*D1theta1/2,((-2*r2*mp4+(-2*d0-2*L)*mp4-L*mL2)*cos(theta2)+((2*mL2+2*mp3+2*mL1+2*mp4)*r1)+((2*L+2*d0)*mp4)+((2*mL2+2*mp3+mL1)*L)+(2*d0*(mL2+mp3+mL1)))*D1r1/2+sin(theta2)*(r1+d0+L)*(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*D1theta2/2+mp4*D1r2*((L+r2+d0)*sin(theta2)^2+((L+r2+d0)*cos(theta2)-r1-d0-L)*cos(theta2)),((L+r2+d0)*cos(theta2)^2+(-r1-d0-L)*cos(theta2)+(L+r2+d0)*sin(theta2)^2)*(D1theta1-D1theta2)*mp4,sin(theta2)*(r1+d0+L)*(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*D1theta1/2-sin(theta2)*(r1+d0+L)*(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*D1theta2/2-mp4*D1r2*((L+r2+d0)*sin(theta2)^2+((L+r2+d0)*cos(theta2)-r1-d0-L)*cos(theta2));
     mp4*sin(theta2)*D1theta1,-(((L+r2+d0)*cos(theta2)^2+(-r1-d0-L)*cos(theta2)+(L+r2+d0)*sin(theta2)^2)*D1theta1-(cos(theta2)^2+sin(theta2)^2)*(L+r2+d0)*D1theta2-sin(theta2)*D1r1)*mp4,0,mp4*(cos(theta2)^2+sin(theta2)^2)*(D1theta1-D1theta2)*(L+r2+d0);
     (2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*cos(theta2)*D1theta1/2,-sin(theta2)*(r1+d0+L)*(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*D1theta1/2+(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*cos(theta2)*D1r1/2-mp4*D1r2*(L+r2+d0),mp4*(-D1theta1+D1theta2)*(L+r2+d0),mp4*D1r2*(L+r2+d0)];


GG=[-g*(mL1+mL2+mp3+mp4)*cos(theta1)+K1*r1;
    -g*((r2*mp4+(L+d0)*mp4+L*mL2/2)*sin(theta1-theta2)-sin(theta1)*((2*mL2+2*mp3+2*mL1+2*mp4)*r1+(2*L+2*d0)*mp4+(2*mL2+2*mp3+mL1)*L+2*d0*(mL2+mp3+mL1))/2);
     g*mp4*cos(theta1-theta2)+K2*r2;
     g*(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*sin(theta1-theta2)/2-GAMMA];

 
 % D2q= M^-1 * (-GG -CC*D1q + U - B.D1q)   --> B= friction and damping
 
 %% define output fpr fukuda controller for 4 DOF :
%                                       theta1 , theta2 , r1, r2
home 
clear
syms mL1 mL2 mp1 mp3 mp4 d0 L g real positive
syms B1s B2s K1 K2 real positive
syms r1 r2 theta1 theta2 real
syms D1r1 D1r2 D1theta1 D1theta2 real
syms GAMMA real

 
r_vir=sqrt((L+d0+r1).^2+(L+d0+r2).^2-2*(L+d0+r1).*(L+d0+r2).*cos(((theta2))));
thetass=asin(((L+d0+r1)./r_vir).*sin(theta2));
theta_vir1= theta1+thetass;
theta_vir=atan2(((r1+L+d0).*sin(theta1)+(r2+L+d0).*sin(theta2-theta1)),((r1+L+d0).*cos(theta1)-(r2+L+d0).*cos(theta2-theta1)));

% with aproximation:  r << L+d0
r_vir_ap=sqrt((L+d0).^2+(L+d0).^2-2*(L+d0).*(L+d0).*cos(((theta2))));
theta_vir1_ap= theta1 + (pi/2 - theta2/2);


h1=theta_vir1_ap;
h2=theta_vir2;
q=[r1 theta1 r2 theta2];

for i=1:length(q)
    Dqh1(1,i)=diff(h1,q(i));
    Dqh2(1,i)=diff(h2,q(i));
end

