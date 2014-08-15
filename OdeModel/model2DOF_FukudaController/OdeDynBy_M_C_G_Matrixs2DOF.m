%% Ode in M, D and C Matrix for my robot for 2 DOF :
%                                       theta1 , theta2
home 
clear
syms mL1 mL2  mp3 mp4 L g real positive
syms theta1 theta2 real
syms D1theta1 D1theta2 real
syms GAMMA real


MM= [-L^2*(-mL1-(4*mL2)+3*mL2*cos(theta2)-(3*mp3)-(6*mp4)+6*mp4*cos(theta2))/3, L^2*(-(2*mL2)+3*mL2*cos(theta2)-(6*mp4)+6*mp4*cos(theta2))/6;
      L^2*(-(2*mL2)+3*mL2*cos(theta2)-(6*mp4)+6*mp4*cos(theta2))/6, L^2*(mL2+3*mp4)/3];

CC= [ L^2*sin(theta2)*D1theta2*(mL2+2*mp4)/2, -L^2*sin(theta2)*(mL2+2*mp4)*(-D1theta1+D1theta2)/2;
     -L^2*sin(theta2)*D1theta1*(mL2+2*mp4)/2, 0];

GG= [L*g*(-2*mp4*sin(theta1-theta2)-mL2*sin(theta1-theta2)+2*mL2*sin(theta1)+2*mp3*sin(theta1)+2*mp4*sin(theta1)+mL1*sin(theta1))/2;
    g*mL2*L*sin(theta1-theta2)/2+g*mp4*L*sin(theta1-theta2)-GAMMA];

whos



 
 %% define output fpr fukuda controLLer for 2 DOF :
%                                       theta1 , theta2
home 
clear
syms mL1 mL2 mp1 mp3 mp4 d0 L g real positive
syms theta1 theta2 real
syms D1theta1 D1theta2 real
syms GAMMA real

 

h=;
q=[theta1 theta2];

for i=1:length(q)
    Dqh(1,i)=diff(h,q(i));
end

