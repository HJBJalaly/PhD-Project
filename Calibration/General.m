clc

Upper=[15 15 15 15 15 5 5]';
Lower=-Upper;
DataLinkagePos=@Data5LinkagePosDeterministic;
PopulationSize=1000;

[x,fval,exitflag,output,population,score] =...
    GaCalibration(Lower,Upper,PopulationSize,DataLinkagePos);


%%

L0=300;
L1=400;
L2=400;
L3=550;
L4=550;

%%
q1= 117.59;
q2= 88.4;
q3= 310.68;
q4= 33.76;
%%
Deltay=( L1*sind(q1)+L3*sind(q1+q3) ) - ( L2*sind(q2)+L4*sind(q2+q4) )
Deltax=( L1*cosd(q1)+L3*cosd(q1+q3)-L0 ) - ( L2*cosd(q2)+L4*cosd(q2+q4) )


%%

E1='L2*cos(q2)+x4-L1*cos(q1)-x3+L0';
E2='L2*sin(q2)+y4-L1*sin(q1)-y3';
E3='x3^2+y3^2-L3^2';
E4='x4^2+y4^2-L4^2';

% E1+E2+E3-->E5
E5='L2^2+L4^2+L1^2+L0^2-L3^2+2*L2*x4*cos(q2)+2*L2*y4*sin(q2)-2*L1*L2*cos(q1-q2)-2*L1*x4*cos(q1)-2*L1*y4*sin(q1)+2*L0*L2*cos(q2)+2*L0*x4-2*L0*L1*cos(q1)';

Sol=solve(E5,E4,'x4','y4')

%%
y4=(L0^2*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + L1^3*sind(q1) - L2^3*sind(q2) + L0^2*L1*sind(q1) + L1*L2^2*sind(q1) - L0^2*L2*sind(q2) - L1*L3^2*sind(q1) - L1^2*L2*sind(q2) + L1*L4^2*sind(q1) + L2*L3^2*sind(q2) - L2*L4^2*sind(q2) + L1^2*cosd(q1)^2*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + L2^2*cosd(q2)^2*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + L1^2*sind(q1)^2*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + L2^2*sind(q2)^2*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) - L0*L1^2*sind(2*q1) - L0*L2^2*sind(2*q2) - 2*L1^2*L2*cosd(q1 - q2)*sind(q1) + 2*L1*L2^2*cosd(q1 - q2)*sind(q2) - 2*L0*L1*cosd(q1)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + 2*L0*L2*cosd(q2)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) - 2*L1*L2*cosd(q1)*cosd(q2)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) - 2*L1*L2*sind(q1)*sind(q2)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + 2*L0*L1*L2*cosd(q1)*sind(q2) + 2*L0*L1*L2*cosd(q2)*sind(q1))/(2*(L0^2 - 2*cosd(q1)*L0*L1 + 2*cosd(q2)*L0*L2 + L1^2 - 2*cosd(q1 - q2)*L1*L2 + L2^2));
x4=-(L1^4*cosd(q1)^2 + L2^4*cosd(q2)^2 + L0^4 + L0^2*L1^2 + L0^2*L2^2 - L0^2*L3^2 + L0^2*L4^2 - 2*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 2*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - L1^3*sind(q1)^3*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + L2^3*sind(q2)^3*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) - 2*L0*L1^3*cosd(q1)^3 + 2*L0*L2^3*cosd(q2)^3 + 5*L0^2*L1^2*cosd(q1)^2 + 5*L0^2*L2^2*cosd(q2)^2 + L1^2*L2^2*cosd(q1)^2 + L1^2*L2^2*cosd(q2)^2 - L1^2*L3^2*cosd(q1)^2 + L1^2*L4^2*cosd(q1)^2 - L2^2*L3^2*cosd(q2)^2 + L2^2*L4^2*cosd(q2)^2 - L1^3*cosd(q1)^2*sind(q1)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + L2^3*cosd(q2)^2*sind(q2)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) - L0^2*L1*sind(q1)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + L0^2*L2*sind(q2)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) - 2*L0*L1*L2^2*cosd(q1) + 2*L0*L1*L3^2*cosd(q1) + 2*L0*L1^2*L2*cosd(q2) - 2*L0*L1*L4^2*cosd(q1) - 2*L0*L2*L3^2*cosd(q2) + 2*L0*L2*L4^2*cosd(q2) - 2*L1^3*L2*cosd(q1 - q2)*cosd(q1)^2 - 2*L1*L2^3*cosd(q1 - q2)*cosd(q2)^2 - 2*L0^2*L1*L2*cosd(q1 - q2) - 2*L1*L2^3*cosd(q1)*cosd(q2) - 2*L1^3*L2*cosd(q1)*cosd(q2) + 2*L0*L1^2*cosd(q1)*sind(q1)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + 2*L0*L2^2*cosd(q2)*sind(q2)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) - 10*L0^2*L1*L2*cosd(q1)*cosd(q2) + 2*L1*L2*L3^2*cosd(q1)*cosd(q2) - 2*L1*L2*L4^2*cosd(q1)*cosd(q2) + 4*L1^2*L2^2*cosd(q1 - q2)*cosd(q1)*cosd(q2) - L1*L2^2*cosd(q2)^2*sind(q1)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + L1^2*L2*cosd(q1)^2*sind(q2)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + 4*L0*L1^2*L2*cosd(q1 - q2)*cosd(q1) - 4*L0*L1*L2^2*cosd(q1 - q2)*cosd(q2) - 3*L1*L2^2*sind(q1)*sind(q2)^2*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) + 3*L1^2*L2*sind(q1)^2*sind(q2)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) - 6*L0*L1*L2^2*cosd(q1)*cosd(q2)^2 + 6*L0*L1^2*L2*cosd(q1)^2*cosd(q2) + 2*L1^2*L2*cosd(q1)*cosd(q2)*sind(q1)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) - 2*L1*L2^2*cosd(q1)*cosd(q2)*sind(q2)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) - 2*L0*L1*L2*cosd(q1)*sind(q2)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2) - 2*L0*L1*L2*cosd(q2)*sind(q1)*(-((L0 - L1*cosd(q1) + L2*cosd(q2))^2*(L0^4 + L1^4 + L2^4 + L3^4 + L4^4 + 4*L0^2*L1^2 + 4*L0^2*L2^2 - 2*L0^2*L3^2 + 4*L1^2*L2^2 - 2*L0^2*L4^2 - 2*L1^2*L3^2 - 2*L1^2*L4^2 - 2*L2^2*L3^2 - 2*L2^2*L4^2 - 2*L3^2*L4^2 + 2*L1^2*L2^2*cosd(2*q1 - 2*q2) - 4*L0*L1^3*cosd(q1) - 4*L0^3*L1*cosd(q1) + 4*L0*L2^3*cosd(q2) + 4*L0^3*L2*cosd(q2) - 4*L1*L2^3*cosd(q1 - q2) - 4*L1^3*L2*cosd(q1 - q2) + 2*L0^2*L1^2*cosd(2*q1) + 2*L0^2*L2^2*cosd(2*q2) + 4*L0*L1^2*L2*cosd(2*q1 - q2) - 4*L0^2*L1*L2*cosd(q1 + q2) - 8*L0*L1*L2^2*cosd(q1) + 4*L0*L1*L3^2*cosd(q1) + 8*L0*L1^2*L2*cosd(q2) + 4*L0*L1*L4^2*cosd(q1) - 4*L0*L2*L3^2*cosd(q2) - 4*L0*L2*L4^2*cosd(q2) - 8*L0^2*L1*L2*cosd(q1 - q2) - 4*L0*L1*L2^2*cosd(q1 - 2*q2) + 4*L1*L2*L3^2*cosd(q1 - q2) + 4*L1*L2*L4^2*cosd(q1 - q2)))/(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2))^2)^(1/2))/((2*L0 - 2*L1*cosd(q1) + 2*L2*cosd(q2))*(L0^2 + L1^2 + L2^2 - 2*L0*L1*cosd(q1) + 2*L0*L2*cosd(q2) - 2*L1*L2*cosd(q1 - q2)));

x3=L2*cosd(q2)+x4-L1*cosd(q1)+L0;
y3=L2*sind(q2)+y4-L1*sind(q1);
%%

tic
CostError(0,0,0,0,0,0,0)
toc
