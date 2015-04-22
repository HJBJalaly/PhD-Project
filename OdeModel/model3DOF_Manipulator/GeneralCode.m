

%% 2DoF 

% envirment
g=9.81;

% parameters of robot

m=1;
L=1;


f=1;
Tres=0.0001;
time=0:Tres:1/f;

Xef=0.5*cos(2*pi/f*time)+.5;
Yef=0.5*sin(2*pi/f*time)+.5;

DXef=-2*pi/f*0.5*sin(2*pi/f*time);
DYef= 2*pi/f*0.5*cos(2*pi/f*time);

D2Xef=-(2*pi/f)^2*0.5*cos(2*pi/f*time);
D2Yef=-(2*pi/f)^2*0.5*sin(2*pi/f*time);

% close
figure('name','Path (2DoF)')
plot(Xef,Yef,'linewidth',2.5,'linestyle','-','color','g')

OdeOpt= odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,5));
InitState=deg2rad([ -29.44 112.02 ,0 0 , 0]);
Pos=[Xef;Yef];
DPos=[DXef;DYef];
D2Pos=[D2Xef;D2Yef];

Kp=10000;
Kd=200;
% ode output variables

q1=InitState(1);
q2=InitState(2);

for i=1:length(time)-1
    lPos=L*[cos(q1(i))+cos(q1(i)+q2(i));
            sin(q1(i))+sin(q1(i)+q2(i))];

    cPos=Pos(:,i+1);
    dx=cPos-lPos;

    JJ=L*[-sin(q1(i))-sin(q1(i)+q2(i)), -sin(q1(i)+q2(i));
           cos(q1(i))+cos(q1(i)+q2(i)),  cos(q1(i)+q2(i))]; 

    dq=JJ^-1*dx;
    q1(i+1)=q1(i)+dq(1);
    q2(i+1)=q2(i)+dq(2);
    
end



RPos=L*[cos(q1)+cos(q1+q2);
        sin(q1)+sin(q1+q2)];
hold all
plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','-.','color','b')
hold off
axis equal


[T,Y] = ode15s(@(t,Y)SirDyn2DoF(t,Y,g,L,m,Pos,DPos,D2Pos,time,Kp,Kd), time,InitState,OdeOpt);
q1=Y(:,1)';
q2=Y(:,2)';
EnergyABS=Y(end,end)

RPos=L*[cos(q1)+cos(q1+q2);
        sin(q1)+sin(q1+q2)];
hold all
plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','--','color','r')
xlabel('x')
ylabel('y')

legend('Desired','Static path planing','dynamic path planing')
    
% AnimBot2DOF(T,Y,[0 0],L)

%% 3DoF

clear
close
home

% envirment
g=9.81*0;

% parameters of robot
m=1;
L=1;

% EF motion
f=1;
A=1;
phi=0;
Tres=0.005;
time=0:Tres:10/f;

Start=ceil(length(time)/4);

% % Circle Motion
% Xef=A*cos(2*pi/f*time)+1.5;
% Yef=A*sin(2*pi/f*time)+1;
% 
% DXef=-(2*pi/f)*A*sin(2*pi/f*time);
% DYef= (2*pi/f)*A*cos(2*pi/f*time);
% 
% D2Xef=-(2*pi/f)^2*A*cos(2*pi/f*time);
% D2Yef=-(2*pi/f)^2*A*sin(2*pi/f*time);
% q1=deg2rad(0);
% q2=deg2rad(8.031);
% q3=deg2rad(51.317);


% Line motion: Horizontal
Xef=A*cos(2*pi*f*time+phi);
Yef=2.5*ones(size(time));
DXef=-(2*pi*f)*A*sin(2*pi*f*time+phi);
DYef= 2*zeros(size(time));
D2Xef=-(2*pi*f)^2*A*cos(2*pi*f*time+phi);
D2Yef= 2*zeros(size(time));
q1=deg2rad( 22.7023); % 41.7023 deg for y=2.5m  ,  22.2 deg for y=2.0m
q2=deg2rad(29.2663); % 18.2663 deg for y=2.5m  ,  29.468 deg for y=2.0m
q3=deg2rad(80.3377); % 44.3377 deg for y=2.5m  ,  71.431 deg for y=2.0m



figure('name','Path (3DoF)')
plot(Xef,Yef,'linewidth',2.5,'linestyle','-','color','g')

Pos=[Xef;Yef];
DPos=[DXef;DYef];
D2Pos=[D2Xef;D2Yef];

%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(time)-1
    lPos=L*[cos(q1(i))+cos(q1(i)+q2(i))+cos(q1(i)+q2(i)+q3(i));
            sin(q1(i))+sin(q1(i)+q2(i))+sin(q1(i)+q2(i)+q3(i))];

    cPos=Pos(:,i);
    dx=cPos-lPos;

    JJ=L*[-sin(q1(i))-sin(q1(i)+q2(i))-sin(q1(i)+q2(i)+q3(i)), -sin(q1(i)+q2(i))-sin(q1(i)+q2(i)+q3(i)), -sin(q1(i)+q2(i)+q3(i));
           cos(q1(i))+cos(q1(i)+q2(i))+cos(q1(i)+q2(i)+q3(i)),  cos(q1(i)+q2(i))+cos(q1(i)+q2(i)+q3(i)),  cos(q1(i)+q2(i)+q3(i))]; 

    dq=JJ'*(JJ*JJ')^-1*dx;

    q1(i+1)=q1(i)+dq(1);
    q2(i+1)=q2(i)+dq(2);
    q3(i+1)=q3(i)+dq(3);
    
end



RPos=L*[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
        sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];
hold all
plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','-.','color','b')

%%%%%%%%%%%%%%%%%%%%%
OdeOpt= odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,7));
InitState=deg2rad([ 0 8.031 51.317 ,0 0 0, 0]);
DPos=[DXef;DYef];
D2Pos=[D2Xef;D2Yef];

Kp=10000;
Kd=200;

[T,Y] = ode15s(@(t,Y)SirDyn(t,Y,g,L,m,Pos,DPos,D2Pos,time,Kp,Kd), time,InitState,OdeOpt);

EnergyABS=Y(end,end)

q1=Y(:,1)';
q2=Y(:,2)';
q3=Y(:,3)';
RPos=L*[cos(q1)+cos(q1+q2)+cos(q1+q2+q3);
        sin(q1)+sin(q1+q2)+sin(q1+q2+q3)];
plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','--','color','r')
xlabel('x')
ylabel('y')
hold off
axis equal
legend('Desired','Static path planing','dynamic path planing')

Torque=TorqueCalculator(T,Y,g,m,m,m,L,L,L,Pos,DPos,D2Pos,time,Kp,Kd);
figure('name','Torque')
plot(time(Start:end),Torque(:,Start:end))
legend('\tau_1','\tau_2','\tau_3')
xlabel('time (s)')
ylabel('\tau')
