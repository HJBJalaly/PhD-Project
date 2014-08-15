function FukudaSimulator()
clear
clc

g=9.81;

mp1=3;
mp2=1;
L1=1.0;
L2=1.0;

ww= 2.25;

Tsim=5;

OdeOpt= odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,4));%,'Events',@(t,X)StopCondition(t,X,d0,L));


theta1=[];
theta2=[];
Torque=[];
Time=[0];
Yout=[];


InitState=[ deg2rad(-45)  deg2rad(-90) , 0  0]; % d=1.41,w=2.25
 InitState=[ deg2rad(-22.5)  deg2rad(-135) , 0  0]; % d= 0.7654 , w=2
% InitState=[ deg2rad(-30)  deg2rad(-120) , 0  0]; % d= 1 , w=2

XY=[0 0 0];

% for i=1:10

    [T,Y] = ode15s(@(t,X)SirDyn(t,X,g,L1,L2,mp1,mp2,ww),[Time(end) Tsim+Time(end)],InitState,OdeOpt);
    Tor=TorCal(T,Y,g,L1,L2,mp1,mp2,ww);


%     theta10=  Y(end,2)-Y(end,4);
%     theta20= 2*pi-Y(end,4);
%     
%     InitState=[ r10 theta10 r20 theta20 0 0 0 0];
%         
     theta1 = [theta1;Y(1:end-1,1)];
     theta2 = [theta2;Y(1:end-1,2)];
     Torque=[Torque,Tor(1:end-1)];
     Time=[Time; T(2:end)];
     
     dtheta=[0; diff(theta2)]';
     Energy=sum(abs(Torque.*dtheta))
%     Yout=[Yout;Y(1:end-1,:)];
%     
%     XY=[XY ; [XY(end,1) + (r1(end)+d0+L)*sin(theta1(end)) + (r2(end)+d0+L)*sin(theta2(end) -theta1(end))  , ...
%                  XY(end,2) +     -(r1(end)+d0+L)*cos(theta1(end)) + (r2(end)+d0+L)*cos(theta2(end) -theta1(end)),...
%                               Time(end)]
%                               ];
%     if((T(end)-T(1))==Tsim)
%         break;
%     end
% 
% end

Time=Time(1:end-1);  % time has two  zero at first
% 
r_cor=sqrt((L1).^2+(L2).^2-2*(L1).*(L2).*cos(((pi-theta2))));
theta_cor= theta1+1/2*theta2;

MyPlotRcor(Time,L1,L2,theta1,theta2,Torque,ww,InitState(1)+InitState(2)/2, 'Trajectory of EF')
PlotRcor(Time,r_cor,theta_cor,'Virtual Link')
% MyPlotRcorLink1(Time,r2,theta2,r_cor,theta_cor,'Trajectory of second link')

% frequency
[xData, yData] = prepareCurveData( Time, theta_cor);
[fitresult, gof] = fit( xData, yData, 'sin1' );
fitresult

AnimBot(T,Y,[0 0],L1,L2);
1;

end


function DX=SirDyn(t,X,g,L1,L2,m1,m2,ww)

theta1=X(1);
theta2=X(2);

D1theta1=X(3);
D1theta2=X(4);

% Output
h=theta1+1/2*theta2;
Dqh=[ 1  1/2];

% Dynamic
MM=[m1*L1^2+m2*L2^2+m2*L1^2+2*m2*L1*L2*cos(theta2) , m2*L2^2+m2*L1*L2*cos(theta2) ;
    m2*L2^2+m2*L1*L2*cos(theta2)                    , m2*L2^2 ];

% this CC, in fact is CC*Dq
CC=[ -D1theta2*m2*L2*sin(theta2)*L1*(2*D1theta1+D1theta2);
      m2*L2*sin(theta2)*L1*D1theta1^2 ];
    
GG= [ g*(m2*L2*sin(theta1+theta2)+sin(theta1)*L1*(m1+m2));
      g*m2*L2*sin(theta1+theta2)];
    
% Input
NN=MM^-1;
GAMMA2=(Dqh * [NN(2,1);NN(2,2)])^-1 * ( -ww^2*h + Dqh*NN*(CC+GG));
% rr=sqrt((L1).^2+(L2).^2-2*(L1).*(L2)*cos(((pi-theta2))));
% GAMMA2=0*(Dqh * [NN(2,1);NN(2,2)])^-1 * ( -g/rr*sin(h) + Dqh*NN*(CC+GG));

%q=[r1theta1r2theta2]';
D1q= [D1theta1 D1theta2]';
D2q= MM^-1*(-CC-GG+ [ 0;GAMMA2]);
DX=[D1q;D2q];
end

function Torque=TorCal(Time,X,g,L1,L2,m1,m2,ww)

for i=1:length(Time)

    theta1=X(i,1);
    theta2=X(i,2);

    D1theta1=X(i,3);
    D1theta2=X(i,4);

    % Output
    h=theta1+1/2*theta2;
    Dqh=[ 1  1/2];

    % Dynamic
    MM=[m1*L1^2+m2*L2^2+m2*L1^2+2*m2*L1*L2*cos(theta2) , m2*L2^2+m2*L1*L2*cos(theta2) ;
        m2*L2^2+m2*L1*L2*cos(theta2)                    , m2*L2^2 ];

    CC=[ -D1theta2*m2*L2*sin(theta2)*L1*(2*D1theta1+D1theta2);
          m2*L2*sin(theta2)*L1*D1theta1^2 ];

    GG= [ g*(m2*L2*sin(theta1+theta2)+sin(theta1)*L1*(m1+m2));
          g*m2*L2*sin(theta1+theta2)];

    % Input
    NN=MM^-1;
    Torque(i)=(Dqh * [NN(2,1);NN(2,2)])^-1 * ( -ww^2*h + Dqh*NN*(CC+GG));
%      h=theta1+1/2*theta2;
%      rr=L1*sqrt(2+2*cos(theta2));
%      Torque(i)=(Dqh * [NN(2,1);NN(2,2)])^-1 * ( -g/rr*sin(h) + Dqh*NN*(CC+GG));

end

end

function [value,isterminal,direction] = StopCondition(t,X,d0,L)
r1=X(1);
theta1=X(2);
r2=X(3);
theta2=X(4);
% r_cor=sqrt((L+d0+r1)^2+(L+d0+r2)^2-2*(L+d0+r1)*(L+d0+r2).*cos(((theta2))));
% thetass=asin(((L+d0+r1)/r_cor)*sin(theta2));
% theta_cor= theta1+thetass-pi/2;
Y = ( -(r1+d0+L)*cos(theta1) + (r2+d0+L)*cos(theta2 -theta1));

isterminal=1* (theta1> 0);
direction=-1;
% value=theta_cor;
value=Y;



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MyPlotRcor(T,r1,r2,theta1,theta2,Torque,ww,Th0,namestr)


[existFlag,figNumber] = figflag(namestr);

% if not, initialize figure
if ~existFlag,
    
    % define figure element
    h0 = figure( ...
        'Tag',                          namestr, ...
        'Name',                         namestr,...
        'NumberTitle',                    'off', ...
        'BackingStore',                   'off'); 
end %if ~existflag

% ----------------------------------
% Reset Figure to Simulation Default
% ----------------------------------

% reset axes to default properties
cla reset;
subplot(4,1,1)
plot(T,theta1)
legend('\theta_1')
grid on
line([0 T(end)],[-pi -pi]/2,'LineStyle','--')
line([0 T(end)],[pi pi]/2,'LineStyle','--')

subplot(4,1,2)
plot(T,theta2)
legend('\theta_2')
grid on
line([0 T(end)],[-pi -pi],'LineStyle','--')
line([0 T(end)],[pi pi],'LineStyle','--')

subplot(4,1,3)
plot(T,theta1+theta2/2)
hold all
plot(T,Th0*cos(T*ww),'-.')
legend('\theta','target dynamic')
grid on
hold off

line([0 T(end)],[-pi -pi]/2,'LineStyle','--')
line([0 T(end)],[pi pi]/2,'LineStyle','--')

subplot(4,1,4)
plot(T,Torque)
legend('\tau')
grid on

end

function PlotRcor(T,rr,theta,namestr)


[existFlag,figNumber] = figflag(namestr);

% if not, initialize figure
if ~existFlag,
    
    % define figure element
    h0 = figure( ...
        'Tag',                          namestr, ...
        'Name',                         namestr,...
        'NumberTitle',                    'off', ...
        'BackingStore',                   'off'); 
end %if ~existflag

% ----------------------------------
% Reset Figure to Simulation Default
% ----------------------------------

% reset axes to default properties
cla reset;
subplot(2,1,1)
plot(T,rr)
legend('r_v_i_r')
grid on

subplot(2,1,2)
plot(T,theta)
legend('\theta_v_i_r')
grid on
hold all
line([0 T(end)],[-pi -pi]/2,'LineStyle','--')
line([0 T(end)],[pi pi]/2,'LineStyle','--')
grid on
hold off

end