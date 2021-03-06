function OdeSirFukudaController2DOF()
fprintf('_________________________________________________________________\a\n')
home
clear

% envirment
g=9.81;

% parameters of robot

mL1=0.2;
mL2=0.2;
mp3=1.5;
mp4=0.2;
L=0.3;

B1t=0;
B2t=0;

ww= 4.5033; % for k=1000 and 0.45

% ode parameters
Tsim=3;
OdeOpt= odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,5),'Events',@(t,X)StopCondition(t,X,L));
theta_rest= InitCondinRest2DOF(0.28,L);
InitState=[ -deg2rad(theta_rest)  2*pi-2*deg2rad(theta_rest) ,...
    0 0,...
    0];

% ode output variables
theta1=-deg2rad(theta_rest);
theta2=2*pi-2*deg2rad(theta_rest);
Time=0;
Yout=InitState;
Torque=[0];

XY=[0 0 0];
EnergyAbs=0;
EnergyPos=0;
EnergyN  =0;

for i=1:10
    
    [T,Y] = ode15s(@(t,X)SirDyn(t,X,g,L,mL1,mL2,mp3,mp4,B1t,B2t,ww),[Time(end) Tsim+Time(end)],InitState,OdeOpt);
    
    %     AnimBot(T,Y,[0 0],L,d0);
    
    
    theta10=  Y(end,1)-Y(end,2);
    theta20= 2*pi-Y(end,2);
    
    InitState=[ theta10 theta20 0 0 0];
    
    theta1 = [theta1;Y(2:end,1)];
    theta2 = [theta2;Y(2:end,2)];
    Time   = [Time; T(2:end)];
    Yout   = [Yout;Y(2:end,:)];
    
    Tor=TorCal(T,Y,g,L,mL1,mL2,mp3,mp4,B1t,B2t,ww)';
    dtheta=[diff(Y(:,2));0];
    EnergyAbs=EnergyAbs+sum(abs(Tor.*dtheta));
    EnergyPos=EnergyPos+sum(( (Tor.*dtheta>0)+0 ).*(Tor.*dtheta));
    EnergyN=EnergyN+Y(end,end);
    Torque=[Torque;Tor(2:end)];
    
    
    XY=[XY ; [XY(end,1) +  (L)*sin(theta1(end)) + (L)*sin(theta2(end) -theta1(end))  ,...
              XY(end,2) + -(L)*cos(theta1(end)) + (L)*cos(theta2(end) -theta1(end))  ,...
              Time(end)]
                        ];
                    
    if(abs((T(end)-T(1))-Tsim)<1e-5)
        break;
    end
    
end

r_vir=sqrt((L).^2+(L).^2-2*(L).*(L).*cos(((theta2))));
theta_vir=theta1+(pi-theta2)/2;

% [xData, yData] = prepareCurveData(Time, theta_vir);
% [fitresult, gof] = fit( xData, yData, 'sin1' );
% fitresult

disp(' ')
disp(['EnergyAbs = ',num2str(EnergyAbs)])
disp(['EnergyPos = ',num2str(EnergyPos)])
disp(['EnergyN = ',num2str(EnergyN)])

MyPlotRcor(Time,theta1,theta2,r_vir,theta_vir,Torque,'Trajectory of EF')
% MyPlotRcorLink(Time,r1*100,theta1,'1','Trajectory of first link')
% MyPlotRcorLink(Time,r2*100,theta2,'2','Trajectory of second link')

AnimBotMultimotuin2DOF(Time,Yout,XY,L);
1;

end

function DX=SirDyn(t,X,g,L,mL1,mL2,mp3,mp4,B1t,B2t,ww)
theta1=X(1);
theta2=X(2);

D1theta1=X(3);
D1theta2=X(4);
D1q=[D1theta1 D1theta2]';

% output
h=theta1+(pi-theta2)/2;
Dqh=[  1 -1/2];


% Dynamic
MM= [-L^2*(-mL1-(4*mL2)+3*mL2*cos(theta2)-(3*mp3)-(6*mp4)+6*mp4*cos(theta2))/3, L^2*(-(2*mL2)+3*mL2*cos(theta2)-(6*mp4)+6*mp4*cos(theta2))/6;
    L^2*(-(2*mL2)+3*mL2*cos(theta2)-(6*mp4)+6*mp4*cos(theta2))/6, L^2*(mL2+3*mp4)/3];

CC= [ L^2*sin(theta2)*D1theta2*(mL2+2*mp4)/2, -L^2*sin(theta2)*(mL2+2*mp4)*(-D1theta1+D1theta2)/2;
    -L^2*sin(theta2)*D1theta1*(mL2+2*mp4)/2, 0];

GG= [L*g*(-2*mp4*sin(theta1-theta2)-mL2*sin(theta1-theta2)+2*mL2*sin(theta1)+2*mp3*sin(theta1)+2*mp4*sin(theta1)+mL1*sin(theta1))/2;
    g*mL2*L*sin(theta1-theta2)/2+g*mp4*L*sin(theta1-theta2)];


% Input: compute torque
GAMMA1=0;
GAMMA2=0;
NN=MM^-1;
GAMMA2=(Dqh * NN(:,2))^-1 * ( -ww^2*h + Dqh*NN*(CC*D1q+GG+[B1t;B2t].*D1q))-0.03;

% r_vir=sqrt((L+d0+r1).^2+(L+d0+r2).^2-2*(L+d0+r1).*(L+d0+r2).*cos(((theta2))));
% GAMMA2=(Dqh * NN(:,4))^-1 * ( -g/r_vir*sin(h) + Dqh*NN*(CC*D1q+GG))-0.005;

% D1q=[D1r1 D1theta1 D1r2 D1theta2]';
De=abs(GAMMA2*D1theta2);
D2q= MM^-1*(-CC*D1q-GG+ [ GAMMA1;GAMMA2] -[B1t;B2t].*D1q);
DX=[D1q;D2q;De];

end

function Torque=TorCal(Time,X,g,L,mL1,mL2,mp3,mp4,B1t,B2t,ww)

for i=1:length(Time)
    
    theta1=X(i,1);
    theta2=X(i,2);
    
    D1theta1=X(i,3);
    D1theta2=X(i,4);
    D1q=[D1theta1 D1theta2]';
    
    % output
    h=theta1+(pi-theta2)/2;
    Dqh=[  1  -1/2];
    
    % Dynamic
    MM= [-L^2*(-mL1-(4*mL2)+3*mL2*cos(theta2)-(3*mp3)-(6*mp4)+6*mp4*cos(theta2))/3, L^2*(-(2*mL2)+3*mL2*cos(theta2)-(6*mp4)+6*mp4*cos(theta2))/6;
        L^2*(-(2*mL2)+3*mL2*cos(theta2)-(6*mp4)+6*mp4*cos(theta2))/6, L^2*(mL2+3*mp4)/3];
    
    CC= [ L^2*sin(theta2)*D1theta2*(mL2+2*mp4)/2, -L^2*sin(theta2)*(mL2+2*mp4)*(-D1theta1+D1theta2)/2;
        -L^2*sin(theta2)*D1theta1*(mL2+2*mp4)/2, 0];
    
    GG= [L*g*(-2*mp4*sin(theta1-theta2)-mL2*sin(theta1-theta2)+2*mL2*sin(theta1)+2*mp3*sin(theta1)+2*mp4*sin(theta1)+mL1*sin(theta1))/2;
        g*mL2*L*sin(theta1-theta2)/2+g*mp4*L*sin(theta1-theta2)];
    
    % Input: compute torque
    NN=MM^-1;
    Torque(i)=(Dqh * NN(:,2))^-1 * ( -ww^2*h + Dqh*NN*(CC*D1q+GG+[B1t;B2t].*D1q))-0.03;
    
end
end

function [value,isterminal,direction] = StopCondition(t,X,L)
theta1=X(1);
theta2=X(2);
% r_cor=sqrt((L+d0+r1)^2+(L+d0+r2)^2-2*(L+d0+r1)*(L+d0+r2).*cos(((theta2))));
% thetass=asin(((L+d0+r2)/r_cor)*sin(theta2));
% theta_cor= theta1+thetass-pi/2;
Y = ( -(L)*cos(theta1) + (L)*cos(theta2 -theta1));

isterminal=1* (theta1> 0);
direction=-1;
% value=theta_cor;
value=Y;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MyPlotRcor(T,theta1,theta2,r_cor,theta_cor,Torque,namestr)


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
hold all
plot(T,theta2)
legend('\theta_1','\theta_2')
grid on
hold off
ylabel('\theta (rad)')
xlim([T(1) T(end)])

subplot(4,1,2)
plot(T,r_cor)
ylabel('r_v_i_r (cm)')
grid on
xlim([T(1) T(end)])

subplot(4,1,3)
plot(T,theta_cor)
line([0 T(end)],[-pi -pi]/2,'LineStyle','--','Color',[1 0 0]);
line([0 T(end)],[ pi  pi]/2,'LineStyle','--','Color',[1 0 0]);
ylabel('\theta_v_i_r (rad)')
grid on
hold off
xlim([T(1) T(end)])

subplot(4,1,4)
plot(T,Torque)
ylabel('Torque (N.m)')
xlabel('Time (s)')
grid on
xlim([T(1) T(end)])

end

function MyPlotRcorLink(T,rr,theta,LinkNUM,namestr)

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
ylabel(['r_',LinkNUM,' (cm)'])
grid on
xlim([T(1) T(end)])

subplot(2,1,2)
plot(T,theta)
ylabel(['\theta_L_',LinkNUM,' (rad)'])
xlabel('Time (s)')
grid on
hold all
line([0 T(end)],[-pi -pi]/2,'LineStyle','--','Color',[0.32 0.32 0.32]);
line([0 T(end)],[ pi  pi]/2,'LineStyle','--','Color',[0.32 0.32 0.32]);
if(max(abs(theta))>pi/2)
    line([0 T(end)],[-pi -pi],'LineStyle','-.','Color',[1 0 0]);
    line([0 T(end)],[ pi  pi],'LineStyle','-.','Color',[1 0 0]);
end
if(max(abs(theta))>pi)
    line([0 T(end)],[-2*pi -2*pi],'LineStyle','-.','Color',[0 0.5 0]);
    line([0 T(end)],[ 2*pi  2*pi],'LineStyle','-.','Color',[0 0.5 0]);
end
hold off
xlim([T(1) T(end)])


end