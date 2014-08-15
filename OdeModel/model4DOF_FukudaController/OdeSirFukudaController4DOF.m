function OdeSirFukudaController4DOF()
fprintf('_________________________________________________________________\a\n')
home
clear

% envirment
g=9.81;

% parameters of robot
K1=1000;
K2=1000;
d0=0.05; % free lenght of spring

mL1=0.2;
mL2=0.2;
mp3=1.5;
mp4=0.2;
L=0.25;

B1s=1; % damping
B2s=1;
B1t=0.0;
B2t=0.0;

ww= 4.3425; % for k=1000
% ww= 4.5630; % for k=100000

% ode parameters
Tsim=3;
OdeOpt= odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,9),'Events',@(t,X)StopCondition(t,X,d0,L));

[rr_rest,theta_rest] = InitCondinRest4DOF(0.45,L,d0,mL1,mp3,K1,g)
InitState=[rr_rest  ,-deg2rad(theta_rest)  , rr_rest , 2*pi-2*deg2rad(theta_rest) ,...
            0 0 0 0,...
            0];


% ode output variables
r1=InitState(1);
theta1=InitState(2);
r2=InitState(3);
theta2=InitState(4);
Time=[0];
Yout=InitState;
Torque=[0];

XY=[0 0 0];
EnergyAbs=0;
EnergyPos=0;
EnergyN  =0;
    


for i=1:10

    [T,Y] = ode15s(@(t,X)SirDyn(t,X,g,d0,L,mL1,mL2,mp3,mp4,K1,K2,B1s,B2s,B1t,B2t,ww),[Time(end) Tsim+Time(end)],InitState,OdeOpt);
    
%     AnimBot(T,Y,[0 0],L,d0);

    r10=Y(end,3);  % r1_0_NewStep= r2_end_LastStep
    r20=Y(end,1);  % r2_0_NewStep= r1_end_LastStep
    
    theta10=  Y(end,2)-Y(end,4);
    theta20= 2*pi-Y(end,4);
    
    InitState=[ r10 theta10 r20 theta20 0 0 0 0 0];
        
    r1     = [r1    ;Y(2:end,1)];
    theta1 = [theta1;Y(2:end,2)];
    r2     = [r2    ;Y(2:end,3)];
    theta2 = [theta2;Y(2:end,4)];
    Time   = [Time; T(2:end)];
    Yout   = [Yout;Y(2:end,:)];
    
    Tor=TorCal(T,Y,g,d0,L,mL1,mL2,mp3,mp4,K1,K2,B1s,B2s,B1t,B2t,ww)';
    dtheta=[0;diff(Y(:,4))];
    EnergyAbs=EnergyAbs+sum(abs(Tor.*dtheta));
    EnergyPos=EnergyPos+sum(( (Tor.*dtheta>0)+0 ).*(Tor.*dtheta));
    EnergyN=EnergyN+Y(end,end);
    
    Torque=[Torque; Tor(2:end)];
          
    
    XY=[XY ; [XY(end,1) + (r1(end)+d0+L)*sin(theta1(end)) + (r2(end)+d0+L)*sin(theta2(end) -theta1(end))  , ...
                 XY(end,2) +     -(r1(end)+d0+L)*cos(theta1(end)) + (r2(end)+d0+L)*cos(theta2(end) -theta1(end)),...
                              Time(end)]
                              ];
                          
    if(abs((T(end)-T(1))-Tsim)<1e-5)
        break;
    end

end

r_vir=sqrt((L+d0+r1).^2+(L+d0+r2).^2-2*(L+d0+r1).*(L+d0+r2).*cos(((theta2))));
thetass=asin(((L+d0+r2)./r_vir).*sin(theta2));
theta_vir= theta1+thetass;
theta_vir_es=theta1+(pi-theta2)/2;


[xData, yData] = prepareCurveData(Time, theta_vir);
[fitresult, gof] = fit( xData, yData, 'sin1' );
fitresult

disp(' ')
disp(['EnergyAbs = ',num2str(EnergyAbs)])
disp(['EnergyPos = ',num2str(EnergyPos)])
disp(['EnergyN = ',num2str(EnergyN)])

MyPlotRcor(Time,r1*100,r2*100,theta1,theta2,r_vir,theta_vir,theta_vir_es,Torque,'Trajectory of EF')
MyPlotRcorLink(Time,r1*100,theta1,'1','Trajectory of first link')
MyPlotRcorLink(Time,r2*100,theta2,'2','Trajectory of second link')

AnimBotMultimotuin4DOF(Time,Yout,XY,L,d0);
1;

end


function DX=SirDyn(t,X,g,d0,L,mL1,mL2,mp3,mp4,K1,K2,B1s,B2s,B1t,B2t,ww)
r1=X(1);
theta1=X(2);
r2=X(3);
theta2=X(4);

D1r1=X(5);
D1theta1=X(6);
D1r2=X(7);
D1theta2=X(8);
D1q=[D1r1 D1theta1 D1r2 D1theta2]';

% output
h=theta1+(pi-theta2)/2;
Dqh=[ 0 1  0 -1/2];



% Dynamic
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
    g*(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*sin(theta1-theta2)/2];


% Input: compute torque
GAMMA1=0;
NN=MM^-1;
GAMMA2=(Dqh * NN(:,4))^-1 * ( -ww^2*h + Dqh*NN*(CC*D1q+GG + [B1s B1t B2s B2t]'.*D1q))-0.05;
% r_vir=sqrt((L+d0+r1).^2+(L+d0+r2).^2-2*(L+d0+r1).*(L+d0+r2).*cos(((theta2))));
% GAMMA2=(Dqh * NN(:,4))^-1 * ( -g/r_vir*sin(h) + Dqh*NN*(CC*D1q+GG))-0.005;

% D1q=[D1r1 D1theta1 D1r2 D1theta2]';
De=abs(GAMMA2*D1theta2);
D2q= MM^-1*(-CC*D1q-GG+ [ 0;GAMMA1;0;GAMMA2] - [B1s B1t B2s B2t]'.*D1q);
DX=[D1q;D2q;De];

end

function Torque=TorCal(Time,X,g,d0,L,mL1,mL2,mp3,mp4,K1,K2,B1s,B2s,B1t,B2t,ww)

for i=1:length(Time)
    
    r1=X(i,1);
    theta1=X(i,2);
    r2=X(i,3);
    theta2=X(i,4);
    
    D1r1=X(i,5);
    D1theta1=X(i,6);
    D1r2=X(i,7);
    D1theta2=X(i,8);
    D1q=[D1r1 D1theta1 D1r2 D1theta2]';
    
    % output
    h=theta1+(pi-theta2)/2;
    Dqh=[ 0 1  0 -1/2];
    
    
    % Dynamic
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
        g*(2*r2*mp4+(2*L+2*d0)*mp4+L*mL2)*sin(theta1-theta2)/2];
    
    
    % Input: compute torque
    NN=MM^-1;
    Torque(i)=(Dqh * NN(:,4))^-1 * ( -ww^2*h + Dqh*NN*(CC*D1q+GG - [B1s B1t B2s B2t]'.*D1q) )-0.05;
    
end
end

function [value,isterminal,direction] = StopCondition(t,X,d0,L)
r1=X(1);
theta1=X(2);
r2=X(3);
theta2=X(4);
% r_cor=sqrt((L+d0+r1)^2+(L+d0+r2)^2-2*(L+d0+r1)*(L+d0+r2).*cos(((theta2))));
% thetass=asin(((L+d0+r2)/r_cor)*sin(theta2));
% theta_cor= theta1+thetass-pi/2;
Y = ( -(r1+d0+L)*cos(theta1) + (r2+d0+L)*cos(theta2 -theta1));

isterminal=1* (theta1> 0);
direction=-1;
% value=theta_cor;
value=Y;



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MyPlotRcor(T,r1,r2,theta1,theta2,r_cor,theta_cor,theta_cor_es,Torque,namestr)


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
subplot(5,1,1)
plot(T,r1)
hold all
plot(T,r2)
legend('r_1','r_2')
grid on
hold off
ylabel('r (cm)')
xlim([T(1) T(end)])

subplot(5,1,2)
plot(T,theta1)
hold all
plot(T,theta2)
legend('\theta_1','\theta_2')
grid on
hold off
ylabel('\theta (rad)')
xlim([T(1) T(end)])

subplot(5,1,3)
plot(T,r_cor)
ylabel('r_v_i_r (cm)')
grid on
xlim([T(1) T(end)])

subplot(5,1,4)
plot(T,theta_cor)
hold all
plot(T,theta_cor_es)
legend('real','estimate')
line([0 T(end)],[-pi -pi]/2,'LineStyle','--','Color',[1 0 0]);
line([0 T(end)],[ pi  pi]/2,'LineStyle','--','Color',[1 0 0]);
ylabel('\theta_v_i_r (rad)')
grid on
hold off
xlim([T(1) T(end)])

subplot(5,1,5)
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