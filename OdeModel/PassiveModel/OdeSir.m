function OdeSir()
fprintf('_________________________________________________________________\a\n')
home
clear

g=9.81;

K1=1000;
K2=1000;
d0=0.05; % free lenght of spring

mL1=0.2;
mL2=0.2;
mp3=1.5;
mp4=0.2;
L=0.25;

B1s=100;
B2s=100;
B3pin3=0;  % for damping in actuator

Ksp=0;
theta2Init=pi;

GAMMAP1=0;
GAMMAP2=0;

Tsim=5;

OdeOpt= odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,8),'Events',@(t,X)StopCondition(t,X,d0,L));


r1=[];
theta1=[];
r2=[];
theta2=[];
Time=[0];
Yout=[];

InitState=[ 0 deg2rad(-50) 0 deg2rad(180+80) , 0 0 0 0];
%or 
%InitState=[ 0 deg2rad(-50) 0 deg2rad(-100) , 0 0 0 0];

XY=[0 0 0];

for i=1:10

    [T,Y] = ode15s(@(t,X)SirDyn(t,X,GAMMAP1,GAMMAP2,g,d0,L,mL1,mL2,mp3,mp4,K1,K2,B1s,B2s,B3pin3,Ksp,theta2Init),[Time(end) Tsim+Time(end)],InitState,OdeOpt);
    
    %AnimBot(T,Y,[0 0],L,d0);

    r10=Y(end,3);  % r1_0_NewStep= r2_end_LastStep
    r20=Y(end,1);  % r2_0_NewStep= r1_end_LastStep
    
    theta10=  Y(end,2)-Y(end,4);
    theta20= 2*pi-Y(end,4);
    
    InitState=[ r10 theta10 r20 theta20 0 0 0 0];
        
    r1     = [r1    ;Y(1:end-1,1)];
    theta1 = [theta1;Y(1:end-1,2)];
    r2     = [r2    ;Y(1:end-1,3)];
    theta2 = [theta2;Y(1:end-1,4)];
    Time=[Time; T(2:end)];
    Yout=[Yout;Y(1:end-1,:)];
    
    XY=[XY ; [XY(end,1) + (r1(end)+d0+L)*sin(theta1(end)) + (r2(end)+d0+L)*sin(theta2(end) -theta1(end))  , ...
                 XY(end,2) +     -(r1(end)+d0+L)*cos(theta1(end)) + (r2(end)+d0+L)*cos(theta2(end) -theta1(end)),...
                              Time(end)]
                              ];
    if((T(end)-T(1))==Tsim)
        break;
    end

end

Time=Time(1:end-1);  % time has two  zero at first

r_cor=sqrt((L+d0+r1).^2+(L+d0+r2).^2-2*(L+d0+r1).*(L+d0+r2).*cos(((theta2))));
thetass=asin(((L+d0+r2)./r_cor).*sin(theta2));
theta_cor= theta1+thetass;
theta_vir=atan2(((r1+L+d0).*sin(theta1)+(r2+L+d0).*sin(theta2-theta1)),((r1+L+d0).*cos(theta1)-(r2+L+d0).*cos(theta2-theta1)));

figure
plot(theta_cor)
hold all
plot(theta_vir)

Torque=zeros(size(Time));

MyPlotRcor(Time,r1*100,r2*100,theta1,theta2,r_cor,theta_cor,Torque,'Trajectory of EF')
MyPlotRcorLink(Time,r1*100,theta1,'1','Trajectory of first link')
MyPlotRcorLink(Time,r2*100,theta2,'2','Trajectory of second link')


AnimBotMultimotuin(Time,Yout,XY,L,d0);
1;

end


function DX=SirDyn(t,X,GAMMAP1,GAMMAP2,g,d0,L,mL1,mL2,mp3,mp4,K1,K2,B1s,B2s,B3pin3,Ksp,theta2Init)
r1=X(1);
theta1=X(2);
r2=X(3);
theta2=X(4);

D1r1=X(5);
% D2r1
D1theta1=X(6);
% D2theta1
D1r2=X(7);
% D2r2
D1theta2=X(8);
% D2theta2


% eq :  D.D2q + C.D1q + g -B.u
% State Form: Dx= D^-1( B.u -C.D1q -g) = D^-1 * Tx
DD2 = [[mL1+mL2+mp3+mp4,-(mL2*L*sin(theta2))/2-mp4*(L+r2+d0)*sin(theta2),-mp4*cos(theta2),(mL2*L*sin(theta2))/2+mp4*(L+r2+d0)*sin(theta2)];
       [-(mL2*L*sin(theta2))/2-mp4*(L+r2+d0)*sin(theta2),mL1*(r1+L/2+d0)^2+(mL1*L^2)/12+(mL2*((L^2*(sin(theta2))^2)/2+2*(r1+d0+L-(L*cos(theta2))/2)^2))/2+(mL2*L^2)/12+mp3*(r1+d0+L)^2+(mp4*(2*(L+r2+d0)^2*(sin(theta2))^2+2*(r1+d0+L-(L+r2+d0)*cos(theta2))^2))/2,(mp4*(2*cos(theta2)*(L+r2+d0)*sin(theta2)+2*sin(theta2)*(r1+d0+L-(L+r2+d0)*cos(theta2))))/2,(mL2*(-(L^2*(sin(theta2))^2)/2+L*cos(theta2)*(r1+d0+L-(L*cos(theta2))/2)))/2-(mL2*L^2)/12+(mp4*(-2*(L+r2+d0)^2*(sin(theta2))^2+2*(L+r2+d0)*cos(theta2)*(r1+d0+L-(L+r2+d0)*cos(theta2))))/2];
       [-mp4*cos(theta2),(mp4*(2*cos(theta2)*(L+r2+d0)*sin(theta2)+2*sin(theta2)*(r1+d0+L-(L+r2+d0)*cos(theta2))))/2,(mp4*(2*(cos(theta2))^2+2*(sin(theta2))^2))/2,0];
       [mp4*sin(theta2)*r2+mp4*sin(theta2)*d0+(mL2*L*sin(theta2))/2+mp4*L*sin(theta2),-(mL2*L^2)/3-mp4*d0^2-mp4*L^2-mp4*(r2)^2+(mL2*L*cos(theta2)*r1)/2+(mL2*L*cos(theta2)*d0)/2+mp4*L*cos(theta2)*r2+mp4*cos(theta2)*r2*d0+mp4*cos(theta2)*r1*r2+mp4*cos(theta2)*r1*d0+mp4*L*cos(theta2)*r1+2*mp4*L*cos(theta2)*d0-2*mp4*r2*d0-2*mp4*L*d0-2*mp4*L*r2+mp4*L^2*cos(theta2)+mp4*cos(theta2)*d0^2+(mL2*L^2*cos(theta2))/2,0,(mL2*L^2)/3+mp4*d0^2+mp4*L^2+mp4*(r2)^2+2*mp4*r2*d0+2*mp4*L*d0+2*mp4*L*r2]];

%%% this equation are caLcuLated symboLic and simpLified one is substitudeafter them:
% Dyn1= ...
% Dyn2= ...
% Dyn3= ...
% Dyn4= ...
%
% Tx1=simpLify( Dyn1 - DD(1,:)*[D2r1 D2theta1 D2r2  D2theta2]');
% Tx2=simpLify( Dyn2 - DD(2,:)*[D2r1 D2theta1 D2r2  D2theta2]');
% Tx3=simpLify( Dyn3 - DD(3,:)*[D2r1 D2theta1 D2r2  D2theta2]');
% Tx4=simpLify( Dyn4 - DD(4,:)*[D2r1 D2theta1 D2r2  D2theta2]');
%%% These One is simpLified version of above equations:
Tx1 = K1*r1+B1s*D1r1-g*mL1*cos(theta1)-g*mL2*cos(theta1)-g*mp3*cos(theta1)-g*mp4*cos(theta1)-(D1theta1^2*L*mL1)/2-D1theta1^2*L*mL2-D1theta1^2*L*mp3-D1theta1^2*L*mp4-D1theta1^2*d0*mL1-D1theta1^2*d0*mL2-D1theta1^2*d0*mp3-D1theta1^2*d0*mp4-D1theta1^2*mL1*r1-D1theta1^2*mL2*r1-D1theta1^2*mp3*r1-D1theta1^2*mp4*r1-2*D1r2*D1theta1*mp4*sin(theta2)+2*D1r2*D1theta2*mp4*sin(theta2)+(D1theta1^2*L*mL2*cos(theta2))/2+(D1theta2^2*L*mL2*cos(theta2))/2+D1theta1^2*L*mp4*cos(theta2)+D1theta2^2*L*mp4*cos(theta2)+D1theta1^2*d0*mp4*cos(theta2)+D1theta2^2*d0*mp4*cos(theta2)+D1theta1^2*mp4*r2*cos(theta2)+D1theta2^2*mp4*r2*cos(theta2)-D1theta1*D1theta2*L*mL2*cos(theta2)-2*D1theta1*D1theta2*L*mp4*cos(theta2)-2*D1theta1*D1theta2*d0*mp4*cos(theta2)-2*D1theta1*D1theta2*mp4*r2*cos(theta2);
Tx2 = (L*g*mL1*sin(theta1))/2-D1theta2^2*L^2*mp4*sin(theta2)-D1theta2^2*d0^2*mp4*sin(theta2)-(D1theta2^2*L^2*mL2*sin(theta2))/2+L*g*mL2*sin(theta1)+L*g*mp3*sin(theta1)+L*g*mp4*sin(theta1)+d0*g*mL1*sin(theta1)+d0*g*mL2*sin(theta1)+d0*g*mp3*sin(theta1)+d0*g*mp4*sin(theta1)+g*mL1*r1*sin(theta1)+g*mL2*r1*sin(theta1)+g*mp3*r1*sin(theta1)+g*mp4*r1*sin(theta1)-(L*g*mL2*sin(theta1-theta2))/2-L*g*mp4*sin(theta1-theta2)+D1r1*D1theta1*L*mL1+2*D1r1*D1theta1*L*mL2+2*D1r1*D1theta1*L*mp3+2*D1r1*D1theta1*L*mp4+2*D1r2*D1theta1*L*mp4-2*D1r2*D1theta2*L*mp4-d0*g*mp4*sin(theta1-theta2)+2*D1r1*D1theta1*d0*mL1+2*D1r1*D1theta1*d0*mL2+2*D1r1*D1theta1*d0*mp3+2*D1r1*D1theta1*d0*mp4+2*D1r2*D1theta1*d0*mp4-2*D1r2*D1theta2*d0*mp4-g*mp4*r2*sin(theta1-theta2)+2*D1r1*D1theta1*mL1*r1+2*D1r1*D1theta1*mL2*r1+2*D1r1*D1theta1*mp3*r1+2*D1r1*D1theta1*mp4*r1+2*D1r2*D1theta1*mp4*r2-2*D1r2*D1theta2*mp4*r2+2*D1theta1*D1theta2*d0^2*mp4*sin(theta2)-(D1theta2^2*L*d0*mL2*sin(theta2))/2-2*D1theta2^2*L*d0*mp4*sin(theta2)-(D1theta2^2*L*mL2*r1*sin(theta2))/2-D1theta2^2*L*mp4*r1*sin(theta2)-D1theta2^2*L*mp4*r2*sin(theta2)-D1theta2^2*d0*mp4*r1*sin(theta2)-D1theta2^2*d0*mp4*r2*sin(theta2)-D1theta2^2*mp4*r1*r2*sin(theta2)-D1r1*D1theta1*L*mL2*cos(theta2)-2*D1r1*D1theta1*L*mp4*cos(theta2)-2*D1r2*D1theta1*L*mp4*cos(theta2)+2*D1r2*D1theta2*L*mp4*cos(theta2)-2*D1r1*D1theta1*d0*mp4*cos(theta2)-2*D1r2*D1theta1*d0*mp4*cos(theta2)+2*D1r2*D1theta2*d0*mp4*cos(theta2)-2*D1r1*D1theta1*mp4*r2*cos(theta2)-2*D1r2*D1theta1*mp4*r1*cos(theta2)+2*D1r2*D1theta2*mp4*r1*cos(theta2)+D1theta1*D1theta2*L^2*mL2*sin(theta2)+2*D1theta1*D1theta2*L^2*mp4*sin(theta2)+2*D1theta1*D1theta2*d0*mp4*r1*sin(theta2)+2*D1theta1*D1theta2*d0*mp4*r2*sin(theta2)+2*D1theta1*D1theta2*mp4*r1*r2*sin(theta2)+D1theta1*D1theta2*L*d0*mL2*sin(theta2)+4*D1theta1*D1theta2*L*d0*mp4*sin(theta2)+D1theta1*D1theta2*L*mL2*r1*sin(theta2)+2*D1theta1*D1theta2*L*mp4*r1*sin(theta2)+2*D1theta1*D1theta2*L*mp4*r2*sin(theta2)-GAMMAP1;
Tx3 = K2*r2+B2s*D1r2+g*mp4*cos(theta1-theta2)-D1theta1^2*L*mp4-D1theta2^2*L*mp4-D1theta1^2*d0*mp4-D1theta2^2*d0*mp4-D1theta1^2*mp4*r2-D1theta2^2*mp4*r2+2*D1r1*D1theta1*mp4*sin(theta2)+D1theta1^2*L*mp4*cos(theta2)+2*D1theta1*D1theta2*L*mp4+D1theta1^2*d0*mp4*cos(theta2)+2*D1theta1*D1theta2*d0*mp4+D1theta1^2*mp4*r1*cos(theta2)+2*D1theta1*D1theta2*mp4*r2;
Tx4 = B3pin3 * D1theta2 +(L*g*mL2*sin(theta1-theta2))/2-(D1theta1^2*L^2*mL2*sin(theta2))/2-D1theta1^2*L^2*mp4*sin(theta2)-D1theta1^2*d0^2*mp4*sin(theta2)-GAMMAP2+L*g*mp4*sin(theta1-theta2)-2*D1r2*D1theta1*L*mp4+2*D1r2*D1theta2*L*mp4+d0*g*mp4*sin(theta1-theta2)-2*D1r2*D1theta1*d0*mp4+2*D1r2*D1theta2*d0*mp4+g*mp4*r2*sin(theta1-theta2)-2*D1r2*D1theta1*mp4*r2+2*D1r2*D1theta2*mp4*r2-(D1theta1^2*L*d0*mL2*sin(theta2))/2-2*D1theta1^2*L*d0*mp4*sin(theta2)-(D1theta1^2*L*mL2*r1*sin(theta2))/2-D1theta1^2*L*mp4*r1*sin(theta2)-D1theta1^2*L*mp4*r2*sin(theta2)-D1theta1^2*d0*mp4*r1*sin(theta2)-D1theta1^2*d0*mp4*r2*sin(theta2)-D1theta1^2*mp4*r1*r2*sin(theta2)+D1r1*D1theta1*L*mL2*cos(theta2)+2*D1r1*D1theta1*L*mp4*cos(theta2)+2*D1r1*D1theta1*d0*mp4*cos(theta2)+2*D1r1*D1theta1*mp4*r2*cos(theta2)-Ksp*(theta2-theta2Init);

%q=[r1theta1r2theta2]';
D1q=[D1r1 D1theta1 D1r2 D1theta2]';
D2q= -DD2^-1*[Tx1, Tx2, Tx3, Tx4]';

DX = [D1q ; D2q];
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
function MyPlotRcor(T,r1,r2,theta1,theta2,r_cor,theta_cor,Torque,namestr)


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