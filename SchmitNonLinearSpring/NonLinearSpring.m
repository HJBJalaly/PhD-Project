function NonLinearSpring()
clear
home
% close all


%% Parameters and Torque Profile

k=1*137;
R=130e-3;

ThetaStep=.05;
ThetaS=deg2rad(0:ThetaStep:270);


k=150;
% R=175e-3;
% q0=300e-3;  % for sinuid
% tau=2*sin(ThetaS/3*4)+3;


% q0=130e-3;  % for cubic
% tau=0.21*(ThetaS-0.75*pi).*(ThetaS-0.25*pi).*(ThetaS-1.25*pi)+2.5;

k=2*170;
q0=190e-3;  % for constant
tau=4*ones(size(ThetaS));
  
% q0=99e-3;  % for exp
% tau=2*(1-exp(-ThetaS))+1;        

%%  Calculate (r,theta) of spool

sum=0;
for i=1:length(ThetaS)
    J(i)=tau(i)/sqrt(2*k*sum+(k*q0)^2);
    sum=sum+tau(i)*deg2rad(ThetaStep);
end

DJ=differential(J,ThetaS,deg2rad(ThetaStep));

r= sqrt( J.^2 + (DJ.^2 .* (R^2 - J.^2))./((DJ+sqrt(R^2-J.^2)).^2) );

if(~(isreal(r)))
    error('complex r is generated.')
    disp('Choose Differet parameters.')
    return;
end

ThetaR=-ThetaS+acos((J.^2+sign(DJ).*sqrt((r.^2-J.^2).*(R^2-J.^2))  )./(R*r));


%% calculate L , Lspool, Fspring and generted torque
LL=sqrt(r.^2+R^2-2*r*R.*cos(ThetaR+ThetaS));
Thetaa=acos((r.^2+LL.^2-R^2)./(2*r.*LL));
L0=LL(1);

F=k*q0;
Lspool=0;
for i=2:length(ThetaR)
    Lspool(i)=Lspool(i-1)-r(i-1)*diff(ThetaR(i-1:i));
    F(i)=k*(Lspool(i)+LL(i)-L0+q0);
end
Tauc=F.*r.*sin(Thetaa);


%% showtime

figure
subplot(1,3,1)
plot(rad2deg(ThetaS),tau)
xlabel('\theta (deg)')
ylabel('\tau')
grid on
subplot(1,3,2)
polar([ThetaR ThetaR(1)],[r r(1)]*100)
subplot(1,3,3)
plot(rad2deg(ThetaS),r*100)
xlabel('\theta (deg)')
ylabel('radius (cm)')
grid on

Xp=700;
Yp=150;
set(gcf,'Units','points', 'Position', [100, 10,100+ Xp,100+ Yp])


figure
DeltaF=max(F)-min(F);
DeltaTau=max(Tauc)-min(Tauc);

set(gcf,'Units','points', 'Position', [100, 300,100+ Xp,100+ 2*Yp])

for i=1:10:length(ThetaS)
    subplot(2,2,[1,3])
    % fake point for ajdust size of polar
    polar(0, R*1.1*100)
    hold on
    % spool
    h2=polar([ThetaR ThetaR(1)]+(ThetaS(i)),[r r(1)]*100);
    set(h2,'linewidth',3,'color','b')
    hold all
    % raduis line
    h3=polar([ThetaR(i)+ThetaS(i) 0],[ r(i) 0]*100,'-.');
    set(h3,'linewidth',2,'color','g')
    % tangent line
    h1=polar([ThetaR(1:i)+ThetaS(i) 0],[r(1:i) R]*100);
    set(h1,'linewidth',2,'color','r')
    % path of tagency point
    h4=polar(ThetaS(1:i)+ThetaR(1:i),r(1:i)*100,'--') ;
    set(h4,'linewidth',1,'color','m')
    hold off
    
    subplot(2,2,2)
    plot(rad2deg(ThetaS(1:i)),F(1:i),'-.')
    hold on
    plot(rad2deg(ThetaS(i)),F(i),'Marker','*','LineStyle','none')
    xlabel('\theta_s (deg)')
    ylabel('F_s_p_r_i_n_g (N)')
    hold off
    axis([0 rad2deg(ThetaS(end)) min(F)-0.1*DeltaF max(F)+0.1*DeltaF])
    grid on

    
    subplot(2,2,4)
    plot(rad2deg(ThetaS(1:i)),Tauc(1:i),'LineStyle','-.','linewidth',2)
    hold on
    plot(rad2deg(ThetaS(i)),Tauc(i),'Marker','*','LineStyle','none')
    plot(rad2deg(ThetaS),tau,'g')
    xlabel('\theta_s (deg)')
    ylabel('\tau (N.m)')
    hold off
    axis([0 rad2deg(ThetaS(end)) min(Tauc)-0.1*DeltaTau max(Tauc)+0.1*DeltaTau])
    grid on
    drawnow
end
