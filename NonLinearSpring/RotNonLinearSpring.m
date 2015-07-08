function RotNonLinearSpring()
%%
ThetaS=deg2rad(0:.1:270)';
tau=.5*(ThetaS-3*pi/4);
% tau=3*ones(size(ThetaS));
% tau=2*(1-exp(-ThetaS))+1;
% tau=(sin((ThetaS-3*pi/4)*2));

% %%
% Joint=1
% rU=4;
% ThetaStep=( (max(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)))  - min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)))) / (floor(size(Q_Opt,2)/2)));
% ThetaRange=min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2))) :ThetaStep:max(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)));
% tauMain=polyval(BetaOptimal_Opt((rU+1)*(Joint-1)+1: (Joint)*(rU+1)),ThetaRange);
% 
% Coef=polyfit(ThetaRange,tauMain,1);
% Flin=polyval(Coef,ThetaRange);
% Residule=tauMain-Flin;
% tau=Residule;
%%
 
 figure
 subplot(2,1,1)
 plot(ThetaS,tau)
% hold all
% plot(ThetaRange,Flin)
% 
% subplot(2,1,2)
% plot(ThetaRange,Residule)



%%
home
close all
K=5000;
R=.02;
l0=.075;
l0Ode=.1;

OdeOpt= odeset('RelTol',1e-4,'AbsTol',1e-4);
[Theta,L_ThetaS]=ode15s(@(theta,l)OdeSolverNonLinTRot(theta,l,tau,ThetaS,K,l0),ThetaS,l0Ode,OdeOpt);

Tan_Phi_Alpha=(tau./(K.*L_ThetaS.*(L_ThetaS-l0)));
Phi_Alpha=atan(tau./(K.*L_ThetaS.*(L_ThetaS-l0)));
Cot_Phi=(L_ThetaS.*cos(Phi_Alpha)-R)./(L_ThetaS.*sin(Phi_Alpha));
Phi=acot(Cot_Phi);
r_thetaR=sin(Phi_Alpha)./sin(Phi).*L_ThetaS;
Sin_Alpha=R*(sin(Phi)./L_ThetaS);
Alpha=asin(Sin_Alpha);
ThetaR=Theta+Alpha;

% figure
% title('force')
% plot(Theta,K*(L_ThetaS-l0))
% figure
% subplot(2,1,1)
% polar(Theta,L_ThetaS)
% subplot(2,1,2)
% plot(Theta,L_ThetaS)
% figure
% subplot(2,1,1)
% plot(Theta,Tan_Phi_Alpha)
% subplot(2,1,2)
% plot(Theta,atan(Tan_Phi_Alpha))
%% Show Time
figure('units','normalized','outerposition',[0 0 1 1])

subplot(3,3,3)
xlim(rad2deg([min(ThetaS) max(ThetaS)]))
ylim(K*([min(L_ThetaS) max(L_ThetaS)]-l0))
hold on
title('F_s')
plot(rad2deg(ThetaS),K*(L_ThetaS-l0),'linewidth',2)
StarFs=plot(rad2deg(ThetaS(1)),K*(L_ThetaS(1)-l0),'linestyle','none','marker','*','markersize',8);
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
ylabel('F_s (N)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
grid on
hold off

subplot(3,3,6)
xlim(rad2deg([min(ThetaS) max(ThetaS)]))
ylim(K*([min(L_ThetaS) max(L_ThetaS)]-l0))
hold on
title('N')
plot(rad2deg(ThetaS),K*(L_ThetaS-l0)./cos(Phi_Alpha),'linewidth',2)
StarN=plot(rad2deg(ThetaS(1)),K*(L_ThetaS(1)-l0)./cos(Phi_Alpha(1)),'linestyle','none','marker','*','markersize',8);
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
ylabel('N (N)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
grid on
hold off

subplot(3,3,9)
xlim(rad2deg([min(ThetaS) max(ThetaS)]))
ylim(([min(tau)-.1 max(tau)+.1]))
hold on
title('\tau')
plot(rad2deg(ThetaS),tau,'linewidth',2)
StarTau=plot(rad2deg(ThetaS(1)),tau(1),'linestyle','none','marker','*','markersize',8);
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
ylabel('\tau (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
grid on
hold off

subplot(3,3,7)
xlim(rad2deg([min(ThetaS) max(ThetaS)]))
ylim(([min([r_thetaR;L_ThetaS]) max([r_thetaR;L_ThetaS])]))
hold on
plot(rad2deg(ThetaS),(r_thetaR),'linewidth',2)
plot(rad2deg(ThetaS),(L_ThetaS),'r','linewidth',2)
legend('r','L','Orientation','horizontal','location','best')
StarR=plot(rad2deg(ThetaS(1)),r_thetaR(1)*100,'linestyle','none','marker','*','markersize',8);
StarL=plot(rad2deg(ThetaS(1)),(L_ThetaS(1))*100,'r','linestyle','none','marker','*','markersize',8);
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
ylabel(' (cm)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
grid on
hold off

subplot(3,3,8)
xlim(rad2deg([min(ThetaS) max(ThetaS)]))
ylim(rad2deg(([min([Phi;  Phi-Alpha]) max([Phi;  Phi-Alpha])])))
hold on
plot(rad2deg(ThetaS),rad2deg(Phi),'linewidth',2)
plot(rad2deg(ThetaS),rad2deg(Phi-Alpha),'r','linewidth',2)
legend('\phi','\phi-\alpha','Orientation','horizontal')
StarPhi=plot(rad2deg(ThetaS(1)),rad2deg(Phi(1)),'linestyle','none','marker','*','markersize',8);
StarPhiAlpha=plot(rad2deg(ThetaS(1)),rad2deg(Phi(1)-Alpha(1)),'r','linestyle','none','marker','*','markersize',8);
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
ylabel('(deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
grid on
hold off

StepT=20;
for i=1:StepT:length(ThetaS)
    
    subplot(3,3,[1,2,4,5])
    hh1=polar(0,max(L_ThetaS)+2*R);
    set(hh1,'linewidth',2)
    hold on
    polar(ThetaR,r_thetaR)
    hold all
    hhh=polar(ThetaS(1:StepT:i),L_ThetaS(1:StepT:i));
    set(hhh,'linestyle','-.')
    viscircles(L_ThetaS(i)*[cos(ThetaS(i)),sin(ThetaS(i)) ], R,'EdgeColor','b');
    plot([0 (L_ThetaS(i)+R)*cos(ThetaS(i))],[0,(L_ThetaS(i)+R)*sin(ThetaS(i))],'color','r','linewidth',2)
    plot([L_ThetaS(i)*cos(ThetaS(i)) (r_thetaR(i))*cos(ThetaR(i))],[L_ThetaS(i)*sin(ThetaS(i)),(r_thetaR(i))*sin(ThetaR(i))],'color','g')
    plot([0 (r_thetaR(i)+2*R)*cos(ThetaR(i))],[0,(r_thetaR(i)+2*R)*sin(ThetaR(i))],'color','m')
    
%     set(gca,'position',[0.08 0.325 .6 .65])
    hold off
    
    subplot(3,3,3)
    set(StarFs, 'XData',rad2deg(ThetaS(i)), 'YData',K*(L_ThetaS(i)-l0));
    
    subplot(3,3,6)
    set(StarN, 'XData',rad2deg(ThetaS(i)), 'YData',K*(L_ThetaS(i)-l0)./cos(Phi_Alpha(i)));
    
    subplot(3,3,9)
    set(StarTau, 'XData',rad2deg(ThetaS(i)), 'YData',tau(i));
    
    subplot(3,3,7)
    set(StarR, 'XData',rad2deg(ThetaS(i)), 'YData',r_thetaR(i));
    set(StarL, 'XData',rad2deg(ThetaS(i)), 'YData',L_ThetaS(i));
    
    subplot(3,3,8)
    set(StarPhi, 'XData',rad2deg(ThetaS(i)), 'YData',rad2deg(Phi(i)));
    set(StarPhiAlpha, 'XData',rad2deg(ThetaS(i)), 'YData',rad2deg(Phi(i)-Alpha(i)));
    
    drawnow
end
    
    drawnow


%%

Tau=interp1(ThetaS,tau,Theta);
Sin2Phi=2*Tau./(K*L_ThetaS.*(L_ThetaS-r0));

subplot(5,1,1)
    plot(Theta,Sin2Phi)
    ylabel('Sin2Phi')
    grid on
subplot(5,1,2)
    plot(Theta,rad2deg( asin( Sin2Phi)/2))
    ylabel('\phi')
    grid on
subplot(5,1,3)
    plot(Theta,L_ThetaS)
    ylabel('r')
    grid on
subplot(5,1,4)
    plot(Theta,K*(L_ThetaS-r0))
    ylabel('F')
    grid on
subplot(5,1,5)
    plot(Theta,tau)
    ylabel('\tau')
    xlabel('\theta')
    grid on
1;
end