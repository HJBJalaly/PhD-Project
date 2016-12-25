function RotNonLinearSpringAppr(ThetaS,tau,K,R,l0,lOde,m,FigName,Xlabel)





%%
home
if(nargin==0)
    ThetaS=deg2rad(0:.05:180)';
% linear
%     tau=.125*(ThetaS)-.1178*2.5; 
% constant
     tau=.3*ones(size(ThetaS));
% exp
%     tau=(2*(1-exp((-ThetaS) )))/2;
% tanh
%      tau=.149*(tanh(-(-ThetaS+3*pi/4)*1));
% Cubic
%     tau=.02*(ThetaS-3*pi/4).^3;
% Sine
      tau=.25*(sin(2*(ThetaS-3*pi/4)));
      tau=2.5*(sin(1*(ThetaS)));
plot(tau)
    m=.0385;
    K=15000;
    R=.015/2;
    l0=.022;
    lOde=.025;
    FigName='Test';
    Xlabel='q';
% Bearing:  http://www.astbearings.com/product.html?product=696H    
end
%%
OdeOpt= odeset('RelTol',1e-6,'AbsTol',1e-6);
[Theta,L_ThetaS]=ode15s(@(theta,l)OdeSolverNonLinTRot(theta,l,tau,ThetaS,K,l0),ThetaS,lOde,OdeOpt);

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
% figure
% hp=polar([ThetaR ;ThetaR(1)],[r_thetaR ;r_thetaR(1)]);
% patch( get(hp,'XData'), get(hp,'YData'), 'g')

figure
set(gcf,'position',[300 25 500 675])
subplot(3,2,[1,2,3,4])
    hold off
    hh1=Mypolar(0,100*(max(L_ThetaS)+R));
    set(hh1,'linewidth',2);
    hold on
    hh=Mypolar([ThetaR; ThetaR(1)],[ 100*r_thetaR ;100*r_thetaR(1)]);
    set(hh,'linewidth',3);
    
    hold all
    hhh=polar(ThetaS,100*L_ThetaS);
    set(hhh,'linestyle','-.');
    viscircles(100*L_ThetaS(1)*[cos(ThetaS(1)),sin(ThetaS(1)) ],100*R,'EdgeColor',[0.31 0.31 0.3]);
    plot([0 100*(L_ThetaS(1)+R)*cos(ThetaS(1))],[0,100*(L_ThetaS(1)+R)*sin(ThetaS(1))],'color','r','linewidth',2)
    plot([100*L_ThetaS(1)*cos(ThetaS(1)) 100*(r_thetaR(1))*cos(ThetaR(1))],[100*L_ThetaS(1)*sin(ThetaS(1)),100*(r_thetaR(1))*sin(ThetaR(1))],'color','g')
    plot([0 100*(r_thetaR(1)+2*R)*cos(ThetaR(1))],[0,100*(r_thetaR(1)+2*R)*sin(ThetaR(1))],'color','m')
    

    
hss=subplot(3,2,5);
    plot(rad2deg(ThetaS),tau,'linewidth',2)
    xlabel([Xlabel,' (deg)'],'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    ylabel('\tau (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10','FontWeight','bold');
    grid on
    set(gca,'FontSize',10,'FontWeight','bold','FontName','mwa_cmb10')
    xlim([0 270])
%     ylim([min(tau) max(tau)]*1.2)
    set(gca,'XTickLabel',{'0','','','90','','','180','','','270'},...
    'XTick',[0 30 60 90 120 150 180 210 240 270])

hss=subplot(3,2,6);
    plot(rad2deg(ThetaS(2:end)),diff(tau(1:end))./diff(ThetaS(1:end)),'linewidth',2)
    ylabel('Stiffness (N.m/deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    grid on
    xlabel([Xlabel,' (deg)'],'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    set(gca,'FontSize',10,'FontWeight','bold','FontName','mwa_cmb10')
    xlim([0 270])
%     ylim([min(diff(tau(1:end))./diff(ThetaS(1:end))) max(diff(tau(1:end))./diff(ThetaS(1:end)))]*1.2)
    set(gca,'XTickLabel',{'0','','','90','','','180','','','270'},...
    'XTick',[0 30 60 90 120 150 180 210 240 270])

%% ShowTime for paper
figure
set(gcf,'position',[0 10 800 660])
subplot(3,4,12)
Fs=K*(L_ThetaS-l0);
plot(rad2deg(ThetaS),Fs,'linewidth',2)
xlim(rad2deg([min(ThetaS) max(ThetaS)]))
ylim(([min(min(Fs)*[1.1 .91]), max(max(Fs)*[1.1 .91])]))
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
ylabel('F_s (N)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
grid on
set(gca,'FontSize',10,'FontWeight','bold','FontName','mwa_cmb10')
set(gca,'XTickLabel',{'0','','','90','','','180','','','270'},...
'XTick',[0 30 60 90 120 150 180 210 240 270])

subplot(3,4,11)
NN=K*(L_ThetaS-l0)./cos(Phi_Alpha);
plot(rad2deg(ThetaS),NN,'linewidth',2)
xlim(rad2deg([min(ThetaS) max(ThetaS)]))
ylim(([min(min(NN)*[1.1 .91]), max(max(NN)*[1.1 .91])]))
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
ylabel('N (N)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
grid on
set(gca,'FontSize',10,'FontWeight','bold','FontName','mwa_cmb10')
set(gca,'XTickLabel',{'0','','','90','','','180','','','270'},...
    'XTick',[0 30 60 90 120 150 180 210 240 270])

subplot(3,4,[1,2,5,6])
plot(rad2deg(ThetaS),tau,'linewidth',2)
xlim(rad2deg([min(ThetaS) max(ThetaS)]))
ylim(([min(min(tau)*[1.1 .91]), max(max(tau)*[1.1 .91])]))
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
ylabel('\tau (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
grid on
set(gca,'FontSize',10,'FontWeight','bold','FontName','mwa_cmb10')
set(gca,'XTickLabel',{'0','','','90','','','180','','','270'},...
    'XTick',[0 30 60 90 120 150 180 210 240 270])
set(gca,'position',[0.175 0.5 0.3 0.4])

subplot(3,4,9)
plot(rad2deg(ThetaS),100*(r_thetaR),'linewidth',2)
hold on
plot(rad2deg(ThetaS),100*(L_ThetaS),'r','linewidth',2)
xlim(rad2deg([min(ThetaS) max(ThetaS)]))
ylim([min([r_thetaR;L_ThetaS])*.91  max([r_thetaR;L_ThetaS])*1.1]*100)
h=legend('$r$','$l$','Orientation','vertical','location','best');
set(h,'interpreter','latex');
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
ylabel('Length (cm)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
grid on
set(gca,'FontSize',10,'FontWeight','bold','FontName','mwa_cmb10')
set(gca,'XTickLabel',{'0','','','90','','','180','','','270'},...
    'XTick',[0 30 60 90 120 150 180 210 240 270])

subplot(3,4,10)
plot(rad2deg(ThetaS),rad2deg(Phi),'linewidth',2)
hold on
plot(rad2deg(ThetaS),rad2deg(Phi-Alpha),'r','linewidth',2)
xlim(rad2deg([min(ThetaS) max(ThetaS)]))
ylim(rad2deg([min(min([Phi;  Phi-Alpha])*[1.2 .82]) max(max([Phi;  Phi-Alpha])*[1.2 .82])]))
legend('\phi','\phi-\alpha','Orientation','vertical','location','best');
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
ylabel('Angle (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
grid on
set(gca,'FontSize',10,'FontWeight','bold','FontName','mwa_cmb10')
set(gca,'XTickLabel',{'0','','','90','','','180','','','270'},...
    'XTick',[0 30 60 90 120 150 180 210 240 270])


subplot(3,4,[3,4,7,8])
hh1=Mypolar(0,100*(max(L_ThetaS)+R));
set(hh1,'linewidth',2);
hold on
hh=Mypolar([ThetaR],[ 100*r_thetaR]);
set(hh,'linewidth',3);
hh=Mypolar([ThetaR(end); ThetaR(1)],[ 100*r_thetaR(end) ;100*r_thetaR(1)]);
set(hh,'linewidth',3,'linestyle','--');

hold all
% hhh=polar(ThetaS,100*L_ThetaS);
% set(hhh,'linestyle','-.');
viscircles(100*L_ThetaS(1)*[cos(ThetaS(1)),sin(ThetaS(1)) ],100*R,'EdgeColor',[0.31 0.31 0.3]);
plot([0 100*(L_ThetaS(1)+R)*cos(ThetaS(1))],[0,100*(L_ThetaS(1)+R)*sin(ThetaS(1))],'color','r','linewidth',2)
plot([100*L_ThetaS(1)*cos(ThetaS(1)) 100*(r_thetaR(1))*cos(ThetaR(1))],[100*L_ThetaS(1)*sin(ThetaS(1)),100*(r_thetaR(1))*sin(ThetaR(1))],'color','g')
plot([0 100*(r_thetaR(1)+2*R)*cos(ThetaR(1))],[0,100*(r_thetaR(1)+2*R)*sin(ThetaR(1))],'color','m')


%%
figure('units','normalized','outerposition',[0 0 1 1])

subplot(3,3,9)
xlim(rad2deg([min(ThetaS) max(ThetaS)]))
ylim(K*([min(L_ThetaS) max(L_ThetaS)]-l0))
hold on
% title('F_s')
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

subplot(3,3,3)
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
ylim(([min([r_thetaR;L_ThetaS]) max([r_thetaR;L_ThetaS])])*100)
hold on
plot(rad2deg(ThetaS),100*(r_thetaR),'linewidth',2)
plot(rad2deg(ThetaS),100*(L_ThetaS),'r','linewidth',2)
legend('r','L','Orientation','horizontal','location','best')
StarR=plot(rad2deg(ThetaS(1)),r_thetaR(1)*100,'linestyle','none','marker','*','markersize',8);
StarL=plot(rad2deg(ThetaS(1)),(L_ThetaS(1))*100,'r','linestyle','none','marker','*','markersize',8);
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
ylabel('Length (cm)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
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
ylabel('Angle (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
grid on
hold off


subplot(3,3,[1,2,4,5])
hh1=Mypolar(0,100*(max(L_ThetaS)+R));
set(hh1,'linewidth',2);
hold on
hh=Mypolar([ThetaR; ThetaR(1)],[ 100*r_thetaR ;100*r_thetaR(1)]);
set(hh,'linewidth',3);
hh=Mypolar([ThetaR(end); ThetaR(1)],[ 100*r_thetaR(end) ;100*r_thetaR(1)]);
set(hh,'linewidth',3,'linestyle','--');


% Create a VideoWriter object to write the video out to a new, different file.
% writerObj = VideoWriter('tanhyper.avi');
% open(writerObj);
% Need to change from the default renderer to zbuffer to get it to work right.
% openGL doesn't work and Painters is way too slow.
set(gcf, 'renderer', 'zbuffer');


StepT=floor(length(ThetaS)/50);
% loops = length(1:StepT:length(ThetaS));
% F(loops) = struct('cdata',[],'colormap',[]);
FramesFolder = './ImageExmaple';
frame=0;

for i=1:StepT:length(ThetaS)
    frame=frame+1;
    subplot(3,3,[1,2,4,5])
    hh1=polar(0,(max(L_ThetaS)+2*R)*100);
    set(hh1,'linewidth',2)
    hold on
    hh=polar(ThetaR,r_thetaR*100);
    set(hh,'linewidth',3)
    hh=Mypolar([ThetaR(end); ThetaR(1)],[ 100*r_thetaR(end) ;100*r_thetaR(1)]);
    set(hh,'linewidth',3,'linestyle','--');

    
    hold all
    hhh=polar(ThetaS(1:StepT:i),L_ThetaS(1:StepT:i)*100);
    set(hhh,'linestyle','-.')
    viscircles(100*L_ThetaS(i)*[cos(ThetaS(i)),sin(ThetaS(i)) ],100*R,'EdgeColor',[0.31 0.31 0.3]);
    plot([0 100*(L_ThetaS(i)+R)*cos(ThetaS(i))],[0,100*(L_ThetaS(i)+R)*sin(ThetaS(i))],'color','r','linewidth',2)
    plot([100*L_ThetaS(i)*cos(ThetaS(i)) 100*(r_thetaR(i))*cos(ThetaR(i))],[100*L_ThetaS(i)*sin(ThetaS(i)),100*(r_thetaR(i))*sin(ThetaR(i))],'color','g')
    plot([0 100*(r_thetaR(i)+2*R)*cos(ThetaR(i))],[0,100*(r_thetaR(i)+2*R)*sin(ThetaR(i))],'color','m')
    
%     set(gca,'position',[0.08 0.325 .6 .65])
    hold off
    
    subplot(3,3,3)
    set(StarFs, 'XData',rad2deg(ThetaS(i)), 'YData',K*(L_ThetaS(i)-l0));
    
    subplot(3,3,6)
    set(StarN, 'XData',rad2deg(ThetaS(i)), 'YData',K*(L_ThetaS(i)-l0)./cos(Phi_Alpha(i)));
    
    subplot(3,3,9)
    set(StarTau, 'XData',rad2deg(ThetaS(i)), 'YData',tau(i));
    
    subplot(3,3,7)
    set(StarR, 'XData',rad2deg(ThetaS(i)), 'YData',100*r_thetaR(i));
    set(StarL, 'XData',rad2deg(ThetaS(i)), 'YData',100*L_ThetaS(i));
    
    subplot(3,3,8)
    set(StarPhi, 'XData',rad2deg(ThetaS(i)), 'YData',rad2deg(Phi(i)));
    set(StarPhiAlpha, 'XData',rad2deg(ThetaS(i)), 'YData',rad2deg(Phi(i)-Alpha(i)));
    
    drawnow
%     print(spintf('-r%d',90*4), '-djpeg', sprintf('%s/t_%03d', FramesFolder,i) );
%     SavePdfFast(sprintf('%s/Constnt_%02d', FramesFolder,frame));
%     F(i) = getframe(gcf);
	% Write this frame out to a new video file.
%   	writeVideo(writerObj, F(i));
	

end

%%
% 1;
% 
% f=.2:.2:1;
% 
% h1=figure('name','fixed velocity');
% 
% for i=1:length(f)
% 
%     Time=linspace(0,1/2/f(i),length(L_ThetaS));
% 
%     % t=0 --> theta =0 & omega= 0
%     % t=T/2 --> theta = middle of range & omega= max
%     % t=T --> theta = end of range & omega= 0
%     % a=-ThetaS(end)/2;
%     % b=pi/Time(end);
%     % c=pi/2;
%     % d=-a;
%     % 
%     % ThetaTime=a*sin(b*Time+c)+d;    % position
%     % OmegaTime=a*b*cos(b*Time+c);    % angular velicty
%     % AlphaTime=-a*b^2*sin(b*Time+c);    % angular acceleration
% 
%     a=ThetaS(end)/Time(end);
%     ThetaTime=a*Time;    % position
%     OmegaTime=a*ones(size(Time));    % angular velicty
%     AlphaTime=0*ones(size(Time));    % angular acceleration
% 
%     a=2*(ThetaS(end))/((Time(end)^2));
%     ThetaTime=1/2*a*Time.^2;    % position
%     OmegaTime=a*(Time);    % angular velicty
%     AlphaTime=a*ones(size(Time));    % angular acceleration
% 
% %     figure
% %     plot(Time,ThetaTime)
% %      hold all
% %      plot(Time,OmegaTime)
% %      plot(Time,AlphaTime)
% 
%     L_ThetaSTime=interp1(ThetaS,smooth(L_ThetaS),ThetaTime);
%     DL_ThetaSTime=differential(smooth(L_ThetaSTime,3)',Time,diff(Time(1:2)));
%     D2L_ThetaSTime=differential(smooth(DL_ThetaSTime,3)',Time,diff(Time(1:2)));
%     r_thetaRTime=interp1(ThetaS,r_thetaR,ThetaTime);
% 
% 
%     ResForce=m*3/2*(DL_ThetaSTime.*D2L_ThetaSTime  +  L_ThetaSTime.^2.*OmegaTime.*AlphaTime  +  L_ThetaSTime.*DL_ThetaSTime.*OmegaTime.^2);
%     tauTime=(m*3/2*(DL_ThetaSTime.*D2L_ThetaSTime  +  L_ThetaSTime.^2.*OmegaTime.*AlphaTime  +  L_ThetaSTime.*DL_ThetaSTime.*OmegaTime.^2)+K*(L_ThetaSTime-l0).*DL_ThetaSTime)./OmegaTime;
% 
% 
%     Phi_Alpha_Time=interp1(ThetaS,smooth(Phi_Alpha),ThetaTime);
% 
%     D_Omega= ( (D2L_ThetaSTime  -L_ThetaSTime.*OmegaTime.^2).*sin(Phi_Alpha_Time) + (L_ThetaSTime.*AlphaTime +2 *DL_ThetaSTime.*OmegaTime).*cos(Phi_Alpha_Time))/R;
%     Ff=m*R*D_Omega/2;
% 
%     figure(h1)
%     subplot(2,2,1)
%     plot(ThetaTime,ResForce)
%     xlabel('\theta_s (deg)')
%     ylabel('Res Force (N)')
%     grid on
%     xlim([ThetaTime(1) ThetaTime(end)])
%     hold all
% 
%     subplot(2,2,2)
%     plot(ThetaTime,D_Omega)
%     xlabel('\theta_s (deg)')
%     ylabel('{\alpha} (rad/s^2)')
%     grid on
%     xlim([ThetaTime(1) ThetaTime(end)])
%     hold all
% 
%     subplot(2,2,3)
%     plot(ThetaTime,Ff)
%     xlabel('\theta_s (deg)')
%     ylabel('f_f (N)')
%     grid on
%     xlim([ThetaTime(1) ThetaTime(end)])
%     hold all
% 
%     subplot(2,2,4)
%     plot(ThetaTime,tauTime)
%     xlabel('\theta_s (deg)')
%     ylabel('tau (N/m)')
%     grid on
%     hold all
%     plot(ThetaS,tau)
%     xlim([ThetaTime(1) ThetaTime(end)])
%     hold all
% end
% %%
% 
% f=0.5:0.5:5;
% 
% 
% % plot(ThetaS,tau)
% % xlabel('\theta_s (deg)')
% % ylabel('tau (N/m)')
% % grid on
% % hold all
% h4=figure('name','fixed acceleration');
% 
% for i=1:length(f);
% 
%     Time=linspace(0,1/2/f(i),200);    
% 
%     a=2*(ThetaS(end))/((Time(end)^2));
%     ThetaTime=1/2*a*Time.^2;    % position
%     OmegaTime=a*(Time);    % angular velicty
%     AlphaTime=a*ones(size(Time));    % angular acceleration
% 
% %     a=ThetaS(end)/Time(end);
% %     ThetaTime=a*Time;    % position
% %     OmegaTime=a*ones(size(Time));    % angular velicty
% %     AlphaTime=0*ones(size(Time));    % angular acceleration
% 
%     L_ThetaSTime=interp1(ThetaS,smooth(L_ThetaS),ThetaTime);
%     DL_ThetaSTime=differential(smooth(L_ThetaSTime,3)',Time,diff(Time(1:2)));
%     D2L_ThetaSTime=differential(smooth(DL_ThetaSTime,3)',Time,diff(Time(1:2)));
% 
% 
%     ResForce=m*3/2*(DL_ThetaSTime.*D2L_ThetaSTime  +  L_ThetaSTime.^2.*OmegaTime.*AlphaTime  +  L_ThetaSTime.*DL_ThetaSTime.*OmegaTime.^2);
%     tauTime=(m*3/2*(DL_ThetaSTime.*D2L_ThetaSTime  +  L_ThetaSTime.^2.*OmegaTime.*AlphaTime  +  L_ThetaSTime.*DL_ThetaSTime.*OmegaTime.^2)+K*(L_ThetaSTime-l0).*DL_ThetaSTime)./OmegaTime;
% 
% 
%     Phi_Alpha_Time=interp1(ThetaS,smooth(Phi_Alpha),ThetaTime);
% 
%     D_Omega= ( (D2L_ThetaSTime  -L_ThetaSTime.*OmegaTime.^2).*sin(Phi_Alpha_Time) + (L_ThetaSTime.*AlphaTime +2 *DL_ThetaSTime.*OmegaTime).*cos(Phi_Alpha_Time))/R;
%     Ff=m*R*D_Omega/2;
% 
%     figure(h4)
%     subplot(2,2,1)
%     plot(ThetaTime,ResForce)
%     xlabel('\theta_s (deg)')
%     ylabel('Res Force (N)')
%     grid on
%     xlim([ThetaTime(1) ThetaTime(end)])
%     hold all
% 
%     subplot(2,2,2)
%     plot(ThetaTime,D_Omega)
%     xlabel('\theta_s (deg)')
%     ylabel('{\alpha} (rad/s^2)')
%     grid on
%     xlim([ThetaTime(1) ThetaTime(end)])
%     hold all
% 
%     subplot(2,2,3)
%     plot(ThetaTime,Ff)
%     xlabel('\theta_s (deg)')
%     ylabel('f_f (N)')
%     grid on
%     xlim([ThetaTime(1) ThetaTime(end)])
%     hold all
% 
%     subplot(2,2,4)
%     plot(ThetaTime,tauTime)
%     xlabel('\theta_s (deg)')
%     ylabel('tau (N/m)')
%     grid on
%     hold all
%     plot(ThetaS,tau)
%     xlim([ThetaTime(1) ThetaTime(end)])
%     hold all
% 
% 
% end
% 
% 
% %%
% 
% f=1/9:2/9:12/9;
% 
% 
% h4=figure('name','Siniude acceleration');
% 
% for i=1:length(f);
% 
%     Time=linspace(0,1/2/f(i),200);    
% 
%     a=-ThetaS(end)/2;
%     b=pi/Time(end);
%     c=pi/2;
%     d=-a;
%     ThetaTime=a*sin(b*Time+c)+d;    % position
%     OmegaTime=a*b*cos(b*Time+c);    % angular velicty
%     AlphaTime=-a*b^2*sin(b*Time+c);    % angular acceleration
% 
% 
% %     a=ThetaS(end)/Time(end);
% %     ThetaTime=a*Time;    % position
% %     OmegaTime=a*ones(size(Time));    % angular velicty
% %     AlphaTime=0*ones(size(Time));    % angular acceleration
% 
%     L_ThetaSTime=interp1(ThetaS,smooth(L_ThetaS),ThetaTime);
%     DL_ThetaSTime=differential(smooth(L_ThetaSTime,3)',Time,diff(Time(1:2)));
%     D2L_ThetaSTime=differential(smooth(DL_ThetaSTime,3)',Time,diff(Time(1:2)));
% 
% 
%     ResForce=m*3/2*(DL_ThetaSTime.*D2L_ThetaSTime  +  L_ThetaSTime.^2.*OmegaTime.*AlphaTime  +  L_ThetaSTime.*DL_ThetaSTime.*OmegaTime.^2);
%     tauTime=(m*3/2*(DL_ThetaSTime.*D2L_ThetaSTime  +  L_ThetaSTime.^2.*OmegaTime.*AlphaTime  +  L_ThetaSTime.*DL_ThetaSTime.*OmegaTime.^2)+K*(L_ThetaSTime-l0).*DL_ThetaSTime)./OmegaTime;
% 
% 
%     Phi_Alpha_Time=interp1(ThetaS,smooth(Phi_Alpha),ThetaTime);
% 
%     D_Omega= ( (D2L_ThetaSTime  -L_ThetaSTime.*OmegaTime.^2).*sin(Phi_Alpha_Time) + (L_ThetaSTime.*AlphaTime +2 *DL_ThetaSTime.*OmegaTime).*cos(Phi_Alpha_Time))/R;
%     Ff=m*R*D_Omega/2;
% 
%     figure(h4)
%     subplot(2,2,1)
%     plot(ThetaTime,ResForce)
%     xlabel('\theta_s (deg)')
%     ylabel('Res Force (N)')
%     grid on
%     xlim([ThetaTime(1) ThetaTime(end)])
%     hold all
% 
%     subplot(2,2,2)
%     plot(ThetaTime,D_Omega)
%     xlabel('\theta_s (deg)')
%     ylabel('{\alpha} (rad/s^2)')
%     grid on
%     xlim([ThetaTime(1) ThetaTime(end)])
%     hold all
% 
%     subplot(2,2,3)
%     plot(ThetaTime,Ff)
%     xlabel('\theta_s (deg)')
%     ylabel('f_f (N)')
%     grid on
%     xlim([ThetaTime(1) ThetaTime(end)])
%     hold all
% 
%     subplot(2,2,4)
%     plot(ThetaTime(3:end-2),tauTime(3:end-2))
%     xlabel('\theta_s (deg)')
%     ylabel('tau (N/m)')
%     grid on
%     hold all
%     plot(ThetaS,tau)
%     xlim([ThetaTime(1) ThetaTime(end)])
%     hold all
% 
% 
% end
% 
% %% ShowTime for paper
% 
% f=2/9:2/9:12/9;
% 
% 
% h4=figure('name','Siniude acceleration');
% h5=figure('name','Siniude acceleration');
% Marker={{'o'},{'square'},{'diamond'},{'v'},{'*'},{'>'}};    
% Color=[      0         0    1.0000
%              0    0.5000         0
%         1.0000         0         0
%              0    0.7500    0.7500
%         0.7500         0    0.7500
%         0.7500    0.7500         0
%         0.2500    0.2500    0.2500];
% 
%     
%     figure(h4)
%     subplot(2,1,1)
%     plot(rad2deg(ThetaS),tau)
%     hold all
%         
% for i=1:length(f)-1;
% 
%     Time=linspace(0,1/2/f(i),200);    
% 
%     a=-ThetaS(end)/2;
%     b=pi/Time(end);
%     c=pi/2;
%     d=-a;
%     ThetaTime=a*sin(b*Time+c)+d;    % position
%     OmegaTime=a*b*cos(b*Time+c);    % angular velicty
%     AlphaTime=-a*b^2*sin(b*Time+c);    % angular acceleration
% 
% 
%     L_ThetaSTime=interp1(ThetaS,smooth(L_ThetaS),ThetaTime);
%     DL_ThetaSTime=differential(smooth(L_ThetaSTime,3)',Time,diff(Time(1:2)));
%     D2L_ThetaSTime=differential(smooth(DL_ThetaSTime,3)',Time,diff(Time(1:2)));
% 
% 
%     tauTime=(m*3/2*(DL_ThetaSTime.*D2L_ThetaSTime  +  L_ThetaSTime.^2.*OmegaTime.*AlphaTime  +  L_ThetaSTime.*DL_ThetaSTime.*OmegaTime.^2)+K*(L_ThetaSTime-l0).*DL_ThetaSTime)./OmegaTime;
% 
%     figure(h4)
%     subplot(2,1,1)
%     plot(rad2deg(ThetaTime(3:end-2)),tauTime(3:end-2))
%     ylim(([min(min(tau)*[1.1 .91]), max(max(tau)*[1.1 .91])]))
%     xlim(rad2deg([ThetaS(1) ThetaS(end)]))
%     xlabel('\theta_s (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     ylabel('\tau (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     grid on
%     set(gca,'FontSize',10,'FontWeight','bold','FontName','mwa_cmb10')
%     set(gca,'XTickLabel',{'0','','','90','','','180','','','270'},...
%             'XTick',[0 30 60 90 120 150 180 210 240 270])
%     hold all
% 
%     subplot(2,1,2)
%     plot(rad2deg(ThetaTime(30:20:end-30)),((interp1(ThetaS,tau,ThetaTime(30:20:end-30))-tauTime(30:20:end-30))),...
%             'linewidth',2,'DisplayName',['\omega=',num2str(270*2*f(i)/360*60)],...
%             'marker',Marker{i}{1},'color',Color(i,:),'markersize',8)
%     ylim(([-.025, .025]))
%     xlim(rad2deg([ThetaS(1) ThetaS(end)]))
%     xlabel('\theta_s (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     ylabel('Deviation (N.m)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     grid on
%     set(gca,'FontSize',10,'FontWeight','bold','FontName','mwa_cmb10')
%     set(gca,'XTickLabel',{'0','','','90','','','180','','','270'},...
%             'XTick',[0 30 60 90 120 150 180 210 240 270])
%     hold all
% 
%     
%     Phi_Alpha_Time=interp1(ThetaS,smooth(Phi_Alpha),ThetaTime);
% 
%     D_Omega= ( (D2L_ThetaSTime  -L_ThetaSTime.*OmegaTime.^2).*sin(Phi_Alpha_Time) + (L_ThetaSTime.*AlphaTime +2 *DL_ThetaSTime.*OmegaTime).*cos(Phi_Alpha_Time))/R;
%     Ff=m*R*D_Omega/2;
% 
%     
%     figure(h5)
%     plot(rad2deg(ThetaTime(20:20:end-20)),Ff(20:20:end-20),...
%             'linewidth',2,'DisplayName',['\omega=',num2str(270*2*f(i)/360*60)],...
%             'marker',Marker{i}{1},'color',Color(i,:),'markersize',8)
%     grid on
%     ylim(([-.2 .2]))
%     xlim(rad2deg([ThetaS(1) ThetaS(end)]))
%     xlabel('\theta_s (deg)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     ylabel('f_f (N)','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%     grid on
%     set(gca,'FontSize',10,'FontWeight','bold','FontName','mwa_cmb10')
%     set(gca,'XTickLabel',{'0','','','90','','','180','','','270'},...
%             'XTick',[0 30 60 90 120 150 180 210 240 270])
%     hold all
%     
% 
% 
% end
% figure(h4)
% gg=legend(gca,'show');
% set(gg,'Orientation','horizontal',...
%     'Location','Best');
% figure(h5)
% gg=legend(gca,'show');
% set(gg,'Orientation','horizontal',...
%     'Location','Best');
% %%
% for i=1:length(f)-1;
% 
%     Time=linspace(0,1/2/f(i),200);    
% 
%     a=-ThetaS(end)/2;
%     b=pi/Time(end);
%     c=pi/2;
%     d=-a;
%     ThetaTime=a*sin(b*Time+c)+d;    % position
%     OmegaTime=a*b*cos(b*Time+c);    % angular velicty
%     AlphaTime=-a*b^2*sin(b*Time+c);    % angular acceleration
% 
% 
%     L_ThetaSTime=interp1(ThetaS,smooth(L_ThetaS),ThetaTime);
%     DL_ThetaSTime=differential(smooth(L_ThetaSTime,3)',Time,diff(Time(1:2)));
%     D2L_ThetaSTime=differential(smooth(DL_ThetaSTime,3)',Time,diff(Time(1:2)));
% 
% 
%     tauTime=(m*3/2*(DL_ThetaSTime.*D2L_ThetaSTime  +  L_ThetaSTime.^2.*OmegaTime.*AlphaTime  +  L_ThetaSTime.*DL_ThetaSTime.*OmegaTime.^2)+K*(L_ThetaSTime-l0).*DL_ThetaSTime)./OmegaTime;
% 
%     figure(h4)
%     subplot(2,1,2)
%     plot(rad2deg(ThetaTime(3:end-2)),((interp1(ThetaS,tau,ThetaTime(3:end-2))-tauTime(3:end-2))),...
%         'linewidth',2,'color',Color(i,:))
%     ylim(([-.025, .025]))
%     xlim(rad2deg([ThetaS(1) ThetaS(end)]))
%     grid on
%     set(gca,'FontSize',10,'FontWeight','bold','FontName','mwa_cmb10')
%     set(gca,'XTickLabel',{'0','','','90','','','180','','','270'},...
%             'XTick',[0 30 60 90 120 150 180 210 240 270])
%     hold all
% 
%     Phi_Alpha_Time=interp1(ThetaS,smooth(Phi_Alpha),ThetaTime);
% 
%     D_Omega= ( (D2L_ThetaSTime  -L_ThetaSTime.*OmegaTime.^2).*sin(Phi_Alpha_Time) + (L_ThetaSTime.*AlphaTime +2 *DL_ThetaSTime.*OmegaTime).*cos(Phi_Alpha_Time))/R;
%     Ff=m*R*D_Omega/2;
% 
%     
%     figure(h5)
%     plot(rad2deg(ThetaTime),Ff,...
%             'linewidth',2,'color',Color(i,:))
%     
% 
% end

end