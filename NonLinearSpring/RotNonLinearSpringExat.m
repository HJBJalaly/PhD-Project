function RotNonLinearSpringExat(NumMemShpAction,tau,K,R,l0,lOde,FigName,Xlabel)



%%
clear all
home
% close all

% if(nargin==0)
    f=1;
    m=.015;
    K=1000*1;
    R=.016*1;
    l0=.08*1/1;
    lOde=.1*1/1;
    FigName='Test';
    Xlabel='q';
%%    Kinematic Specifications
Time=linspace(0,2/2/f,10000);
    
% t=0 --> theta =0 & omega= 0
% t=T/2 --> theta = middle of range & omega= max
% t=T --> theta = end of range & omega= 0

%     % Theta= [0 : 3/2*pi]
%     a=-(3/4*pi);
%     b=pi/Time(end);
%     c=pi/2;
%     d=-a;
%     ThetaTime=a*sin(b*Time+c)+d;    % position
%     OmegaTime=a*b*cos(b*Time+c);    % angular velicty
%     AlphaTime=-a*b^2*sin(b*Time+c);    % angular acceleration

%     a=(3/2*pi)/Time(end);
%     ThetaTime=a*Time;    % position
%     OmegaTime=a*ones(size(Time));    % angular velicty
%     AlphaTime=0*ones(size(Time));    % angular acceleration

    a=2*(3/2*pi)/((Time(end)^2));
    ThetaTime=1/2*a*Time.^2;    % position
    OmegaTime=a*(Time);    % angular velicty
    AlphaTime=a*ones(size(Time));    % angular acceleration
%%    Torque Profile
%     ThetaS=deg2rad(0:.05:270)';
%       tau=.05*(ThetaTime(1:end/2))+.1;
%       tau=[tau -.05*(ThetaTime(end/2+1:end))-.1];
%      tau=-tau*.01;
       tau=-.0005*ones(size(ThetaTime));
%      tau=2*(1-exp(-ThetaTime))+1;
%       tau=.1*(sin((ThetaTime-3*pi/4)*2));

%% 
figure
subplot(4,1,1)
plot(Time,rad2deg(ThetaTime))
xlabel('time (s)')
ylabel('\theta_s (deg)')
grid on
subplot(4,1,2)
plot(Time,rad2deg(OmegaTime))
xlabel('time (s)')
ylabel('\omega_s (deg/s)')
grid on
subplot(4,1,3)
plot(Time,rad2deg(AlphaTime))
xlabel('time (s)')
ylabel('\alpha_s (deg/s^2)')
grid on
subplot(4,1,4)
plot(rad2deg(ThetaTime),tau)
xlabel('\theta_s (deg)')
ylabel('\tau (N/m)')
grid on

% end
%% Matlab Ode (variable time step)
OdeOpt= odeset('RelTol',1e-6,'AbsTol',1e-6);
% L0DL=.01;% for linear
L0DL=-.02;% for constant
% L0DL=.02;% for exp
% L0DL=.5;% for sin
[time,Xs]=ode15s(@(t,l)OdeSolverNonLinTRotExact(t,l,tau,OmegaTime,AlphaTime,Time,K,l0,m),Time,[lOde L0DL],OdeOpt);
L_ThetaS=Xs(:,1)';
DL_ThetaS=Xs(:,2)';
disp('here')
%% Manual ODE ( Fixed Step)

DeltaT=diff(Time(1:2));
ODETime=Time(1):DeltaT:Time(end);
x1(1)=lOde;
x2(1)=L0DL;
for i=1:length(Time)-1
   
%     Tau     =interp1(Time,tau,ODETime(i-1));
%     DTheta  =interp1(Time,OmegaTime,ODETime(i-1));
%     D2Theta =interp1(Time,AlphaTime,ODETime(i-1));
    Tau     =tau(i);
    DTheta  =OmegaTime(i);
    D2Theta =AlphaTime(i);


    x1(i+1)=x1(i) + x2(i)*DeltaT;
    x2(i+1)=x2(i) +...
        ((2/3/m*(-K*(x1(i)-l0)*x2(i) + Tau* DTheta)       -  x1(i)^2*DTheta*D2Theta - x1(i)*x2(i)*DTheta^2)    /x2(i))*DeltaT;

 
end
% L_ThetaS =interp1(ODETime,L_ode,Time);
% DL_ThetaS=interp1(ODETime,DL_ode,Time);
L_ThetaS =x1;
DL_ThetaS=x2;
disp('there')
%%

% Tan_Phi_Alpha=(tau./(K.*L_ThetaS.*(L_ThetaS-l0)));
for i=1:length(time)
    if(OmegaTime(i)<1e-5)
        Tan_Phi_Alpha(i)=(tau(i)./(K.*L_ThetaS(i).*(L_ThetaS(i)-l0))); % approximation
    else
        Tan_Phi_Alpha(i)=(DL_ThetaS(i)./(L_ThetaS(i).*OmegaTime(i)));
    end
end
% Phi_Alpha=atan(tau./(K.*L_ThetaS.*(L_ThetaS-l0)));
Phi_Alpha=atan(Tan_Phi_Alpha);
Cot_Phi=(L_ThetaS.*cos(Phi_Alpha)-R)./(L_ThetaS.*sin(Phi_Alpha));
Phi=acot(Cot_Phi);
r_thetaR=sin(Phi_Alpha)./sin(Phi).*L_ThetaS;
Sin_Alpha=R*(sin(Phi)./L_ThetaS);
Alpha=asin(Sin_Alpha);
ThetaR=ThetaTime(1:length(time))+Alpha;

    %% Show Time
figure
hp=polar([ThetaR' ;ThetaR(1)],[r_thetaR' ;r_thetaR(1)]);
patch( get(hp,'XData'), get(hp,'YData'), 'g')

figure
subplot(1,3,[1,2])
    hh1=polar(0,100*(max(L_ThetaS)+R));
    set(hh1,'linewidth',2);
    hold on
    hh=polar([ThetaR(1)*.9; ThetaR(1)*.9; ThetaR'; ThetaR(end)*1.1 ;ThetaR(end)*1.1],[100*r_thetaR(1); 100*r_thetaR(1); 100*r_thetaR' ;100*r_thetaR(end); 100*r_thetaR(end)]);
    set(hh,'linewidth',3);
        
    hold all
    hhh=polar(ThetaTime,100*L_ThetaS);
    set(hhh,'linestyle','-.');
    viscircles(100*L_ThetaS(1)*[cos(ThetaTime(1)),sin(ThetaTime(1)) ],100*R,'EdgeColor',[0.31 0.31 0.3]);
    plot([0 100*(L_ThetaS(1)+R)*cos(ThetaTime(1))],[0,100*(L_ThetaS(1)+R)*sin(ThetaTime(1))],'color','r','linewidth',2)
    plot([100*L_ThetaS(1)*cos(ThetaTime(1)) 100*(r_thetaR(1))*cos(ThetaR(1))],[100*L_ThetaS(1)*sin(ThetaTime(1)),100*(r_thetaR(1))*sin(ThetaR(1))],'color','g')
    plot([0 100*(r_thetaR(1)+2*R)*cos(ThetaR(1))],[0,100*(r_thetaR(1)+2*R)*sin(ThetaR(1))],'color','m')
    th = findall(gca,'Type','text');
    for i = 1:length(th),
      set(th(i),'FontSize',18)
    end

    
hss=subplot(1,3,3);
    p = get(hss, 'pos');
    p(3) = p(3) + 0.075;
    set(hss, 'pos', p);
    hold on
%     title('\tau')
    plot(rad2deg(ThetaR),tau,'linewidth',2)
    StarTau=plot(rad2deg(ThetaTime(1)),tau(1),'linestyle','none','marker','*','markersize',8);
    xlabel([Xlabel,' (deg)'],'FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    grid on
    set(gca,'FontSize',18)
    legend(FigName)


figure('units','normalized','outerposition',[0 0 1 1])

subplot(3,3,3)
xlim(rad2deg([min(ThetaTime) max(ThetaTime)]))
ylim(K*([min(L_ThetaS) max(L_ThetaS)]-l0))
hold on
title('F_s')
F_s=K*(L_ThetaS-l0);
plot(rad2deg(ThetaTime),F_s,'linewidth',2)
StarFs=plot(rad2deg(ThetaTime(1)),K*(L_ThetaS(1)-l0),'linestyle','none','marker','*','markersize',8);
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
ylabel('F_s (N)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
grid on
hold off

subplot(3,3,6)
xlim(rad2deg([min(ThetaTime) max(ThetaTime)]))
ylim(K*([min(L_ThetaS) max(L_ThetaS)]-l0))
hold on
title('N')
% N= F_s.*cos(Phi_Alpha)+tau./L_ThetaS.*sin(Phi_Alpha)+m*(  )
plot(rad2deg(ThetaTime),K*(L_ThetaS-l0)./cos(Phi_Alpha),'linewidth',2)
StarN=plot(rad2deg(ThetaTime(1)),K*(L_ThetaS(1)-l0)./cos(Phi_Alpha(1)),'linestyle','none','marker','*','markersize',8);
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
ylabel('N (N)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
grid on
hold off

subplot(3,3,9)
xlim(rad2deg([min(ThetaTime) max(ThetaTime)]))
ylim(([min(tau)-.1 max(tau)+.1]))
hold on
title('\tau')
plot(rad2deg(ThetaTime),tau,'linewidth',2)
StarTau=plot(rad2deg(ThetaTime(1)),tau(1),'linestyle','none','marker','*','markersize',8);
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
ylabel('\tau (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
grid on
hold off

subplot(3,3,7)
xlim(rad2deg([min(ThetaTime) max(ThetaTime)]))
ylim(([min([r_thetaR';L_ThetaS']) max([r_thetaR';L_ThetaS'])])*100)
hold on
plot(rad2deg(ThetaTime),100*(r_thetaR),'linewidth',2)
plot(rad2deg(ThetaTime),100*(L_ThetaS),'r','linewidth',2)
legend('r','L','Orientation','horizontal','location','best')
StarR=plot(rad2deg(ThetaTime(1)),r_thetaR(1)*100,'linestyle','none','marker','*','markersize',8);
StarL=plot(rad2deg(ThetaTime(1)),(L_ThetaS(1))*100,'r','linestyle','none','marker','*','markersize',8);
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
ylabel(' (cm)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
grid on
hold off

subplot(3,3,8)
xlim(rad2deg([min(ThetaTime) max(ThetaTime)]))
ylim(rad2deg(([min([Phi';  Phi'-Alpha']) max([Phi';  Phi'-Alpha'])])))
hold on
plot(rad2deg(ThetaTime),rad2deg(Phi),'linewidth',2)
plot(rad2deg(ThetaTime),rad2deg(Phi-Alpha),'r','linewidth',2)
legend('\phi','\phi-\alpha','Orientation','horizontal')
StarPhi=plot(rad2deg(ThetaTime(1)),rad2deg(Phi(1)),'linestyle','none','marker','*','markersize',8);
StarPhiAlpha=plot(rad2deg(ThetaTime(1)),rad2deg(Phi(1)-Alpha(1)),'r','linestyle','none','marker','*','markersize',8);
xlabel('\theta_s (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
ylabel('(deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
grid on
hold off

StepT=floor(length(ThetaTime)/20);
for i=1:StepT:length(ThetaTime)
    
    subplot(3,3,[1,2,4,5])
    hh1=polar(0,(max(L_ThetaS)+2*R)*100);
    set(hh1,'linewidth',2)
    hold on
    hh=polar(ThetaR,r_thetaR*100);
    set(hh,'linewidth',3)
    
    hold all
    hhh=polar(ThetaTime(1:StepT:i),L_ThetaS(1:StepT:i)*100);
    set(hhh,'linestyle','-.')
    viscircles(100*L_ThetaS(i)*[cos(ThetaTime(i)),sin(ThetaTime(i)) ],100*R,'EdgeColor',[0.31 0.31 0.3]);
    plot([0 100*(L_ThetaS(i)+R)*cos(ThetaTime(i))],[0,100*(L_ThetaS(i)+R)*sin(ThetaTime(i))],'color','r','linewidth',2)
    plot([100*L_ThetaS(i)*cos(ThetaTime(i)) 100*(r_thetaR(i))*cos(ThetaR(i))],[100*L_ThetaS(i)*sin(ThetaTime(i)),100*(r_thetaR(i))*sin(ThetaR(i))],'color','g')
    plot([0 100*(r_thetaR(i)+2*R)*cos(ThetaR(i))],[0,100*(r_thetaR(i)+2*R)*sin(ThetaR(i))],'color','m')
    
%     set(gca,'position',[0.08 0.325 .6 .65])
    hold off
    
    subplot(3,3,3)
    set(StarFs, 'XData',rad2deg(ThetaTime(i)), 'YData',K*(L_ThetaS(i)-l0));
    
    subplot(3,3,6)
    set(StarN, 'XData',rad2deg(ThetaTime(i)), 'YData',K*(L_ThetaS(i)-l0)./cos(Phi_Alpha(i)));
    
    subplot(3,3,9)
    set(StarTau, 'XData',rad2deg(ThetaTime(i)), 'YData',tau(i));
    
    subplot(3,3,7)
    set(StarR, 'XData',rad2deg(ThetaTime(i)), 'YData',100*r_thetaR(i));
    set(StarL, 'XData',rad2deg(ThetaTime(i)), 'YData',100*L_ThetaS(i));
    
    subplot(3,3,8)
    set(StarPhi, 'XData',rad2deg(ThetaTime(i)), 'YData',rad2deg(Phi(i)));
    set(StarPhiAlpha, 'XData',rad2deg(ThetaTime(i)), 'YData',rad2deg(Phi(i)-Alpha(i)));
    
    drawnow
end

%%
1;



% plot(Time,ThetaTime)
% hold all
% plot(Time,OmegaTime)
% plot(Time,AlphaTime)

L_ThetaSTime=interp1(ThetaTime,smooth(L_ThetaS),ThetaTime);
DL_ThetaSTime=differential(smooth(L_ThetaSTime,3)',Time,diff(Time(1:2)));
D2L_ThetaSTime=differential(smooth(DL_ThetaSTime,3)',Time,diff(Time(1:2)));


ResForce=m*3/2*(DL_ThetaSTime.*D2L_ThetaSTime  +  L_ThetaSTime.^2.*OmegaTime.*AlphaTime  +  L_ThetaSTime.*DL_ThetaSTime.*OmegaTime.^2);


Phi_Alpha_Time=interp1(ThetaTime,smooth(Phi_Alpha),ThetaTime);

D_Omega= ( (D2L_ThetaSTime  -L_ThetaSTime.*OmegaTime.^2).*sin(Phi_Alpha_Time) + (L_ThetaSTime.*AlphaTime +2 *DL_ThetaSTime.*OmegaTime).*cos(Phi_Alpha_Time))/R;
Ff=m*R*D_Omega/2;

figure
subplot(2,2,1)
plot(time,ResForce)
xlabel('time (s)')
ylabel('Res Force (N)')
grid on
xlim([time(1) time(end)])

subplot(2,2,2)
plot(time,D_Omega)
xlabel('time (s)')
ylabel('{\alpha} (rad/s^2)')
grid on
xlim([time(1) time(end)])

subplot(2,2,3)
plot(time,Ff)
xlabel('time (s)')
ylabel('f_f (N)')
grid on
xlim([time(1) time(end)])

subplot(2,2,4)
plot(time,cumsum(D_Omega))
xlabel('time (s)')
ylabel(' \omega (N)') 
grid on
xlim([time(1) time(end)])

%%

fileID = fopen('Data.txt','w');
Data=[r_thetaR(1:2:end).*cos(ThetaR(1:2:end)), r_thetaR(1:2:end).*sin(ThetaR(1:2:end)), zeros(size(ThetaR(1:2:end)))];
fprintf(fileID,'%6.5f %6.5f %6.5f \r',Data');
fclose(fileID);


end