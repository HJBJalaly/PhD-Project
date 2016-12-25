% %%
% load Data
% aalpha=[];
% LL=[];
% pphi=[];
% LL(1)=lOde;
% for i=1:length(Time)-1
%     Phiaalpha(i)=acos((LL(i)^2-r_thetaRTime(i)^2+R^2)/2/R/LL(i))*sign(tauTime(i));
% %     if(abs(tauTime(i))<1e-4)
% %         Phiaalpha(i)=0;
% %     end
%         
% %     if(LL(i)/R*sin(aalpha(i))> .001)
% %         pphi(i)=asin(LL(i)/R*sin(aalpha(i)));
% %     else
% %         pphi(i)=(LL(i)/R*sin(aalpha(i)));
% %     end
%     LL(i+1)=LL(i)+diff(ThetaTime(i:i+1))*LL(i)*tan(real(Phiaalpha(i)));
% end
% 

%%  Set Parameters and Torque Profile
home
close all
clear

% Paramters
    m=.0385;
    K=5000;
    R=.015/2;
    l0=.04;
    lOde=.05;

% Profile
    ThetaS=deg2rad(0:.05:270)';
% linear
%     tau=.125*(ThetaS)-.1178*2.5; 
% constant
%      tau=.3*ones(size(ThetaS));
% exp
%     tau=(2*(1-exp((-ThetaS) )))/2;
% tanh
%      tau=.149*(tanh(-(-ThetaS+3*pi/4)*1));
% Cubic
%     tau=.02*(ThetaS-3*pi/4).^3;
% Sine
    tau=.1*(sin(2*(ThetaS-3*pi/4)));
    plot(tau)

%% Calculate the profile of Cam for the given Torque Profile
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
disp('here')

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
    xlabel(['q',' (deg)'],'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
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
    xlabel(['q',' (deg)'],'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
    set(gca,'FontSize',10,'FontWeight','bold','FontName','mwa_cmb10')
    xlim([0 270])
%     ylim([min(diff(tau(1:end))./diff(ThetaS(1:end))) max(diff(tau(1:end))./diff(ThetaS(1:end)))]*1.2)
    set(gca,'XTickLabel',{'0','','','90','','','180','','','270'},...
    'XTick',[0 30 60 90 120 150 180 210 240 270])


%% Specify the Speed Profile

f=.1;
Time=linspace(0,1/2/f,length(L_ThetaS));

% t=0 --> theta =0 & omega= 0
% t=T/2 --> theta = middle of range & omega= max
% t=T --> theta = end of range & omega= 0
% a=-ThetaS(end)/2;
% b=pi/Time(end);
% c=pi/2;
% d=-a;
% ThetaTime=a*sin(b*Time+c)+d;    % position
% OmegaTime=a*b*cos(b*Time+c);    % angular velicty
% AlphaTime=-a*b^2*sin(b*Time+c);    % angular acceleration

a=ThetaS(end)/Time(end);
ThetaTime=a*Time;    % position
OmegaTime=a*ones(size(Time));    % angular velicty
AlphaTime=0*ones(size(Time));    % angular acceleration

% a=2*(ThetaS(end))/((Time(end)^2));
% ThetaTime=1/2*a*Time.^2;    % position
% OmegaTime=a*(Time)+.01;    % angular velicty
% AlphaTime=a*ones(size(Time));    % angular acceleration

L_ThetaSTime=interp1(ThetaS,smooth(L_ThetaS),ThetaTime);
DL_ThetaSTime=differential(smooth(L_ThetaSTime,3)',Time,diff(Time(1:2)));
D2L_ThetaSTime=differential(smooth(DL_ThetaSTime,3)',Time,diff(Time(1:2)));

DL_ThetaS =differential(smooth(L_ThetaS,3)',ThetaS',diff(ThetaS(1:2)));
D2L_ThetaS=differential(smooth(DL_ThetaS,3)',ThetaS',diff(ThetaS(1:2)));

tauTime=(m*3/2*(DL_ThetaSTime.*D2L_ThetaSTime  +  L_ThetaSTime.^2.*OmegaTime.*AlphaTime  +  L_ThetaSTime.*DL_ThetaSTime.*OmegaTime.^2)+K*(L_ThetaSTime-l0).*DL_ThetaSTime)./OmegaTime;

Lextra=linspace(0,l0/5,5);
for i=1:length(Lextra)
    tauTime(end+1,:)=(m*3/2*(DL_ThetaSTime.*D2L_ThetaSTime  +  L_ThetaSTime.^2.*OmegaTime.*AlphaTime  +  L_ThetaSTime.*DL_ThetaSTime.*OmegaTime.^2)+K*(L_ThetaSTime-l0+Lextra(i)).*DL_ThetaSTime)./OmegaTime;
end

figure
    subplot(2,1,1)
        plot(ThetaTime,tauTime(1,:),'linewidth',2,'DisplayName','main')
        xlabel(['q',' (deg)'],'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ylabel(['\tau',' (N/m)'],'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        hold all
        for i=2:size(tauTime,1)
            plot(ThetaTime,tauTime(i,:),'linewidth',2,'DisplayName',['l0=',num2str(l0+Lextra(i-1))])
        end
        grid on
        hold off
        legend(gca,'show');

    subplot(2,1,2)
        xlabel(['q',' (deg)'],'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        ylabel('Sclae','FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        hold all
        for i=2:size(tauTime,1)
            plot(ThetaTime,tauTime(i,:)./tauTime(1,:),'linewidth',2,'DisplayName',['l0=',num2str(l0+Lextra(i-1))])
        end
        grid on
        hold off
        legend(gca,'show');
%% Calculate the Generated Torque 

OdeOpt= odeset('RelTol',1e-4,'AbsTol',1e-4,'Events',@(t,theta)StopCond(t,theta));
[TimeOde,Xs]=ode15s(@(t,theta)OdeSolverNonLinTRotExact_Inverser(t,theta,1.05*tauTime,ThetaS,L_ThetaS,DL_ThetaS,D2L_ThetaS,Time,K,l0,m),Time,[0 .1],OdeOpt);
ThetaOde=Xs(:,1)';
disp('there')
clf
subplot(2,1,1)
plot(TimeOde,ThetaOde)
hold all
plot(Time,ThetaTime)

subplot(2,1,2)
plot(rad2deg(ThetaOde),interp1(Time,tauTime,TimeOde))
hold all
plot(rad2deg(ThetaS),tau)
