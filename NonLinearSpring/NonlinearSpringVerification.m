%%
load Data
aalpha=[];
LL=[];
pphi=[];
LL(1)=lOde;
for i=1:length(Time)-1
    Phiaalpha(i)=acos((LL(i)^2-r_thetaRTime(i)^2+R^2)/2/R/LL(i))*sign(tauTime(i));
%     if(abs(tauTime(i))<1e-4)
%         Phiaalpha(i)=0;
%     end
        
%     if(LL(i)/R*sin(aalpha(i))> .001)
%         pphi(i)=asin(LL(i)/R*sin(aalpha(i)));
%     else
%         pphi(i)=(LL(i)/R*sin(aalpha(i)));
%     end
    LL(i+1)=LL(i)+diff(ThetaTime(i:i+1))*LL(i)*tan(real(Phiaalpha(i)));
end


%%
home
close all
clear

% if(nargin==0)
    ThetaS=deg2rad(0:.1:270)';
%     tau=.05*(ThetaS)+.1;
%     tau=.1*ones(size(ThetaS));
%      tau=+.2*(1-exp(-ThetaS))+.1;
    tau=-.2*(sin((ThetaS-3*pi/4)*2))-.25;
    tau=tau;

    m=.015;
    K=7500;
    R=.016*1;
    l0=.08*1/1;
    lOde=.1*1/1;
    FigName='Test';
    Xlabel='q';
% end
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
disp('here')
subplot(1,3,[1,2])
    hh1=polar(0,100*(max(L_ThetaS)+R));
    set(hh1,'linewidth',2);
    hold on
    hh=polar([ThetaR(1)*.9; ThetaR(1)*.9; ThetaR; ThetaR(end)*1.1 ;ThetaR(end)*1.1],[0; 100*r_thetaR(1); 100*r_thetaR ;100*r_thetaR(end); 0]);
    set(hh,'linewidth',3);
    
    
    hold all
    hhh=polar(ThetaS,100*L_ThetaS);
    set(hhh,'linestyle','-.');
    viscircles(100*L_ThetaS(1)*[cos(ThetaS(1)),sin(ThetaS(1)) ],100*R,'EdgeColor',[0.31 0.31 0.3]);
    plot([0 100*(L_ThetaS(1)+R)*cos(ThetaS(1))],[0,100*(L_ThetaS(1)+R)*sin(ThetaS(1))],'color','r','linewidth',2)
    plot([100*L_ThetaS(1)*cos(ThetaS(1)) 100*(r_thetaR(1))*cos(ThetaR(1))],[100*L_ThetaS(1)*sin(ThetaS(1)),100*(r_thetaR(1))*sin(ThetaR(1))],'color','g')
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
    plot(rad2deg(ThetaS),tau,'linewidth',2)
    plot(rad2deg(ThetaS(1)),tau(1),'linestyle','none','marker','*','markersize',8);
    xlabel([Xlabel,' (deg)'],'FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    ylabel('u (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
    grid on
    set(gca,'FontSize',18)
    legend(FigName)



%%
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
% tauTime=interp1(ThetaS,tau,ThetaTime);


%%
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
