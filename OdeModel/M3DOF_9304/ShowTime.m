function [TorqueDesire,Q,D1Q,D2Q,BetaOptimal,IntU2,IntUdq,IntAbsUdq,IntAbsUdqDesire,CostSlope,RMS]=...
                ShowTime(Alpha,Time,Tres,Degree,Weight,Landa,QQ,B,Xef,Yef,m,L,g,MinSinValue,ShowFlag,Period,Mode,Name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% close all

% ShowFlag =  Show   or  DntShow
% Period   =  1cycle or  2Cycle
% Mode     =  CostA  or  CastB
%% Regerate Trjectory from Coef

if(strcmp(Mode,'CostA'))
    rQ=Degree(1);
    SubplotNUM=1;
elseif(strcmp(Mode,'CostB'))
    nn=Degree(1);
    rQ=Degree(2) ;
    rU=Degree(3) ;
    SubplotNUM=2;
elseif(strcmp(Mode,'CostC'))
    nn=Degree(1);
    rQ=Degree(2) ;
    rU=Degree(3) ;
    SubplotNUM=2;    
end


Alpha_Q1=Alpha(1:(rQ+1));
Alpha_Q2=Alpha((rQ+1)+1:2*(rQ+1));
Alpha_Q3=Alpha(2*(rQ+1)+1:3*(rQ+1));
Alpha_D1Q1=Alpha_Q1(1:end-1).*(rQ:-1:1);
Alpha_D1Q2=Alpha_Q2(1:end-1).*(rQ:-1:1);
Alpha_D1Q3=Alpha_Q3(1:end-1).*(rQ:-1:1);
Alpha_D2Q1=Alpha_D1Q1(1:end-1).*(rQ-1:-1:1);
Alpha_D2Q2=Alpha_D1Q2(1:end-1).*(rQ-1:-1:1);
Alpha_D2Q3=Alpha_D1Q3(1:end-1).*(rQ-1:-1:1);

Q1=polyval(Alpha_Q1,Time);
Q2=polyval(Alpha_Q2,Time);
Q3=polyval(Alpha_Q3,Time);
D1Q1=polyval(Alpha_D1Q1,Time);
D1Q2=polyval(Alpha_D1Q2,Time);
D1Q3=polyval(Alpha_D1Q3,Time);
D2Q1=polyval(Alpha_D2Q1,Time);
D2Q2=polyval(Alpha_D2Q2,Time);
D2Q3=polyval(Alpha_D2Q3,Time);

QVal=[Q1;Q2;Q3];

%% EF

% % For Circle
% Xef=A*cos(2*pi/f*time)+1.5;
% Yef=A*sin(2*pi/f*time)+1;

% % For Line
% Xef=A*cos(2*pi/(1*f)*time);
% Yef=2*ones(size(time));

Pos=[Xef;Yef];

RPos=L*[cos(Q1)+cos(Q1+Q2)+cos(Q1+Q2+Q3);
        sin(Q1)+sin(Q1+Q2)+sin(Q1+Q2+Q3)];

RMS=sqrt(sum(sum((RPos-Pos).^2))*Tres/(Time(end)-Time(1)));

if(strcmp(ShowFlag,'Show'))
    
    figure('name',['WorkSapce : ',Name])
        plot(Pos(1,:),Pos(2,:),'linewidth',2,'linestyle','-','color','b')
        title('WorkSpace','FontWeight','bold','FontName','mwa_cmb10');
        hold on
        plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','-','color','r')
        hold off
        xlabel('x (m)','fontsize',12,'FontName','mwa_cmb10')
        ylabel('y (m)','fontsize',12,'FontName','mwa_cmb10')

    
    AnimBot3DOF(Time,QVal',L);
end
%% Trajectory
if(strcmp(ShowFlag,'Show'))
    
    figure('name',['Joints trajectory : ',Name])
        subplot(3,1,1)
        plot(Time,Q1,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),Q1,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        title('Jonits Trajectory','FontWeight','bold','FontName','mwa_cmb10');
        grid on
        xlabel('time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('q_1 (rad)','fontsize',14,'FontName','mwa_cmb10');

        subplot(3,1,2)
        plot(Time,Q2,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),Q2,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        grid on
        xlabel('time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('q_2 (rad)','fontsize',14,'FontName','mwa_cmb10');

        subplot(3,1,3)
        plot(Time,Q3,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold all
            plot(Time+Time(end),Q3,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        grid on
        xlabel('time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('q_3 (rad)','fontsize',14,'FontName','mwa_cmb10');

end
%% Torque

mL1=m;
mL2=m;
mL3=m;
LL1=L;
LL2=L;
LL3=L;

TorqueDesire=zeros(3,length(Time));
IntU2=0;
IntAbsUdq=0;
IntAbsUdqDesire=[];
IntUdq=0;
CostSlope=0;

for i=1:length(Time)
    q1=Q1(i);
    q2=Q2(i);
    q3=Q3(i);
    D1q1=D1Q1(i);
    D1q2=D1Q2(i);
    D1q3=D1Q3(i);
    D2q1=D2Q1(i);
    D2q2=D2Q2(i);
    D2q3=D2Q3(i);

    
    MM=[mL3*LL1*LL3*cos(q2+q3)+LL1*LL2*cos(q2)*mL2+2*LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL3*LL1^2+LL1^2*mL2+mL1*LL1^2/3+mL2*LL2^2/3,mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6;
        mL3*LL1*LL3*cos(q2+q3)/2+LL1*LL2*cos(q2)*mL2/2+LL1*LL2*cos(q2)*mL3+mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL2*LL3*cos(q3)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6;
        mL3*LL3*(3*LL1*cos(q2+q3)+2*LL3+3*LL2*cos(q3))/6,(2*LL3+3*LL2*cos(q3))*LL3*mL3/6,mL3*LL3^2/3];

    CC= [-mL3*LL1*LL3*sin(q2+q3)*(D1q2+D1q3)/2-LL2*(LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL1*LL3*(D1q1+D1q2+D1q3)*sin(q2+q3)/2-LL2*(sin(q2)*D1q1*LL1*(mL2+2*mL3)+LL1*sin(q2)*(mL2+2*mL3)*D1q2+mL3*LL3*sin(q3)*D1q3)/2,-mL3*LL3*(LL1*sin(q2+q3)+LL2*sin(q3))*(D1q1+D1q2+D1q3)/2;
          mL3*LL1*D1q1*LL3*sin(q2+q3)/2+(sin(q2)*D1q1*LL1*(mL2+2*mL3)-mL3*LL3*sin(q3)*D1q3)*LL2/2,-mL3*LL2*LL3*sin(q3)*D1q3/2,-mL3*LL3*sin(q3)*(D1q1+D1q2+D1q3)*LL2/2;
          mL3*(LL1*D1q1*sin(q2+q3)+LL2*sin(q3)*(D1q1+D1q2))*LL3/2,LL3*LL2*sin(q3)*(D1q1+D1q2)*mL3/2,0];

    GG= [g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2+2*cos(q1)*LL1*mL3+2*cos(q1)*LL1*mL2+cos(q1)*LL1*mL1)/2;
         g*(mL3*LL3*cos(q1+q2+q3)+2*LL2*cos(q1+q2)*mL3+LL2*cos(q1+q2)*mL2)/2;
         g*mL3*LL3*cos(q1+q2+q3)/2];


    TorqueDesire(:,i) = MM*[D2q1;D2q2;D2q3] + CC*[D1q1;D1q2;D1q3] + GG;
end

if (strcmp(Mode,'CostA')) % for CF1 % for CF1 % for CF1 % for CF1

    IntU2=sum(sum(TorqueDesire.^2,2).*Weight)*Tres/sum(Weight);
    IntAbsUdq=sum(sum(abs(TorqueDesire.*[D1Q1;D1Q2;D1Q3]),2).*Weight)*Tres/sum(Weight);
    IntUdq=sum(((sum((TorqueDesire.*[D1Q1;D1Q2;D1Q3]),2)).^2).*Weight)*Tres/sum(Weight);
    for Joint=1:3
        Till=floor(size(QVal,2)/1);
        ThetaShift=QVal(Joint, 1:Till)-min(QVal(Joint, 1:Till));
        ThetaShiftScale = ThetaShift* (deg2rad(270) /  max(ThetaShift));
        tau=TorqueDesire(Joint, 1:Till)-min(TorqueDesire(Joint, 1:Till));
        tauShiftScale= tau /max(tau);
        DTa=diff(tauShiftScale)./diff(ThetaShiftScale)*0.75;
        CostSlope=CostSlope+sum((DTa*4/3).^2)*Weight(Joint)/sum(Weight);
    end
    
    BetaOptimal=[];

elseif(strcmp(Mode,'CostB'))   % for CF2
    % Omega matrix
    Omega = B * (B'*B)^-1 *diag(Weight) * (B'*B)^-1 * B';
    % Integral Matrix
    Iu=0;
    Iq=zeros(nn*(rU+1),nn*(rU+1));
    Iuq=zeros(1,(rU+1)*nn);
    for tt=1:length(Time)
        QQ_conc=zeros((rU+1)*nn,nn);% \underline{\underline{\mathcal{Q}}}^{r_u}
        for joint=1:nn
            QQ_rU_Joint=QVal(joint,tt).^(rU:-1:0)';
            QQ_conc(1+(joint-1)*(rU+1):(joint)*(rU+1),joint)= QQ_rU_Joint;
        end
        Iu = Iu + TorqueDesire(:,tt)'*Omega*TorqueDesire(:,tt);
        Iq = Iq + QQ_conc*Omega*QQ_conc';
        Iuq= Iuq+ TorqueDesire(:,tt)'*Omega*QQ_conc';
    end
    Iu=Iu*Tres;
    Iq=Iq*Tres;
    Iuq=Iuq*Tres;
    Idq_conc=zeros((rU+1)*nn,(rU+1)*nn);
    for Joint=1:nn
        c_hat =(max(QVal(Joint,:)) - min(QVal(Joint,:)))* 2 / 3/pi;
        d_hat = min(QVal(Joint,:));
        for kk=1:rU+1
            Psi(kk,:) = [zeros(1,kk-1), (c_hat^(rU+1-(kk)))* poly(-d_hat/c_hat*ones(1,rU+1-(kk)))];
        end
        Idq = Weight(Joint)*c_hat* Psi*QQ*Psi';
        Idq_conc((Joint-1)*(rU+1)+1:(Joint)*(rU+1),(Joint-1)*(rU+1)+1:(Joint)*(rU+1)) = Idq;
    end
    
    SVDsol=SVDBlockInvertor((Landa*Iq+(1-Landa)*Idq_conc),nn,rU+1,MinSinValue);
%     BetaOptimal=  Landa*(Landa*Iq + (1-Landa)*Idq_conc )^-1 *Iuq';
    BetaOptimal=  Landa*SVDsol*Iuq';
     
    IntU2=1/2*Iu+1/2*BetaOptimal'*Iq*BetaOptimal-Iuq*BetaOptimal;
    CostSlope=1/2*BetaOptimal'*Idq_conc*BetaOptimal;
%     Cost = 1/2* Landa*Iu - 1/2*Landa^2*Iuq*(Landa*Iq + (1-Landa)*Idq_conc )^-1 *Iuq'
%     IntU2*Landa+(1-Landa)*CostSlope    
%     Cost = 1/2* Landa*Iu - 1/2*Landa^2*Iuq*SVDsol *Iuq';

     
    TorquePassiveQ1valOptimal=polyval(BetaOptimal(1: rU+1),Q1);
    TorquePassiveQ2valOptimal=polyval(BetaOptimal(rU+1+1: 2*(rU+1)),Q2);
    TorquePassiveQ3valOptimal=polyval(BetaOptimal(2*(rU+1)+1:3*(rU+1)),Q3);
    TorquePassiveValOptimal=[TorquePassiveQ1valOptimal; TorquePassiveQ2valOptimal; TorquePassiveQ3valOptimal];
    TorqueActive=TorqueDesire-TorquePassiveValOptimal;
    
    IntAbsUdqDesire=(sum(abs(TorqueDesire.*[D1Q1;D1Q2;D1Q3]),2).*Weight)*Tres/sum(Weight);
    IntAbsUdq=sum(sum(abs(TorqueActive.*[D1Q1;D1Q2;D1Q3]),2).*Weight)*Tres/sum(Weight);
    
  
elseif(strcmp(Mode,'CostC'))   % for CF3
    BetaOptimal=[];
    Cost=0;
    CostSub=0;
    for i=1:nn
        QQ=[];
        DQ=[];
        CoefBLSI = LSParamPoly(QVal(i,:),TorqueDesire(i,:)',rU,Landa);    
        for j=1:length(QVal(i,:))
             QQ(j,:) = QVal(i,j).^(rU:-1:0)';
             DQ(j,:) = ([QVal(i,j).^(rU-1:-1:0) 0].*(rU:-1:0))';
        end
        
        IntU2=IntU2 + ...
              Weight(i)* 1/2*(TorqueDesire(i,:)' - QQ*CoefBLSI )'* (TorqueDesire(i,:)' - QQ*CoefBLSI )*Tres;
                          
        CostSlope = CostSlope + Weight(i)* 1/2*CoefBLSI'*(DQ'*DQ)*CoefBLSI*Tres;

        BetaOptimal=[BetaOptimal;CoefBLSI];
    end


    TorquePassiveQ1valOptimal=polyval(BetaOptimal(0*(rU+1)+1:1*(rU+1)),Q1);
    TorquePassiveQ2valOptimal=polyval(BetaOptimal(1*(rU+1)+1:2*(rU+1)),Q2);
    TorquePassiveQ3valOptimal=polyval(BetaOptimal(2*(rU+1)+1:3*(rU+1)),Q3);
    TorquePassiveValOptimal=[TorquePassiveQ1valOptimal; TorquePassiveQ2valOptimal; TorquePassiveQ3valOptimal];
    TorqueActive=TorqueDesire-TorquePassiveValOptimal;
    
    IntAbsUdqDesire=(sum(abs(TorqueDesire.*[D1Q1;D1Q2;D1Q3]),2).*Weight)*Tres/sum(Weight);
    IntAbsUdq=sum(sum(abs(TorqueActive.*[D1Q1;D1Q2;D1Q3]),2).*Weight)*Tres/sum(Weight);

end


if(strcmp(ShowFlag,'Show'))
    
    figure('name',[' Desired Torque vs Time : ',Name])
        plot(Time,TorqueDesire,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),TorqueDesire,'linewidth',2,'linestyle','-.')
            hold off
        end
        title(' Desired Torque','FontWeight','bold','FontName','mwa_cmb10');
        legend('\tau_d_1','\tau_d_2','\tau_d_3')
        xlabel('time (s)','fontsize',14,'FontName','mwa_cmb10');
        ylabel('\tau','fontsize',14,'FontName','mwa_cmb10');
        grid on
        
    figure('name',[' Desired Torque*\dot{q} vs time : ',Name])
        subplot(3,1,1,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(Time,TorqueDesire(1,:).*D1Q1,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),TorqueDesire(1,:).*D1Q1,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        title('${\tau * \dot q}$ vs Time','FontWeight','bold', 'interpreter','latex','fontsize',18)
        xlabel('time','fontsize',12,'FontName','mwa_cmb10');
        ylabel('${\tau_d_1 * \dot q_1}$', 'interpreter','latex','fontsize',12,'FontName','mwa_cmb10');
        grid on
        set(gca,'YMinorGrid','on')


        subplot(3,1,2,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(Time,TorqueDesire(2,:).*D1Q2,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),TorqueDesire(2,:).*D1Q2,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        xlabel('time','fontsize',12,'FontName','mwa_cmb10');
        ylabel('${\tau_d_2 * \dot q_2}$', 'interpreter','latex','fontsize',12,'FontName','mwa_cmb10');
        grid on
        set(gca,'YMinorGrid','on')


        subplot(3,1,3,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(Time,TorqueDesire(3,:).*D1Q3,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),TorqueDesire(3,:).*D1Q3,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        xlabel('time','fontsize',12,'FontName','mwa_cmb10');
        ylabel('${\tau_d_3 * \dot q_3}$', 'interpreter','latex','fontsize',12,'FontName','mwa_cmb10');
        grid on
        set(gca,'YMinorGrid','on')

   figure('name',['Torque vs Angle : ',Name])
        subplot(3,SubplotNUM,0*(SubplotNUM-1)+1,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(Q1,TorqueDesire(1,:),'linewidth',2)
        title('Desired/Passive Torque Angle Profile','FontWeight','bold','FontName','mwa_cmb10');
        xlabel('q_1 (rad)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('\tau_1','fontsize',14,'FontName','mwa_cmb10');
        hold off
        grid on
        set(gca,'YMinorGrid','on')
        if(strcmp(Mode,'CostB') || strcmp(Mode,'CostC'))
            hold on
            plot(Q1,TorquePassiveQ1valOptimal,'linewidth',2,'color','g','linestyle','-.')
            hold off
            legend('\tau_d','\tau_p','Orientation','horizontal')
            legend BOXOFF
            subplot(3,2,2)
            plot(Q1,TorqueActive(1,:),'linewidth',2,'color','r')
            title('Actuator Torque','FontName','mwa_cmb10');
            xlabel('q_1 (rad)','fontsize',12,'FontName','mwa_cmb10');
            ylabel('\tau_a_1','fontsize',14,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            subplot(3,2,1)
        end
        hold on
        plot(Q1(1),TorqueDesire(1,1),'linewidth',2,'linestyle','none','marker','*','markersize',10)
        plot(Q1(10),TorqueDesire(1,10),'linewidth',2,'linestyle','none','marker','*','markersize',6,'markeredgecolor','r')
        hold off
        
        subplot(3,SubplotNUM,1*(SubplotNUM-1)+2,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(Q2,TorqueDesire(2,:),'linewidth',2)
        grid on
        set(gca,'YMinorGrid','on')
        if(strcmp(Mode,'CostB') || strcmp(Mode,'CostC'))
            hold on
            plot(Q2,TorquePassiveQ2valOptimal,'linewidth',2,'color','g','linestyle','-.')
            hold off
            legend('\tau_d','\tau_p','Orientation','horizontal')
            legend BOXOFF
            subplot(3,2,4)
            plot(Q2,TorqueActive(2,:),'linewidth',2,'color','r')
            title('Actuator Torque','FontName','mwa_cmb10');
            xlabel('q_2 (rad)','fontsize',12,'FontName','mwa_cmb10');
            ylabel('\tau_a_2','fontsize',14,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            subplot(3,2,3)
        end
        hold on
        plot(Q2(1),TorqueDesire(2,1),'linewidth',2,'linestyle','none','marker','*','markersize',10)
        plot(Q2(10),TorqueDesire(2,10),'linewidth',2,'linestyle','none','marker','*','markersize',6,'markeredgecolor','r')
        xlabel('q_2 (rad)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('\tau_2','fontsize',14,'FontName','mwa_cmb10');
        hold off
        
        subplot(3,SubplotNUM,2*(SubplotNUM-1)+3,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
        plot(Q3,TorqueDesire(3,:),'linewidth',2)
        grid on
        set(gca,'YMinorGrid','on')
        if(strcmp(Mode,'CostB') || strcmp(Mode,'CostC'))
            hold on
            plot(Q3,TorquePassiveQ3valOptimal,'linewidth',2,'color','g','linestyle','-.')
            hold off
            legend('\tau_d','\tau_p','Orientation','horizontal')
            legend BOXOFF
            subplot(3,2,6)
            plot(Q3,TorqueActive(3,:),'linewidth',2,'color','r')
            title('Actuator Torque','FontName','mwa_cmb10');
            xlabel('q_3 (rad)','fontsize',12,'FontName','mwa_cmb10');
            ylabel('\tau_a_3','fontsize',14,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            subplot(3,2,5)
        end
        hold on
        plot(Q3(1),TorqueDesire(3,1),'linewidth',2,'linestyle','none','marker','*','markersize',10)
        plot(Q3(10),TorqueDesire(3,10),'linewidth',2,'linestyle','none','marker','*','markersize',6,'markeredgecolor','r')
        xlabel('q_3 (rad)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('\tau_3','fontsize',14,'FontName','mwa_cmb10');
        hold off

end


Q=[Q1;Q2;Q3];
D1Q=[D1Q1;D1Q2;D1Q3];
D2Q=[D2Q1;D2Q2;D2Q3];
end