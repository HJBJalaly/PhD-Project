function [TorqueDesire,TorqueActive_n,TorqueActive_o,TorqueMono_n,TorqueBiceps_n,Qq,D1Qq,D2Qq, IntU2_n,IntAbsUdq_n, IntU2_o,IntAbsUdq_o,IntAbsUdqDesire,RMS]=...
                ShowTimePathOptimizer(Alpha,Beta,Theta,Time,Tres,Degree,Weight,Landa,SampleRate,Xef,Yef,m,L,g,ShowFlag,Period,Name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% close all

% ShowFlag =  Show   or  DntShow
% Period   =  1cycle or  2Cycle
% Mode     =  CostA  or  CastB
%% Regerate Trjectory from Coef

nn=Degree(1);
rQ=Degree(2);
rM=Degree(3);
rB=Degree(4);
SubplotNUM=3;    


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

Qhat1=Q1+Q2;
Qhat2=Q2+Q3;
QhatJ=[Qhat1;Qhat2];

Qq=[Q1;Q2;Q3];
D1Qq=[D1Q1;D1Q2;D1Q3];
D2Qq=[D2Q1;D2Q2;D2Q3];
%% EF


Pos=[Xef;Yef];

RPos=L*[cos(Q1)+cos(Q1+Q2)+cos(Q1+Q2+Q3);
        sin(Q1)+sin(Q1+Q2)+sin(Q1+Q2+Q3)];

ErrorPower2=  (sum((RPos-Pos).^2));  
RMS=sqrt(sum(ErrorPower2)*Tres/(Time(end)-Time(1)));

if(strcmp(ShowFlag,'Show'))
    
    figure('name',['WorkSapce : ',Name])
        plot(Pos(1,:),Pos(2,:),'linewidth',2,'linestyle','-.',...
            'Color',[0.87058824300766 0.490196079015732 0])
        title('Task Space','FontWeight','bold','FontSize',16);
        hold on
        plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','-','color','b')
        hold off
        xlabel('x (m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('y (m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        axis equal
        legend('Desired Path','Optimal Travelled Path')
    figure('name',['WorkSapce Error: ',Name])
        plot(Time,sqrt(ErrorPower2)*100,'linewidth',2)
        ylim([min(sqrt(ErrorPower2))*100, max(sqrt(ErrorPower2))*100])
        title('Task Space Error','FontWeight','bold','FontSize',16);
        xlabel('Time (s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('Error (cm)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');


    
    AnimBot3DOF(Time,Qq',L);
end
%% Trajectory
if(strcmp(ShowFlag,'Show'))
    
    figure('name',['Joints trajectory : ',Name])
        subplot(3,1,1)
        plot(Time,InRangeShifter(Q1),'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),InRangeShifter(Q1),'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        title('Jonits Trajectory','FontWeight','bold','FontName','mwa_cmb10');
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('q_1 (rad)','fontsize',14,'FontName','mwa_cmb10');

        subplot(3,1,2)
        plot(Time,InRangeShifter(Q2),'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),InRangeShifter(Q2),'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('q_2 (rad)','fontsize',14,'FontName','mwa_cmb10');

        subplot(3,1,3)
        plot(Time,InRangeShifter(Q3),'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold all
            plot(Time+Time(end),InRangeShifter(Q3),'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('q_3 (rad)','fontsize',14,'FontName','mwa_cmb10');

    figure('name',['Joints Velocity : ',Name])
        subplot(3,1,1)
        plot(Time,D1Q1,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),D1Q1,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        title('Jonits Velocity','FontWeight','bold','FontName','mwa_cmb10');
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('Dq_1 (rad/s)','fontsize',14,'FontName','mwa_cmb10');

        subplot(3,1,2)
        plot(Time,D1Q2,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),D1Q2,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('Dq_2 (rad/s)','fontsize',14,'FontName','mwa_cmb10');

        subplot(3,1,3)
        plot(Time,D1Q3,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold all
            plot(Time+Time(end),D1Q3,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('Dq_3 (rad/s)','fontsize',14,'FontName','mwa_cmb10');

end
%% Torque


mL1=m;
mL2=m;
mL3=m;
LL1=L;
LL2=L;
LL3=L;

TorqueDesire=zeros(3,length(Time));
% IntU2=0;
% IntAbsUdq=0;
% IntAbsUdqDesire=[];
% IntUdq=0;

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

rMp=rM+1;
rBp=rB+1;
% New Actuation Cost With OLD designed compliance
for i=1:nn-1
    TorqueMono_n(i,:)=polyval(Beta((i-1)*(rMp)+1:i*(rMp)),Qq(i,:));
    TorqueBiceps_n(i,:)=polyval(Theta((i-1)*(rBp)+1:(i)*(rBp)),QhatJ(i,:));
end
i=nn;
TorqueMono_n(i,:)=polyval(Beta((i-1)*(rMp)+1:i*(rMp)),Qq(i,:));

TorqueActive_n=TorqueDesire-TorqueMono_n-[TorqueBiceps_n;zeros(size(Q1))]-[zeros(size(Q1));TorqueBiceps_n];

IntU2_n=0;
for i=1:nn
    IntU2_n=IntU2_n + ...
          Weight(i)*( 1/2*(TorqueActive_n(i,:)')'*TorqueActive_n(i,:)')*Tres;
end


IntAbsUdqDesire=(sum(abs(TorqueDesire.*[D1Q1;D1Q2;D1Q3]),2))*Tres;
IntAbsUdq_n=sum(sum(abs(TorqueActive_n.*[D1Q1;D1Q2;D1Q3]),2))*Tres;


% Optimal Actuation Cost With Optimal designed compliance for new Task
[BetaOptimal_o,ThetaOptimal_o,IntU2_o,CostSlopeD2Q_o,CostParamReg_o,TorqueMonoOptimal_o,TorqueBicepsOptimal_o]=...
                OptimalParam(Qq,QhatJ,TorqueDesire,nn,rM,rB,Landa,Weight,SampleRate);
TorqueActive_o=TorqueDesire-TorqueMonoOptimal_o-[TorqueBicepsOptimal_o;zeros(size(Q1))]-[zeros(size(Q1));TorqueBicepsOptimal_o];

IntAbsUdq_o=sum(sum(abs(TorqueActive_o.*[D1Q1;D1Q2;D1Q3]),2))*Tres;


if(strcmp(ShowFlag,'Show'))
    
        Q1=InRangeShifter(Q1);
        Q2=InRangeShifter(Q2);
        Q3=InRangeShifter(Q3);

        figure('name','Torque vs Time')
            subplot(3,1,1)
            plot(Time,TorqueDesire(1,:),'linewidth',2,'linestyle','-','color','b')
            hold on
            plot(Time,TorqueMono_n(1,:),'linewidth',2,'linestyle','-.','color','g')
            plot(Time,TorqueBiceps_n(1,:),'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
            plot(Time,TorqueActive_n(1,:),'linewidth',2,'linestyle','-','color','r')
            hold off
            grid on
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            Llg=legend('u_r_1','u_m_1','u_b_1','u_a_1');
            set(Llg,'orientation','horizontal')
            xlim([Time(1) Time(end)])
            %
            subplot(3,1,2)
            plot(Time,TorqueDesire(2,:),'linewidth',2,'linestyle','-','color','b')
            hold on
            plot(Time,TorqueMono_n(2,:),'linewidth',2,'linestyle','-.','color','g')
            plot(Time,TorqueBiceps_n(1,:),'linewidth',2,'linestyle','--','Color',[0.75 0 0.75])
            plot(Time,TorqueBiceps_n(2,:),'linewidth',2,'linestyle','--','Color',[0.68 0.45 0])
            plot(Time,TorqueActive_n(2,:),'linewidth',2,'linestyle','-','color','r')
            hold off
            grid on
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            Llg=legend('u_r_2','u_m_2','u_b_1','u_b_2','u_a_2');
            set(Llg,'orientation','horizontal')
            xlim([Time(1) Time(end)])
            %
            subplot(3,1,3)
            plot(Time,TorqueDesire(3,:),'linewidth',2,'linestyle','-','color','b')
            hold on
            plot(Time,TorqueMono_n(3,:),'linewidth',2,'linestyle','-.','color','g')
            plot(Time,TorqueBiceps_n(2,:),'linewidth',2,'linestyle','--','Color',[0.68 0.45 0])
            plot(Time,TorqueActive_n(3,:),'linewidth',2,'linestyle','-','color','r')
            hold off
            grid on
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            Llg=legend('u_r_3','u_m_3','u_b_2','u_a_3');
            set(Llg,'orientation','horizontal')
            xlim([Time(1) Time(end)])
            xlabel('time (s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');


        figure('name',[' Desired Torque*\dot{q} vs time : ',Name])
            subplot(3,1,1,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            plot(Time,TorqueDesire(1,:).*D1Q1,'linewidth',2,'displayname','power_r_1')
            hold on
            plot(Time,TorqueActive_n(1,:).*D1Q1,'linewidth',2,'color','r','displayname','power_a_1')
            hold off
            legend(gca,'show','Orientation','horizontal')
            if(strcmp(Period,'2Cycle'))
                hold on
                plot(Time+Time(end),TorqueDesire(1,:).*D1Q1,'linewidth',2,'linestyle','-.')
                plot(Time+Time(end),TorqueActive_n(1,:).*D1Q1,'linewidth',2,'color','r','linestyle','-.')
                hold off
            end
            %title('${u * \dot q}$ vs Time','FontWeight','bold', 'interpreter','latex','fontsize',18)
            ylabel('${u_1 * \dot q_1}$', 'interpreter','latex','fontsize',12,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            DeltaTorque=max(TorqueDesire(1,:).*D1Q1)-min(TorqueDesire(1,:).*D1Q1);
            YLIM1=([min(TorqueDesire(1,:).*D1Q1)-0.1*DeltaTorque,max(TorqueDesire(1,:).*D1Q1)+0.1*DeltaTorque]);
            DeltaTorque=max(TorqueActive_n(1,:).*D1Q1)-min(TorqueActive_n(1,:).*D1Q1);
            YLIM2=([min(TorqueActive_n(1,:).*D1Q1)-0.1*DeltaTorque,max(TorqueActive_n(1,:).*D1Q1)+0.1*DeltaTorque]);
            ylim([min([YLIM1(1),YLIM2(1)])   max([YLIM1(2),YLIM2(2)])  ])
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            %
            subplot(3,1,2,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            plot(Time,TorqueDesire(2,:).*D1Q2,'linewidth',2,'displayname','power_r_2')
            hold on
            plot(Time,TorqueActive_n(2,:).*D1Q2,'linewidth',2,'color','r','displayname','power_a_2')
            hold off
            legend(gca,'show','Orientation','horizontal')
            if(strcmp(Period,'2Cycle'))
                hold on
                plot(Time+Time(end),TorqueDesire(2,:).*D1Q2,'linewidth',2,'linestyle','-.')
                plot(Time+Time(end),TorqueActive_n(2,:).*D1Q2,'linewidth',2,'color','r','linestyle','-.')
                hold off
            end
            ylabel('${u_2 * \dot q_2}$', 'interpreter','latex','fontsize',12,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            DeltaTorque=max(TorqueDesire(2,:).*D1Q2)-min(TorqueDesire(2,:).*D1Q2);
            YLIM1=([min(TorqueDesire(2,:).*D1Q2)-0.1*DeltaTorque,max(TorqueDesire(2,:).*D1Q2)+0.1*DeltaTorque]);
            DeltaTorque=max(TorqueActive_n(2,:).*D1Q2)-min(TorqueActive_n(2,:).*D1Q2);
            YLIM2=([min(TorqueActive_n(2,:).*D1Q2)-0.1*DeltaTorque,max(TorqueActive_n(2,:).*D1Q2)+0.1*DeltaTorque]);
            ylim([min([YLIM1(1),YLIM2(1)])   max([YLIM1(2),YLIM2(2)])  ])
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            %
            subplot(3,1,3,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            plot(Time,TorqueDesire(3,:).*D1Q3,'linewidth',2,'displayname','power_r_1')
            hold on
            plot(Time,TorqueActive_n(3,:).*D1Q3,'linewidth',2,'color','r','displayname','power_a_3')
            hold off
            legend(gca,'show','Orientation','horizontal')
            if(strcmp(Period,'2Cycle'))
                hold on
                plot(Time+Time(end),TorqueDesire(3,:).*D1Q3,'linewidth',2,'linestyle','-.')
                plot(Time+Time(end),TorqueActive_n(3,:).*D1Q3,'linewidth',2,'color','r','linestyle','-.')
                hold off
            end
            xlabel('Time','fontsize',12,'FontName','mwa_cmb10');
            ylabel('${u_3 * \dot q_3}$', 'interpreter','latex','fontsize',12,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            DeltaTorque=max(TorqueDesire(3,:).*D1Q3)-min(TorqueDesire(3,:).*D1Q3);
            YLIM1=([min(TorqueDesire(3,:).*D1Q3)-0.1*DeltaTorque,max(TorqueDesire(3,:).*D1Q3)+0.1*DeltaTorque]);
            DeltaTorque=max(TorqueActive_n(3,:).*D1Q3)-min(TorqueActive_n(3,:).*D1Q3);
            YLIM2=([min(TorqueActive_n(3,:).*D1Q3)-0.1*DeltaTorque,max(TorqueActive_n(3,:).*D1Q3)+0.1*DeltaTorque]);
            DeltaTorque=max(TorqueDesire(1,:).*D1Q1)-min(TorqueDesire(1,:).*D1Q1);
            YLIM1=([min(TorqueDesire(1,:).*D1Q1)-0.1*DeltaTorque,max(TorqueDesire(1,:).*D1Q1)+0.1*DeltaTorque]);
            DeltaTorque=max(TorqueActive_n(1,:).*D1Q1)-min(TorqueActive_n(1,:).*D1Q1);
            YLIM2=([min(TorqueActive_n(1,:).*D1Q1)-0.1*DeltaTorque,max(TorqueActive_n(1,:).*D1Q1)+0.1*DeltaTorque]);
            ylim([min([YLIM1(1),YLIM2(1)])   max([YLIM1(2),YLIM2(2)])  ])
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

    figure('name',['Torque vs Angle : ',Name])
            ap=get(gca,'position');
            subplot(3,3,1,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            plot(rad2deg( Q1),TorqueDesire(1,:),'linewidth',2)
            %title('Optimal Required and Compliance Torque-Angle Profile','FontSize',16);
            xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            hold on
            plot(rad2deg(Q1),TorqueActive_n(1,:),'linewidth',2,'color','r','linestyle','-.')
            hold off
            legend('u_r_1','u_a_1','Orientation','horizontal')
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            %
            subplot(2+(rB>0),3,2,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            plot(rad2deg( Q2),TorqueDesire(2,:),'linewidth',2)
            xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            hold on
            plot(rad2deg(Q2),TorqueActive_n(2,:),'linewidth',2,'color','r','linestyle','-.')
            hold off
            legend('u_r_2','u_a_2','Orientation','horizontal')
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            %
            subplot(3,3,3,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            plot(rad2deg(Q3),TorqueDesire(3,:),'linewidth',2)
            xlabel('q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            hold on
            plot(rad2deg(Q3),TorqueActive_n(3,:),'linewidth',2,'color','r','linestyle','-.')
            hold off
            legend('u_r_3','u_a_3','Orientation','horizontal')
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

            subplot(3,3,4,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            %title('Optimal Required and Compliance Torque-Angle Profile','FontSize',16);
            plot(rad2deg(Q1),TorqueMono_n(1,:),'linewidth',3,'color','g')
            xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            hold on
            plot(rad2deg( Q1),TorqueDesire(1,:)-TorqueBiceps_n(1,:),'linewidth',2,'linestyle','-.')
            hold off
            legend('u_m_1','u_r_1-u_b_1','Orientation','horizontal')
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            %
            subplot(3,3,5,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            plot(rad2deg(Q2),TorqueMono_n(2,:),'linewidth',3,'color','g')
            xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            hold on
            plot(rad2deg( Q2),TorqueDesire(2,:)-TorqueBiceps_n(1,:)-TorqueBiceps_n(2,:),'linewidth',2,'linestyle','-.')
            hold off
            legend('u_m_2','u_r_2-u_b_1-u_b_2','Orientation','horizontal')
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            %
            subplot(3,3,6,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            plot(rad2deg(Q3),TorqueMono_n(3,:),'linewidth',3,'color','g')
            xlabel('q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            hold on
            plot(rad2deg(Q3),TorqueDesire(3,:)-TorqueBiceps_n(2,:),'linewidth',2,'linestyle','-.')
            hold off
            legend('u_m_3','u_r_3-u_b_2','Orientation','horizontal')
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');


            if(rB>0)

                sh1=subplot(337);
                sp1=get(sh1,'position');
                set(sh1,'position',[sp1(1)+.25*(ap(3)-sp1(3)),sp1(2:end)]); 
                plot(rad2deg(Q1+Q2),2*TorqueBiceps_n(1,:),'linewidth',3,'Color',[0.75 0 0.75])
                xlabel('q_1+q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                ylabel('u_1+u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                hold on
                plot(rad2deg(Q1+Q2),TorqueDesire(1,:)+TorqueDesire(2,:)-TorqueMono_n(1,:)-TorqueMono_n(2,:)-TorqueBiceps_n(2,:),'linewidth',2,'linestyle','-.')
                hold off
                legend('2\times u_b_1','u_r_1+u_r_2-u_m_1-u_m_2-u_b_2','Orientation','horizontal')
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');

                sh2=subplot(339);
                sp2=get(sh2,'position');
                plot(rad2deg(Q2+Q3),2*TorqueBiceps_n(2,:),'linewidth',3,'Color',[0.68 0.45 0])
                set(sh2,'position',[sp2(1)-.25*(ap(3)-sp2(3)),sp2(2:end)]); 
                xlabel('q_2+q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                ylabel('u_2+u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                hold on
                plot(rad2deg(Q2+Q3),TorqueDesire(2,:)+TorqueDesire(3,:)-TorqueMono_n(2,:)-TorqueMono_n(3,:)-TorqueBiceps_n(1,:),'linewidth',2,'linestyle','-.')
                hold off
                legend('2\times u_b_2','u_r_2+u_r_3-u_m_2-u_m_3-u_b_1','Orientation','horizontal')
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');


            end


        figure('name',['Torque vs Angle : ',Name])
            ap=get(gca,'position');
            set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
            plot(rad2deg( Q1+Q2+Q3),...
                TorqueActive_n(1,:)+TorqueActive_n(2,:)+TorqueActive_n(3,:),...
                'linewidth',2)
            %title('Optimal Required and Compliance Torque-Angle Profile','FontSize',16);
            xlabel('q_1+q_2+q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            ylabel('u_a_1+u_a_2+u_a_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
            grid on
            set(gca,'YMinorGrid','on')
            hold on
                
                        

end



end