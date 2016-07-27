function [TorqueDesire,TorqueActive,TorqueMonoOptimal,TorqueBicepsOptimal,Qq,D1Qq,D2Qq,BetaOptimal,ThetaOptimal, IntU2,IntUdq,IntAbsUdq,IntAbsUdqDesire,CostSlopeD1Q,CostSlopeD2Q,CostParamReg ,RMS]=...
                ShowTime_4R_3D(Alpha,Time,Tres,Degree,Weight,Landa,SampleRate,Sat,QQ,B,Xef,Yef,Zef,mL1,mL2,mL3,mL4,LL1,LL2,LL3,LL4,g,MinSinValue,ShowFlag,Period,Mode,Name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% close all

% ShowFlag =  Show   or  DntShow
% Period   =  1cycle or  2Cycle
% Mode     =  CostA  or  CastB
%% Regerate Trjectory from Coef

if(strcmp(Mode,'CostCc'))
    nn=Degree(1);
    rQ=Degree(2);
    rU=Degree(3);
    rB=Degree(4);
    SubplotNUM=3;    
else
    error('Wrong mode')
    return
end


Alpha_Q1=Alpha(1:(rQ+1));
Alpha_Q2=Alpha((rQ+1)+1:2*(rQ+1));
Alpha_Q3=Alpha(2*(rQ+1)+1:3*(rQ+1));
Alpha_Q4=Alpha(3*(rQ+1)+1:4*(rQ+1));
Alpha_D1Q1=Alpha_Q1(1:end-1).*(rQ:-1:1);
Alpha_D1Q2=Alpha_Q2(1:end-1).*(rQ:-1:1);
Alpha_D1Q3=Alpha_Q3(1:end-1).*(rQ:-1:1);
Alpha_D1Q4=Alpha_Q4(1:end-1).*(rQ:-1:1);
Alpha_D2Q1=Alpha_D1Q1(1:end-1).*(rQ-1:-1:1);
Alpha_D2Q2=Alpha_D1Q2(1:end-1).*(rQ-1:-1:1);
Alpha_D2Q3=Alpha_D1Q3(1:end-1).*(rQ-1:-1:1);
Alpha_D2Q4=Alpha_D1Q4(1:end-1).*(rQ-1:-1:1);

Q1=polyval(Alpha_Q1,Time);
Q2=polyval(Alpha_Q2,Time);
Q3=polyval(Alpha_Q3,Time);
Q4=polyval(Alpha_Q4,Time);
D1Q1=polyval(Alpha_D1Q1,Time);
D1Q2=polyval(Alpha_D1Q2,Time);
D1Q3=polyval(Alpha_D1Q3,Time);
D1Q4=polyval(Alpha_D1Q4,Time);
D2Q1=polyval(Alpha_D2Q1,Time);
D2Q2=polyval(Alpha_D2Q2,Time);
D2Q3=polyval(Alpha_D2Q3,Time);
D2Q4=polyval(Alpha_D2Q4,Time);

Qhat1=zeros(size(Q1));
Qhat2=Q2+Q3;
Qhat3=Q3+Q4;
QhatJ=[Qhat1;Qhat2;Qhat3];

QVal=[Q1;Q2;Q3;Q4];
Qq=[Q1;Q2;Q3;Q4];
D1Qq=[D1Q1;D1Q2;D1Q3;D1Q4];
D2Qq=[D2Q1;D2Q2;D2Q3;D2Q4];
%% EF


Pos=[Xef;Yef;Zef];

RPos=FK_RzRyRyRy_3D(Q1,Q2,Q3,Q4,LL1,LL2,LL3,LL4);

ErrorPower2=  (sum((RPos-Pos).^2));  
RMS=sqrt(sum(ErrorPower2)*Tres/(Time(end)-Time(1)));

if(strcmp(ShowFlag,'Show'))
    
    figure('name',['WorkSapce : ',Name])
        plot3(Pos(1,:),Pos(2,:),Pos(3,:),'linewidth',2,'linestyle','-.',...
            'Color',[0.87058824300766 0.490196079015732 0])
        title('Task Space','FontWeight','bold','FontSize',16);
        hold on
        plot3(RPos(1,:),RPos(2,:),RPos(3,:),'linewidth',2,'linestyle','-','color','b')
        hold off
        xlabel('x (m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('y (m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        zlabel('z (m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        axis equal
        legend('Desired Path','Optimal Travelled Path')
    figure('name',['WorkSapce Error: ',Name])
        plot(Time,sqrt(ErrorPower2)*100,'linewidth',2)
        ylim([min(sqrt(ErrorPower2))*100, max(sqrt(ErrorPower2))*100])
        title('Task Space Error','FontWeight','bold','FontSize',16);
        xlabel('Time (s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
        ylabel('Error (cm)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');


    
    AnimBot4DOF_4R_3D(Time,QVal',LL1,LL2,LL3,LL4);
end
%% Trajectory
if(strcmp(ShowFlag,'Show'))
    
    figure('name',['Joints trajectory : ',Name])
        subplot(4,1,1)
        plot(Time,Q1,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),Q1,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        title('Jonits Trajectory','FontWeight','bold','FontName','mwa_cmb10');
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('q_1 (rad)','fontsize',14,'FontName','mwa_cmb10');

        subplot(4,1,2)
        plot(Time,Q2,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),Q2,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('q_2 (rad)','fontsize',14,'FontName','mwa_cmb10');

        subplot(4,1,3)
        plot(Time,Q3,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold all
            plot(Time+Time(end),Q3,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('q_3 (rad)','fontsize',14,'FontName','mwa_cmb10');
        
        subplot(4,1,4)
        plot(Time,Q4,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold all
            plot(Time+Time(end),Q4,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('q_4 (rad)','fontsize',14,'FontName','mwa_cmb10');

    figure('name',['Joints Velocity : ',Name])
        subplot(4,1,1)
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

        subplot(4,1,2)
        plot(Time,D1Q2,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold on
            plot(Time+Time(end),D1Q2,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('Dq_2 (rad/s)','fontsize',14,'FontName','mwa_cmb10');

        subplot(4,1,3)
        plot(Time,D1Q3,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold all
            plot(Time+Time(end),D1Q3,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('Dq_3 (rad/s)','fontsize',14,'FontName','mwa_cmb10');

        subplot(4,1,4)
        plot(Time,D1Q4,'linewidth',2)
        if(strcmp(Period,'2Cycle'))
            hold all
            plot(Time+Time(end),D1Q4,'linewidth',2,'color','r','linestyle','-.')
            hold off
        end
        grid on
        xlabel('Time (s)','fontsize',12,'FontName','mwa_cmb10');
        ylabel('Dq_4 (rad/s)','fontsize',14,'FontName','mwa_cmb10');

end
%% Torque


IntU2=0;
IntAbsUdq=0;
IntAbsUdqDesire=[];
IntUdq=0;
CostSlopeD1Q=0;
CostSlopeD2Q=0;

TorqueDesire=TorqueCalculator_4R_3D(D2Qq,D1Qq,Qq,g,mL1,mL2,mL3,mL4,LL1,LL2,LL3,LL4);

if(strcmp(Mode,'CostCc'))   % for CF3c
    BetaOptimal=[];
    ThetaOptimal=[];
    
    [BetaOptimal,ThetaOptimal,IntU2,CostSlopeD2Q,CostParamReg,TorqueMonoOptimal,TorqueBicepsOptimal]=...
                    OptimalParam_4R_3D(QVal,QhatJ,TorqueDesire,nn,rU,rB,Landa,Weight,SampleRate);
    

    TorqueActive=TorqueDesire-TorqueMonoOptimal-[TorqueBicepsOptimal;zeros(size(Q1))]-[zeros(size(Q1));TorqueBicepsOptimal];

    
    IntAbsUdqDesire=(sum(abs(TorqueDesire.*[D1Q1;D1Q2;D1Q3;D1Q4]),2))*Tres;
    IntAbsUdq=sum(sum(abs(TorqueActive.*[D1Q1;D1Q2;D1Q3;D1Q4]),2))*Tres;

end


if(strcmp(ShowFlag,'Show'))
    if(strcmp(Mode,'CostCc'))
    
            Q1=InRangeShifter(Q1);
            Q2=InRangeShifter(Q2);
            Q3=InRangeShifter(Q3);
            Q4=InRangeShifter(Q4);

            figure('name','Torque vs Time')
                subplot(4,1,1)
                plot(Time,TorqueDesire(1,:),'linewidth',2,'linestyle','-','color','b')
                hold on
                plot(Time,TorqueMonoOptimal(1,:),'linewidth',2,'linestyle','-.','color','g')
                plot(Time,TorqueBicepsOptimal(1,:),'linewidth',2,'linestyle','-.','Color',[0.75 0 0.75])
                plot(Time,TorqueActive(1,:),'linewidth',2,'linestyle','-','color','r')
                hold off
                grid on
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                Llg=legend('u_r_1','u_m_1','u_b_1','u_a_1');
                set(Llg,'orientation','horizontal')
                xlim([Time(1) Time(end)])
                %
                subplot(4,1,2)
                plot(Time,TorqueDesire(2,:),'linewidth',2,'linestyle','-','color','b')
                hold on
                plot(Time,TorqueMonoOptimal(2,:),'linewidth',2,'linestyle','-.','color','g')
                plot(Time,TorqueBicepsOptimal(1,:),'linewidth',2,'linestyle','--','Color',[0.75 0 0.75])
                plot(Time,TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','--','Color',[0.68 0.45 0])
                plot(Time,TorqueActive(2,:),'linewidth',2,'linestyle','-','color','r')
                hold off
                grid on
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                Llg=legend('u_r_2','u_m_2','u_b_1','u_b_2','u_a_2');
                set(Llg,'orientation','horizontal')
                xlim([Time(1) Time(end)])
                %
                subplot(4,1,3)
                plot(Time,TorqueDesire(3,:),'linewidth',2,'linestyle','-','color','b')
                hold on
                plot(Time,TorqueMonoOptimal(3,:),'linewidth',2,'linestyle','-.','color','g')
                plot(Time,TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','--','Color',[0.75 0 0.75])
                plot(Time,TorqueBicepsOptimal(3,:),'linewidth',2,'linestyle','--','Color',[0.68 0.45 0])
                plot(Time,TorqueActive(3,:),'linewidth',2,'linestyle','-','color','r')
                hold off
                grid on
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                Llg=legend('u_r_3','u_m_3','u_b_2','u_b_3','u_a_3');
                set(Llg,'orientation','horizontal')
                xlim([Time(1) Time(end)])
                %
                subplot(4,1,4)
                plot(Time,TorqueDesire(4,:),'linewidth',2,'linestyle','-','color','b')
                hold on
                plot(Time,TorqueMonoOptimal(4,:),'linewidth',2,'linestyle','-.','color','g')
                plot(Time,TorqueBicepsOptimal(3,:),'linewidth',2,'linestyle','--','Color',[0.68 0.45 0])
                plot(Time,TorqueActive(4,:),'linewidth',2,'linestyle','-','color','r')
                hold off
                grid on
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                ylabel('u_4 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                Llg=legend('u_r_4','u_m_4','u_b_3','u_a_4');
                set(Llg,'orientation','horizontal')
                xlim([Time(1) Time(end)])
                xlabel('time (s)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');

    
            figure('name',[' Desired Torque*\dot{q} vs time : ',Name])
                subplot(4,1,1,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                plot(Time,TorqueDesire(1,:).*D1Q1,'linewidth',2,'displayname','power_r_1')
                hold on
                plot(Time,TorqueActive(1,:).*D1Q1,'linewidth',2,'color','r','displayname','power_a_1')
                hold off
                legend(gca,'show','Orientation','horizontal')
                if(strcmp(Period,'2Cycle'))
                    hold on
                    plot(Time+Time(end),TorqueDesire(1,:).*D1Q1,'linewidth',2,'linestyle','-.')
                    plot(Time+Time(end),TorqueActive(1,:).*D1Q1,'linewidth',2,'color','r','linestyle','-.')
                    hold off
                end
                ylabel('${u_1 * \dot q_1}$', 'interpreter','latex','fontsize',12,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                DeltaTorque=max(TorqueDesire(1,:).*D1Q1)-min(TorqueDesire(1,:).*D1Q1);
                YLIM1=([min(TorqueDesire(1,:).*D1Q1)-0.1*DeltaTorque,max(TorqueDesire(1,:).*D1Q1)+0.1*DeltaTorque]);
                DeltaTorque=max(TorqueActive(1,:).*D1Q1)-min(TorqueActive(1,:).*D1Q1);
                YLIM2=([min(TorqueActive(1,:).*D1Q1)-0.1*DeltaTorque,max(TorqueActive(1,:).*D1Q1)+0.1*DeltaTorque]);
                ylim([min([YLIM1(1),YLIM2(1)])   max([YLIM1(2),YLIM2(2)])  ])
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                %
                subplot(4,1,2,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                plot(Time,TorqueDesire(2,:).*D1Q2,'linewidth',2,'displayname','power_r_2')
                hold on
                plot(Time,TorqueActive(2,:).*D1Q2,'linewidth',2,'color','r','displayname','power_a_2')
                hold off
                legend(gca,'show','Orientation','horizontal')
                if(strcmp(Period,'2Cycle'))
                    hold on
                    plot(Time+Time(end),TorqueDesire(2,:).*D1Q2,'linewidth',2,'linestyle','-.')
                    plot(Time+Time(end),TorqueActive(2,:).*D1Q2,'linewidth',2,'color','r','linestyle','-.')
                    hold off
                end
                ylabel('${u_2 * \dot q_2}$', 'interpreter','latex','fontsize',12,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                DeltaTorque=max(TorqueDesire(2,:).*D1Q2)-min(TorqueDesire(2,:).*D1Q2);
                YLIM1=([min(TorqueDesire(2,:).*D1Q2)-0.1*DeltaTorque,max(TorqueDesire(2,:).*D1Q2)+0.1*DeltaTorque]);
                DeltaTorque=max(TorqueActive(2,:).*D1Q2)-min(TorqueActive(2,:).*D1Q2);
                YLIM2=([min(TorqueActive(2,:).*D1Q2)-0.1*DeltaTorque,max(TorqueActive(2,:).*D1Q2)+0.1*DeltaTorque]);
                ylim([min([YLIM1(1),YLIM2(1)])   max([YLIM1(2),YLIM2(2)])  ])
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                %
                subplot(4,1,3,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                plot(Time,TorqueDesire(3,:).*D1Q3,'linewidth',2,'displayname','power_r_1')
                hold on
                plot(Time,TorqueActive(3,:).*D1Q3,'linewidth',2,'color','r','displayname','power_a_3')
                hold off
                legend(gca,'show','Orientation','horizontal')
                if(strcmp(Period,'2Cycle'))
                    hold on
                    plot(Time+Time(end),TorqueDesire(3,:).*D1Q3,'linewidth',2,'linestyle','-.')
                    plot(Time+Time(end),TorqueActive(3,:).*D1Q3,'linewidth',2,'color','r','linestyle','-.')
                    hold off
                end
                xlabel('Time','fontsize',12,'FontName','mwa_cmb10');
                ylabel('${u_3 * \dot q_3}$', 'interpreter','latex','fontsize',12,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                DeltaTorque=max(TorqueDesire(3,:).*D1Q3)-min(TorqueDesire(3,:).*D1Q3);
                YLIM1=([min(TorqueDesire(3,:).*D1Q3)-0.1*DeltaTorque,max(TorqueDesire(3,:).*D1Q3)+0.1*DeltaTorque]);
                DeltaTorque=max(TorqueActive(3,:).*D1Q3)-min(TorqueActive(3,:).*D1Q3);
                YLIM2=([min(TorqueActive(3,:).*D1Q3)-0.1*DeltaTorque,max(TorqueActive(3,:).*D1Q3)+0.1*DeltaTorque]);
                DeltaTorque=max(TorqueDesire(3,:).*D1Q3)-min(TorqueDesire(3,:).*D1Q3);
                YLIM1=([min(TorqueDesire(3,:).*D1Q3)-0.1*DeltaTorque,max(TorqueDesire(3,:).*D1Q3)+0.1*DeltaTorque]);
                DeltaTorque=max(TorqueActive(3,:).*D1Q3)-min(TorqueActive(3,:).*D1Q3);
                YLIM2=([min(TorqueActive(3,:).*D1Q3)-0.1*DeltaTorque,max(TorqueActive(3,:).*D1Q3)+0.1*DeltaTorque]);
                ylim([min([YLIM1(1),YLIM2(1)])   max([YLIM1(2),YLIM2(2)])  ])
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                %
                subplot(4,1,4,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                plot(Time,TorqueDesire(4,:).*D1Q4,'linewidth',2,'displayname','power_r_1')
                hold on
                plot(Time,TorqueActive(4,:).*D1Q4,'linewidth',2,'color','r','displayname','power_a_3')
                hold off
                legend(gca,'show','Orientation','horizontal')
                if(strcmp(Period,'2Cycle'))
                    hold on
                    plot(Time+Time(end),TorqueDesire(4,:).*D1Q4,'linewidth',2,'linestyle','-.')
                    plot(Time+Time(end),TorqueActive(4,:).*D1Q4,'linewidth',2,'color','r','linestyle','-.')
                    hold off
                end
                xlabel('Time','fontsize',12,'FontName','mwa_cmb10');
                ylabel('${u_4 * \dot q_4}$', 'interpreter','latex','fontsize',12,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                DeltaTorque=max(TorqueDesire(4,:).*D1Q4)-min(TorqueDesire(4,:).*D1Q4);
                YLIM1=([min(TorqueDesire(4,:).*D1Q4)-0.1*DeltaTorque,max(TorqueDesire(4,:).*D1Q4)+0.1*DeltaTorque]);
                DeltaTorque=max(TorqueActive(4,:).*D1Q4)-min(TorqueActive(4,:).*D1Q4);
                YLIM2=([min(TorqueActive(4,:).*D1Q4)-0.1*DeltaTorque,max(TorqueActive(4,:).*D1Q4)+0.1*DeltaTorque]);
                DeltaTorque=max(TorqueDesire(4,:).*D1Q4)-min(TorqueDesire(4,:).*D1Q4);
                YLIM1=([min(TorqueDesire(4,:).*D1Q4)-0.1*DeltaTorque,max(TorqueDesire(4,:).*D1Q4)+0.1*DeltaTorque]);
                DeltaTorque=max(TorqueActive(4,:).*D1Q4)-min(TorqueActive(4,:).*D1Q4);
                YLIM2=([min(TorqueActive(4,:).*D1Q4)-0.1*DeltaTorque,max(TorqueActive(4,:).*D1Q4)+0.1*DeltaTorque]);
                ylim([min([YLIM1(1),YLIM2(1)])   max([YLIM1(2),YLIM2(2)])  ])
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
           
        figure('name',['Torque vs Angle : ',Name])
                ap=get(gca,'position');
                subplot(2+(rB>0),4,1,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                plot(rad2deg( Q1),TorqueDesire(1,:),'linewidth',2)
                %title('Optimal Required and Compliance Torque-Angle Profile','FontSize',16);
                xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                hold on
                plot(rad2deg(Q1),TorqueActive(1,:),'linewidth',2,'color','r','linestyle','-.')
                hold off
                legend('u_r_1','u_a_1','Orientation','horizontal')
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                %
                subplot(2+(rB>0),4,2,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                plot(rad2deg( Q2),TorqueDesire(2,:),'linewidth',2)
                xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                hold on
                plot(rad2deg(Q2),TorqueActive(2,:),'linewidth',2,'color','r','linestyle','-.')
                hold off
                legend('u_r_2','u_a_2','Orientation','horizontal')
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                %
                subplot(2+(rB>0),4,3,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                plot(rad2deg(Q3),TorqueDesire(3,:),'linewidth',2)
                xlabel('q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                hold on
                plot(rad2deg(Q3),TorqueActive(3,:),'linewidth',2,'color','r','linestyle','-.')
                hold off
                legend('u_r_3','u_a_3','Orientation','horizontal')
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                %
                subplot(2+(rB>0),4,4,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                plot(rad2deg(Q4),TorqueDesire(4,:),'linewidth',2)
                xlabel('q_4 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                ylabel('u_4 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                hold on
                plot(rad2deg(Q4),TorqueActive(4,:),'linewidth',2,'color','r','linestyle','-.')
                hold off
                legend('u_r_4','u_a_4','Orientation','horizontal')
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                
                subplot(2+(rB>0),4,5,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                plot(rad2deg(Q1),TorqueMonoOptimal(1,:),'linewidth',3,'color','g')
                xlabel('q_1 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                ylabel('u_1 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                hold on
                plot(rad2deg( Q1),TorqueDesire(1,:)-TorqueBicepsOptimal(1,:),'linewidth',2,'linestyle','-.')
                hold off
                legend('u_m_1','u_r_1','Orientation','horizontal')
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                %
                subplot(2+(rB>0),4,6,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                plot(rad2deg(Q2),TorqueMonoOptimal(2,:),'linewidth',3,'color','g')
                xlabel('q_2 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                ylabel('u_2 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                hold on
                plot(rad2deg( Q2),TorqueDesire(2,:)-TorqueBicepsOptimal(1,:)-TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-.')
                hold off
                legend('u_m_2','u_r_2-u_b_2','Orientation','horizontal')
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                %
                subplot(2+(rB>0),4,7,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                plot(rad2deg(Q3),TorqueMonoOptimal(3,:),'linewidth',3,'color','g')
                xlabel('q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                ylabel('u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                hold on
                plot(rad2deg(Q3),TorqueDesire(3,:)-TorqueBicepsOptimal(2,:)-TorqueBicepsOptimal(3,:),'linewidth',2,'linestyle','-.')
                hold off
                legend('u_m_3','u_r_3-u_b_2-u_b_3','Orientation','horizontal')
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                %
                subplot(2+(rB>0),4,8,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                plot(rad2deg(Q4),TorqueMonoOptimal(4,:),'linewidth',3,'color','g')
                xlabel('q_4 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                ylabel('u_4 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                grid on
                set(gca,'YMinorGrid','on')
                hold on
                plot(rad2deg(Q4),TorqueDesire(4,:)-TorqueBicepsOptimal(3,:),'linewidth',2,'linestyle','-.')
                hold off
                legend('u_m_4','u_r_4-u_b_3','Orientation','horizontal')
                set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                

                if(rB>0)
                    
                    sh1=subplot(3,4,10);
                    sp1=get(sh1,'position');
                    plot(rad2deg(Q2+Q3),2*TorqueBicepsOptimal(2,:),'linewidth',3,'Color',[0.75 0 0.75])
                    set(sh1,'position',[sp1(1)+.1,sp1(2:end)]); 
                    xlabel('q_2+q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                    ylabel('u_2+u_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                    grid on
                    set(gca,'YMinorGrid','on')
                    hold on
                    plot(rad2deg(Q2+Q3),TorqueDesire(2,:)+TorqueDesire(3,:)-TorqueMonoOptimal(2,:)-TorqueMonoOptimal(3,:)-TorqueBicepsOptimal(3,:),'linewidth',2,'linestyle','-.')
                    hold off
                    legend('2\times u_b_2','u_r_2+u_r_3-u_m_2-u_m_3-u_b_3','Orientation','horizontal')
                    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                    
                    sh2=subplot(3,4,12);
                    sp2=get(sh2,'position');
                    set(sh2,'position',[sp2(1)-.1,sp2(2:end)]); 
                    plot(rad2deg(Q3+Q4),2*TorqueBicepsOptimal(3,:),'linewidth',3,'Color',[0.68 0.45 0])
                    xlabel('q_3+q_4 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                    ylabel('u_3+u_4 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
                    grid on
                    set(gca,'YMinorGrid','on')
                    hold on
                    plot(rad2deg(Q3+Q4),TorqueDesire(3,:)+TorqueDesire(4,:)-TorqueMonoOptimal(3,:)-TorqueMonoOptimal(4,:)-TorqueBicepsOptimal(2,:),'linewidth',2,'linestyle','-.')
                    hold off
                    legend('2\times u_b_3','u_r_3+u_r_4-u_m_3-u_m_4-u_b_2','Orientation','horizontal')
                    set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
                

                end
                
                
%     figure('name',['Torque vs Angle : ',Name])
%                 ap=get(gca,'position');
%                 set(gca,'FontWeight','bold','FontSize',12,'FontName','mwa_cmb10');
%                 plot(rad2deg( Q1+Q2+Q3),...
%                     TorqueActive(1,:)+TorqueActive(2,:)+TorqueActive(3,:),...
%                     'linewidth',2)
%                 %title('Optimal Required and Compliance Torque-Angle Profile','FontSize',16);
%                 xlabel('q_1+q_2+q_3 (deg)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
%                 ylabel('u_a_1+u_a_2+u_a_3 (N.m)','FontWeight','bold','FontSize',14,'FontName','mwa_cmb10');
%                 grid on
%                 set(gca,'YMinorGrid','on')
%                 hold on
%                 
                        
    end
            

end



end