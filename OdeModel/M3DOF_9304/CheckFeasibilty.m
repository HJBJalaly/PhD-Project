function CheckFeasibilty()
%%
close all
clear all
load testCc_Lin_1_WithOptimizedTrajectory3_ForPaper_WindowsOnly_Uni_rQ_20
% load Task_eli_withJointLimitation(ForPaper_JustOnWindows)_rQ=10.mat
% load testCc_cir_14good.mat
%%
QJ=q(:,Middle1:Middle2);
DQJ=Dq(:,Middle1:Middle2);
E=[];

for i=1:size(QJ,2)
    q1f=QJ(1,i);
    q2f=QJ(2,i);
    q3f=QJ(3,i);
    

    MM=[mL3*LL1*LL3*cos(q2f+q3f)+LL1*LL2*cos(q2f)*mL2+2*LL1*LL2*cos(q2f)*mL3+mL3*LL2*LL3*cos(q3f)+mL3*LL2^2+mL3*LL3^2/3+mL3*LL1^2+LL1^2*mL2+mL1*LL1^2/3+mL2*LL2^2/3,mL3*LL1*LL3*cos(q2f+q3f)/2+LL1*LL2*cos(q2f)*mL2/2+LL1*LL2*cos(q2f)*mL3+mL3*LL2*LL3*cos(q3f)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL3*(3*LL1*cos(q2f+q3f)+2*LL3+3*LL2*cos(q3f))/6;
        mL3*LL1*LL3*cos(q2f+q3f)/2+LL1*LL2*cos(q2f)*mL2/2+LL1*LL2*cos(q2f)*mL3+mL3*LL2*LL3*cos(q3f)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,mL3*LL2*LL3*cos(q3f)+mL3*LL2^2+mL3*LL3^2/3+mL2*LL2^2/3,(2*LL3+3*LL2*cos(q3f))*LL3*mL3/6;
        mL3*LL3*(3*LL1*cos(q2f+q3f)+2*LL3+3*LL2*cos(q3f))/6,(2*LL3+3*LL2*cos(q3f))*LL3*mL3/6,mL3*LL3^2/3];

    E(i)=DQJ(:,i)'*MM*DQJ(:,i);
end
%%
figure(1)
subplot(2,1,1)
plot(QJ(1,:),E)
    hold all
    plot(QJ(2,:),E)
    plot(QJ(3,:),E)
    legend('E(q_1)','E(q_2)','E(q_3)','orientation','horizontal')
    grid on
subplot(2,1,2)
    plot(QJ(1,:),E)
    hold all
    plot(QJ(2,:)+QJ(1,:),E)
    plot(QJ(3,:)+QJ(1,:),E)
    legend('E(q_1)','E(q_1+q_2)','E(q_1+q_3)','orientation','horizontal')
    grid on
figure(2)
    subplot(2,2,1)
    plot3(QJ(1,:),QJ(2,:),E)
    xlabel('q_1')
    ylabel('q_2')
    grid on
    box on
    subplot(2,2,2)
    plot3(QJ(1,:),QJ(3,:),E)
    xlabel('q_1')
    ylabel('q_3')
    grid on
    box on
    subplot(2,2,3)
    plot3(QJ(3,:),QJ(2,:),E)
    xlabel('q_3')
    ylabel('q_2')
    grid on
    box on
%%
    f1=diff(E)./diff(QJ(1,:));
    f2=diff(E)./diff(QJ(2,:));
    figure(3)
    subplot(2,2,1)
    plot(QJ(1,1:end-1),f1)
    subplot(2,2,2)
    plot3(QJ(1,1:end-1),QJ(2,1:end-1),f1)
    subplot(2,2,2)
    plot3(QJ(1,1:end-1),QJ(2,1:end-1),f1)
    grid on
    subplot(2,2,3)
    plot(QJ(2,1:end-1),f2)
    subplot(2,2,4)
    plot3(QJ(1,1:end-1),QJ(2,1:end-1),f2)
    grid on
    