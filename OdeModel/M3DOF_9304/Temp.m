% figure('name','Torque vs Angle۲')
        subplot(3,1,1)
        plot(Q_Opt(1,1:49),Torque_Opt(1,1:49),'r','linewidth',2)
        title('Torque Angle Profile','FontWeight','bold')
        hold on
        plot(Q_Opt(1,49:100),Torque_Opt(1,49:100),'b','linewidth',2)
        
        plot(Q_Opt(1,1),Torque_Opt(1,1),'linewidth',2,'linestyle','none','marker','*','markersize',6)
        xlabel('q_1 (rad)','fontsize',12)
        ylabel('\tau_1','fontsize',14)
        hold off
        grid on

        subplot(3,1,2)
        plot(Q_Opt(2,1:49),Torque_Opt(2,1:49),'r','linewidth',2)
        hold on
        plot(Q_Opt(2,49:100),Torque_Opt(2,49:100),'b','linewidth',2)
        plot(Q_Opt(2,1),Torque_Opt(2,1),'linewidth',2,'linestyle','none','marker','*','markersize',6)
        xlabel('q_2 (rad)','fontsize',12)
        ylabel('\tau_2','fontsize',14)
        hold off
        grid on

        subplot(3,1,3)
        plot(Q_Opt(3,1:49),Torque_Opt(3,1:49),'r','linewidth',2)
        hold on
        plot(Q_Opt(3,49:100),Torque_Opt(3,49:100),'b','linewidth',2)
        
        plot(Q_Opt(3,1),Torque_Opt(3,1),'linewidth',2,'linestyle','none','marker','*','markersize',6)
        xlabel('q_3 (rad)','fontsize',12)
        ylabel('\tau_3','fontsize',14)
        hold off
        grid on

        
%%



ThetaStep=deg2rad( 1);
ThetaS=deg2rad(0:1:270);


tau=0.21*(ThetaS-0.75*pi).*(ThetaS-0.25*pi).*(ThetaS-1.25*pi)+2.5;

k0=1000;
R0=1;
q00=1;  % for cubic

nvars=3;
lb=[0 0 0];
PopInitRange=[10 5e-2 5e-2; 1000 1 1];
PopulationSize=200;
InitialPopulation=[k0*rand(PopulationSize,1) R0*rand(PopulationSize,1) q00*rand(PopulationSize,1)];
CostParam=@(Param)FindBestParamCost(Param,ThetaStep,ThetaS,tau);

[ParamA,fval,exitflag,output,population,score] = ...
    Ga_FindParamOfNonLinearSpring(CostParam,nvars,lb,PopInitRange,PopulationSize,InitialPopulation);
disp(output.message)

k=ParamA(1);
R=ParamA(2);
q0=ParamA(3);

NonLinearSpring(ThetaStep,ThetaS,tau,k,R,q0)


