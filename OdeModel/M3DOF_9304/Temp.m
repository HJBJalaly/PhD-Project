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


%%
% home
x=2:0.001:3;
y1=sqrt(-x+5)+4;
y2=y1;%sqrt(-x+5)+4;
X=[x, x(end:-1:1)];
Y=[y1,y2(end:-1:1)];
DY=diff(Y)./diff(X);

Joint=2;

ThetaStep=( (max(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)))  - min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2))))/200);
ThetaS=min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2))) :ThetaStep:max(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)));
ThetaShift=ThetaS-min(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)));
ThetaShiftScale = ThetaShift* (deg2rad(270) /  max(ThetaShift));
ThetaStepscale  = ThetaStep* (deg2rad(270) /  max(ThetaShift));

tau=interp1(Q_Opt(Joint, 1: floor(size(Q_Opt,2)/2)),Torque_Opt(Joint, 1: floor(size(Q_Opt,2)/2)),ThetaS);
TauMin=abs(min(tau));
tauShift=tau+TauMin;
tauShiftScale=(tauShift)/max(tauShift);

Range=5:25; % range of fiting
SubRange1=[1:Range(1)]; % for first :  range of outlayer
% SubRange1=[Range(end)+1:length(tau)]; % for end :  range of outlayer

SubRange2=[Range(1)+1:length(tau)]; % for first : unchanged value
% SubRange2=[1:Range(end)]; % for end : unchanged value

Curve=fit( ThetaShiftScale(Range)', tauShiftScale(Range)', 'poly2');
NewTau=feval(Curve,ThetaShiftScale(SubRange1))';
tauShiftScaleVar= [NewTau tauShiftScale(SubRange2)]; % for first


DTa=differential(tauShiftScale,ThetaShiftScale,(ThetaStepscale));

DTaVar=differential(tauShiftScaleVar,ThetaShiftScale,(ThetaStepscale));
    
Sum(1)=tauShiftScale(1);
% SumVar(1)=tauShiftScaleVar(1);
for i=2:length(DTa)
    Sum(i)=Sum(i-1)+(DTa(i))^4*ThetaStepscale;
%     SumVar(i)=SumVar(i-1)+(DTaVar(i))*ThetaStepscale*0;
end


subplot(3,1,1)
plot(Q_Opt(Joint,:),Torque_Opt(Joint,:),'linewidth',2)
grid on
xlabel('$q$','interpreter','latex','fontsize',14)
ylabel('$\tau$','interpreter','latex','fontsize',28)
subplot(3,1,2)
plot(ThetaShiftScale,tauShiftScale,'linewidth',2)
grid on
xlabel('$scaled\; and\; shifted\; q$','interpreter','latex','fontsize',14)
ylabel('$scaled\; and\; shifted\;\tau$','interpreter','latex','fontsize',14)
% ylim([-.5 1.5])
hold all
plot(ThetaShiftScale,Sum,'linewidth',2)
% plot(ThetaShiftScale,SumVar,'linewidth',2)
hold off

subplot(3,1,3)
plot(ThetaShiftScale,DTa,'linewidth',2)
grid on
hold on
plot([0 5],[1 1],'color','r','linestyle','-.','linewidth',2)
plot([0 5],[-1 -1],'color','r','linestyle','-.','linewidth',2)
xlabel('$q$','interpreter','latex','fontsize',14)
ylabel('${\frac{{\partial\tau}}{\partial{q}}}$ ','interpreter','latex','fontsize',38)
hold off
ylim([-4,1.5])