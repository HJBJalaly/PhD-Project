disp('  ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('waiting ....')

clear all
load t9conserv.mat
load phase9conserv.mat

load t9Damper.mat
load phase9Damper.mat

gamma=.000000001;
omega=9;
delta=1e-6;
SampleRate=500;
rp=50;
scale=2*pi;

[BetaOptimal,CostAct,CostD2Q,Up]=...
                    OptimalRotaryCompliance(phase9,t9,rp,gamma,omega,delta,SampleRate,scale);               
Ua=t9-Up;
%% 
disp('  ')
figure
plot(phase9,t9,'linewidth',2)
hold all
plot(phase9,Up,'linewidth',2)
plot(phase9,Ua,'r','linewidth',2)
grid on
xlabel('\phi (rad)','fontsize',12,'FontWeight','bold')
ylabel('\tau (N.m)','fontsize',12,'FontWeight','bold')
Le1=legend('u_r','u_p','u_a');
set(Le1,'FontSize',12,'Orientation','horizontal')
xlim([0, phase9(end)])

SumCompliance=sum(Up)*2*pi/length(phase9);

SumReq=sum(t9)*2*pi/length(phase9);
WorkReq=sum(abs(t9))*2*pi/length(phase9);
CostReq=1/2*t9'*t9;

SumAct=sum(Ua)*2*pi/length(phase9);
WorkAct=sum(abs(Ua))*2*pi/length(phase9);

Title=sprintf('%22s  % 11s %10s % 9s % 10s'  ,   'IntU2','C.S.(D2)','Total','Work','Area');
Result_Req=sprintf('%-11s %11g %11.2e %11.2e % 10.2e  % 10.2e  ',   'Required: ',CostReq, 0 , CostReq   , WorkReq , SumReq);
Result_Act=sprintf('%-11s %11g %11.2e %11.2e % 10.2e  % 10.2e\n',  'Actuation:',CostAct,CostD2Q,CostAct+gamma*CostD2Q,WorkAct,SumAct);
disp(Title)
disp(Result_Req)
disp(Result_Act)
disp(['Conservative Condition: ',num2str(SumCompliance,'%3.2e')])
%% Bernsteine

clear all
load t9conserv.mat
load phase9conserv.mat

Bv=@(x,n,v) nchoosek(n,v)*(x^v)*((1-x)^(n-v));

clf
plot(phase9,t9,'linewidth',2)
hold on
for nn=1000;
    nn
    x=[];
    Up=[];
    
    for i=0:nn
        fvn(i+1)=interp1(phase9,t9,i/nn*phase9(end));
    end
    for j=1:500:length(phase9)
        x(end+1)=phase9(j);
        Up(end+1)=0;
        for i=0:nn
            Up(end)=Up(end)+fvn(i+1)*Bv(phase9(j)/2/pi,nn,i);
        end
    end
    plot(x,Up,'r','linewidth',2)
    drawnow
end
grid on 
Le1=legend('main','estimation');
set(Le1,'FontSize',12,'Orientation','horizontal')
xlim([0, phase9(end)])
xlabel('\phi (rad)','fontsize',12,'FontWeight','bold')
ylabel('\tau (N.m)','fontsize',12,'FontWeight','bold')
%% Sine --> Fourier Series

clear all
% load t9conserv.mat
% load phase9conserv.mat

load t9Damper.mat
load phase9Damper.mat

figure
plot(phase9,t9,'linewidth',2)
hold on
for nn=5:5:20
    a=[];
    b=[];
    DeltaPhi=diff(phase9(1:2));
    for i=1:nn
        a(i)=sum(t9.*cos(i*phase9))*DeltaPhi/pi;
        b(i)=sum(t9.*sin(i*phase9))*DeltaPhi/pi;
    end
        
    x=phase9(1:end);
    Up=zeros(size(x));
    for i=1:nn
        Up=Up+a(i)*cos(i*x)+b(i)*sin(i*x);
    end
    plot(x,Up,'r','linewidth',2)
    grid on
    drawnow
end


Ua=t9-Up;

disp('  ')
figure
plot(phase9,t9,'linewidth',2)
hold all
plot(phase9,Up,'linewidth',2)
plot(phase9,Ua,'r','linewidth',2)
grid on
xlabel('\phi (rad)','fontsize',12,'FontWeight','bold')
ylabel('\tau (N.m)','fontsize',12,'FontWeight','bold')
Le1=legend('u_r','u_p','u_a');
set(Le1,'FontSize',12,'Orientation','horizontal')
xlim([0, phase9(end)])

SumCompliance=sum(Up)*2*pi/length(phase9);

SumReq=sum(t9)*2*pi/length(phase9);
WorkReq=sum(abs(t9))*2*pi/length(phase9);
CostReq=1/2*t9'*t9;

SumAct=sum(Ua)*2*pi/length(phase9);
WorkAct=sum(abs(Ua))*2*pi/length(phase9);
CostAct=1/2*Ua'*Ua;
Title=sprintf('%22s  % 11s %10s % 9s % 10s'  ,   'IntU2','Work','Area');
Result_Req=sprintf('%-11s %11g  % 10.2e  % 10.2e  ',   'Required: ',CostReq  , WorkReq , SumReq);
Result_Act=sprintf('%-11s %11g  % 10.2e  % 10.2e\n',  'Actuation:',CostAct,WorkAct,SumAct);
disp(Title)
disp(Result_Req)
disp(Result_Act)
disp(['Conservative Condition: ',num2str(SumCompliance,'%3.2e')])
%% Sine --> LS based
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('waiting ....')

clear all
% load t9conserv.mat
% load phase9conserv.mat

load t9Damper.mat
load phase9Damper.mat

SampleRate=20;
nn=20;
gamma=2e-7;
omega=9;

figure
plot(phase9,t9,'linewidth',2)
hold all

PHI=[];
PHIc=[];
for i=1:SampleRate:length(phase9)
    PHI(:,end+1) = [cos((1:nn)*phase9(i)) sin((1:nn)*phase9(i))]';
end
for i=1:length(phase9)
    PHIc(:,end+1) = [cos((1:nn)*phase9(i)) sin((1:nn)*phase9(i))]';
end
R=diag([(nn:-1:1).^2 (nn:-1:1).^2]);

Beta=(PHI*PHI'+omega*gamma*R*(PHI*PHI')*R)\(PHI*t9(1:SampleRate:length(phase9)));

Up=PHIc'*Beta;
Ua=t9-Up;

SumCompliance=sum(Up)*2*pi/length(Up);

SumReq=sum(t9)*2*pi/length(t9);
WorkReq=sum(abs(t9))*2*pi/length(t9);
CostReq=1/2*t9'*t9;

SumAct=sum(Ua)*2*pi/length(phase9);
WorkAct=sum(abs(Ua))*2*pi/length(phase9);
CostAct=1/2*Ua'*Ua;
CostD2Q=(PHI'*R*Beta)'*(PHI'*R*Beta);

Title=sprintf('%22s  % 11s %10s % 9s % 10s'  ,   'IntU2','C.S.(D2)','Total','Work','Area');
Result_Req=sprintf('%-11s %11g %11.2e %11.2e % 10.2e  % 10.2e  ',   'Required: ',CostReq, 0 , CostReq   , WorkReq , SumReq);
Result_Act=sprintf('%-11s %11g %11.2e %11.2e % 10.2e  % 10.2e\n',  'Actuation:',CostAct,CostD2Q,CostAct+gamma*CostD2Q,WorkAct,SumAct);
disp(Title)
disp(Result_Req)
disp(Result_Act)
disp(['Conservative Condition: ',num2str(SumCompliance,'%3.2e')])


plot(phase9,Up,'linewidth',2)
plot(phase9,Ua,'linewidth',2)
grid on
grid on
xlabel('\phi (rad)','fontsize',12,'FontWeight','bold')
ylabel('\tau (N.m)','fontsize',12,'FontWeight','bold')
Le1=legend('u_r','u_p','u_a');
set(Le1,'FontSize',12,'Orientation','horizontal')
xlim([0, phase9(end)])
