%% Calibration
clc

Upper=[15 15 15 15 15 7.5 7.5]'*1;
% Upper=[10 10 10 10 10 10 10]'*1;
Lower=-Upper;
DataLinkagePos=@Data5LinkagePosBNoisy_Var_10;
PopulationSize=2000;

[x,fval,exitflag,output,population,score] =...
    GaCalibration(Lower,Upper,PopulationSize,DataLinkagePos);


%% Estimation Forward Kinematics

clear all

% Calibrated Parameters
x=[4.4879    2.1247    1.1852    8.7946   -1.5202    1.0727   -6.8509];
L0=300+x(1);
L1=400+x(2);
L2=400+x(3);
L3=550+x(4);
L4=550+x(5);

Data5LinkageQReal;
Data5LinkagePosReal;
[Xr,Yr]=DataGeneration(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
Q1(114)=[]; Q2(114)=[];   Xm(114)=[];Ym(114)=[];    Xr(114)=[];Yr(114)=[]; 

%
figure
plot3(Q1,Q2,Xr,'linestyle','none','marker','*')
grid on
set(gca,'fontsize',12)
set(gca,'fontsize',12,'fontweight','bold')
xlabel('q_1 (deg)','fontsize',14,'fontweight','bold')
ylabel('q_2 (deg)','fontsize',14,'fontweight','bold')
zlabel('x (mm)','fontsize',14,'fontweight','bold')

figure
plot3(Q1,Q2,Yr,'linestyle','none','marker','*')
grid on
set(gca,'fontsize',12)
set(gca,'fontsize',12,'fontweight','bold')
xlabel('q_1 (deg)','fontsize',14,'fontweight','bold')
ylabel('q_2 (deg)','fontsize',14,'fontweight','bold')
zlabel('y (mm)','fontsize',14,'fontweight','bold')




for i=1:5
    for j=1:5
        [fx,gofx]=CreateFit(Q1+x(6),Q2+x(7),Xr,['Poly',num2str(i),num2str(j)],'X','ShowOff');
        XfitRes{i,j}=fx;
        XfitError{i,j}=gofx;
        XfitRMS(i,j)=gofx.rmse;
        
        [fy,gofy]=CreateFit(Q1+x(6),Q2+x(7),Yr,['Poly',num2str(i),num2str(j)],'Y','ShowOff');
        YfitRes{i,j}=fy;
        YfitError{i,j}=gofy;
        YfitRMS(i,j)=gofy.rmse;
        
        PosfitMaxError(i,j) =max(sqrt(  (fx(Q1+x(6),Q2+x(7))-Xr).^2+ (fy(Q1+x(6),Q2+x(7))-Yr).^2 ));
        PosfitmeanError(i,j)=mean(sqrt( (fx(Q1+x(6),Q2+x(7))-Xr).^2+ (fy(Q1+x(6),Q2+x(7))-Yr).^2 ));
    end
end

mesh(XfitRMS)

figure
plot(XfitRes{4,4},[Q1+x(6),Q2+x(7)],Xr)
grid on
set(gca,'fontsize',12)
set(gca,'fontsize',12,'fontweight','bold')
xlabel('q_1 (deg)','fontsize',14,'fontweight','bold')
ylabel('q_2 (deg)','fontsize',14,'fontweight','bold')
zlabel('x (mm)','fontsize',14,'fontweight','bold')

figure
plot(YfitRes{4,4},[Q1+x(7),Q2+x(7)],Yr)
grid on
set(gca,'fontsize',12)
set(gca,'fontsize',12,'fontweight','bold')
xlabel('q_1 (deg)','fontsize',14,'fontweight','bold')
ylabel('q_2 (deg)','fontsize',14,'fontweight','bold')
zlabel('y (mm)','fontsize',14,'fontweight','bold')

%% Estimation Inverse Kinematics

clear all

% Calibrated Parameters
x=[4.4879    2.1247    1.1852    8.7946   -1.5202    1.0727   -6.8509];
L0=300+x(1);
L1=400+x(2);
L2=400+x(3);
L3=550+x(4);
L4=550+x(5);

Data5LinkageQReal;
Data5LinkagePosReal;
[Xr,Yr]=DataGeneration(x(1),x(2),x(3),x(4),x(5),x(6),x(7));
Q1(114)=[]; Q2(114)=[];   Xm(114)=[];Ym(114)=[];    Xr(114)=[];Yr(114)=[]; 
Q1=Q1+x(6);
Q2=Q2+x(7);


Q2In=atan2d(Yr,Xr)-acosd((L2^2+Xr.^2+Yr.^2-L4^2)./(2*sqrt(Xr.^2+Yr.^2)*L2));
Q2-Q2In;

figure
plot3(Xr,Yr,Q1,'linestyle','none','marker','*')
grid on
set(gca,'fontsize',12)
set(gca,'fontsize',12,'fontweight','bold')
xlabel('x (mm)','fontsize',14,'fontweight','bold')
ylabel('y (mm)','fontsize',14,'fontweight','bold')
zlabel('q_1 (deg)','fontsize',14,'fontweight','bold')

figure
plot3(Xr,Yr,Q2,'linestyle','none','marker','*')
grid on
set(gca,'fontsize',12)
set(gca,'fontsize',12,'fontweight','bold')
xlabel('x (mm)','fontsize',14,'fontweight','bold')
ylabel('y (mm)','fontsize',14,'fontweight','bold')
zlabel('q_2 (deg)','fontsize',14,'fontweight','bold')




for i=1:5
    for j=1:5
        [fq1,gofq1]=CreateFit(Xr+L0,Yr,Q1,['Poly',num2str(i),num2str(j)],'q_1','ShowOff');
        Q1fitRes{i,j}=fq1;
        Q1fitError{i,j}=gofq1;
        Q1fitRMS(i,j)=gofq1.rmse;
        
        [fq2,gofq2]=CreateFit(Xr,Yr,Q2,['Poly',num2str(i),num2str(j)],'q_2','ShowOff');
        Q2fitRes{i,j}=fq2;
        Q2fitError{i,j}=gofq2;
        Q2fitRMS(i,j)=gofq2.rmse;
        
        Q1fitMaxError(i,j) =max(abs(fq1(Xr+L0,Yr)-Q1 ));
        Q2fitMaxError(i,j) =max(abs(fq2(Xr,Yr)-Q2 ));
    end
end

figure
plot(Q1fitRes{4,4},[Xr+L0 Yr],Q1)
grid on
set(gca,'fontsize',12)
set(gca,'fontsize',12,'fontweight','bold')
xlabel('x(mm)','fontsize',14,'fontweight','bold')
ylabel('y (mm)','fontsize',14,'fontweight','bold')
zlabel('q_1 (deg)','fontsize',14,'fontweight','bold')

figure
plot(Q2fitRes{4,4},[Xr Yr],Q2)
grid on
set(gca,'fontsize',12)
set(gca,'fontsize',12,'fontweight','bold')
xlabel('x(mm)','fontsize',14,'fontweight','bold')
ylabel('y (mm)','fontsize',14,'fontweight','bold')
zlabel('q_2 (deg)','fontsize',14,'fontweight','bold')