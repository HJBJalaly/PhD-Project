

%% 2DoF 

g=9.81;


m=1;
L=1;

K1=50;
K2=30;


time=[0 1];
OdeOpt= odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,5));
InitState=deg2rad([ -29.44 112.02 ,0 0 , 0]);

[T,Y] = ode15s(@(t,Y)SirDyn2DoFPassive(t,Y,g,L,m,K1,K2), time,InitState,OdeOpt);
q1=Y(:,1)';
q2=Y(:,2)';
EnergyABS=Y(end,end)

RPos=L*[cos(q1)+cos(q1+q2);
        sin(q1)+sin(q1+q2)];

plot(RPos(1,1), RPos(2,1),'linestyle','none','marker','o','markersize',6,'markeredgecolor','r','markerfacecolor','r')
hold all
plot(RPos(1,:),RPos(2,:),'linewidth',2,'linestyle','-','color','b')
xlabel('x')
ylabel('y')
grid on
ezplot(['x^2+y^2-',num2str((2*L)^2)])
hold off
% legend('Desired','Static path planing','dynamic path planing')
    
mL1=m;
mL2=m;

LL1=L;
LL2=L;

for i=1:size(Y,1)
    q1=Y(i,1);
    q2=Y(i,2);
    
    D1q1=Y(i,3);
    D1q2=Y(i,4);
    
    D1q=[D1q1 D1q2]';


    % Dynamic
    MM=[mL2*LL1*LL2*cos(q2)+mL2*LL1^2+mL2*LL2^2/3+mL1*LL1^2/3, mL2*LL2*(3*LL1*cos(q2)+2*LL2)/6;
        mL2*LL2*(3*LL1*cos(q2)+2*LL2)/6, mL2*LL2^2/3];


    CC= [-mL2*LL2*sin(q2)*D1q2*LL1/2, -mL2*LL2*sin(q2)*(D1q1+D1q2)*LL1/2; 
          mL2*LL2*LL1*D1q1*sin(q2)/2, 0];


    GG= [g*(mL2*LL2*cos(q1+q2)+LL1*cos(q1)*mL1+2*LL1*cos(q1)*mL2)/2 + K1*q1;
         g*mL2*LL2*cos(q1+q2)/2 + K2*q2];

    %  Control Law
    Gamma=[0;0];
    % Apply Torqe
    D2q(i,:)= (MM^-1*(-CC*D1q-GG+Gamma))';
    Test(:,i)=GG;
    
end
