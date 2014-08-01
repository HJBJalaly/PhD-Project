function OdeDoubleMassSpringSimple()

figure
home

m1=1;
m2=1;
B1=0.5;
B2=0.5;


OdeOpt= odeset('RelTol',1e-5,'AbsTol',1e-5*ones(1,4));
InitState=[ 0.1 0.1  0.05 0.05];
u1=0;
u2=0;
[T,X] = ode15s(@(t,X)SirDyn(t,X,m1,m2,B1,B2,u1,u2),[0 100],InitState,OdeOpt);

subplot(3,1,1)
plot(T,X(:,1)*100,'linewidth',2)
hold all
plot(T,X(:,2)*100,'linewidth',2)
hold off
grid on
legend({'$$r_1$$','$$r_2$$'},'Interpreter','latex',...
        'EdgeColor',[1 1 1],'Orientation','horizontal','Location','Best')
ylabel('state varibales (cm)')

subplot(3,1,2)
plot(T,X(:,3)*100,'linewidth',2)
hold all
plot(T,X(:,4)*100,'linewidth',2)
hold off
grid on
legend({'$$\dot{r}_1$$','$$\dot{r}_2$$'},'Interpreter','latex',...
        'EdgeColor',[1 1 1],'Orientation','horizontal','Location','Best')
ylabel('state varibales (cm/s)')

subplot(3,1,3)
plot(T,X(:,1)*100,'linewidth',2)
hold all
plot(T,X(:,2)*100+X(:,1)*100,'linewidth',2)
hold off
grid on
legend({'$$x_{m_1}$$','$$x_{m_2}$$'},'Interpreter','latex',...
        'EdgeColor',[1 1 1],'Orientation','horizontal','Location','Best')
xlabel('T (s)')
ylabel('mass position (cm)')


1;

end

function Dx=SirDyn(t,X,m1,m2,B1,B2,u1,u2)

Dx=zeros(4,1);
M=[m1+m2 , m2 ; 0,  m2];
C=[B1 0; 0 B2];

% x=Tr
% Dx= Ax+Bu
A = [zeros(2)   eye(2)  ; zeros(2)  -M^-1 * C];
B=[zeros(2);eye(2)];

Dx=A*X+B*[u1 ; u2];

end

