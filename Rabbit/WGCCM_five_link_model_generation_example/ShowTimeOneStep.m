function  ShowTimeOneStep(Tc,SolC,MotionData,c,Theta_plus,Theta_minus,Alfa,L_fem, L_tib, L_torso,mu,Ma)


q1=SolC(:,1)' ;q2=SolC(:,2)' ;q3=SolC(:,3)' ;q4=SolC(:,4)' ;q5=SolC(:,5)';
dq1=SolC(:,6)';dq2=SolC(:,7)';dq3=SolC(:,8)';dq4=SolC(:,9)';dq5=SolC(:,10)';

p_tib2 = [ + L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5)
           - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + L_tib*cos(q2 + q4 + q5)];

       
v_hip =[dq1.*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) - dq2.*(L_fem*cos(q2 + q5) + L_tib*cos(q2 + q4 + q5)) + dq5.*(L_fem*cos(q1 + q5) - L_fem*cos(q2 + q5) + L_tib*cos(q1 + q3 + q5) - L_tib*cos(q2 + q4 + q5)) + L_tib*dq3.*cos(q1 + q3 + q5) - L_tib*dq4.*cos(q2 + q4 + q5)
        dq1.*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) - dq2.*(L_fem*sin(q2 + q5) + L_tib*sin(q2 + q4 + q5)) + dq5.*(L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5)) + L_tib*dq3.*sin(q1 + q3 + q5) - L_tib*dq4.*sin(q2 + q4 + q5)];

figure(6)
    plot(Tc,MotionData(1:4,:)')
    legend('u_1','u_2','u_3','u_4')
    grid on
    xlabel('Time (s)')
    ylabel('Torque (N.m)')

    

figure(7)
    Thetas=c*[q1;q2;q3;q4;q5];
    Ss=(Thetas-Thetas(1))/(Thetas(end)-Thetas(1));
    Hds=zeros(4,1);
    DHds=zeros(4,1);
    Hdr=zeros(4,1);
    DHdr=zeros(4,1);
    Thetar=linspace(Theta_plus,Theta_minus,length(Tc));
    Sr=(Thetar-Theta_plus)/(Theta_minus-Theta_plus);
    for i=1:length(Thetar)
        [Hds(:,end+1),DHds(:,end+1)]=BezierFunction(Ma,Thetas(i),Alfa,Theta_plus,Theta_minus);
        [Hdr(:,end+1),DHdr(:,end+1)]=BezierFunction(Ma,Thetar(i),Alfa,Theta_plus,Theta_minus);
    end
    Hds(:,1)=[];
    DHds(:,1)=[];
    Hdr(:,1)=[];
    DHdr(:,1)=[];
    subplot(3,2,1)
    plot(Ss,q1)
    hold all
    plot(Ss,Hds(1,:),'-.')
    plot(Sr,Hdr(1,:),'-.r')
    legend1=legend('sim','ref','des');
    set(legend1,'Orientation','horizontal');
    ylabel('q_1')
    hold off
    subplot(3,2,2)
    plot(Ss,q2)
    hold all
    plot(Ss,Hds(2,:),'-.')
    plot(Sr,Hdr(2,:),'-.r')
    legend1=legend('sim','ref','des');
    set(legend1,'Orientation','horizontal');
    hold off
    ylabel('q_2')
    subplot(3,2,3)
    plot(Ss,q3)
    hold all
    plot(Ss,Hds(3,:),'-.')
    plot(Sr,Hdr(3,:),'-.r')
    legend1=legend('sim','ref','des');
    set(legend1,'Orientation','horizontal');
    hold off
    ylabel('q_3')
    subplot(3,2,4)
    plot(Ss,q4)
    hold all
    plot(Ss,Hds(4,:),'-.')
    plot(Sr,Hdr(4,:),'-.r')
    legend1=legend('sim','ref','des');
    set(legend1,'Orientation','horizontal');
    ylabel('q_4')
    hold off
    subplot(3,2,5)
    plot(Ss,q5)
    ylabel('q_5')
    subplot(3,2,6)
    plot(Ss,Thetas)
    ylabel('\theta')
    
figure(8)
    subplot(3,2,1)
    plot(Ss,dq1)
    set(legend1,'Orientation','horizontal');
    ylabel('Dq_1')
    
    subplot(3,2,2)
    plot(Ss,dq2)
    ylabel('Dq_2')
    
    subplot(3,2,3)
    plot(Ss,dq3)
    ylabel('dq_3')
    
    subplot(3,2,4)
    plot(Ss,dq4)
    ylabel('dq_4')
    
    subplot(3,2,5)
    plot(Ss,dq5)
    ylabel('dq_5')
    
    

figure(9)
    plot(Tc,MotionData(5,:))%Ft
    hold all
    plot(Tc,MotionData(6,:))%Fn
    plot(Tc,MotionData(6,:)*mu,'--')%mu*Fn
    plot(Tc,-MotionData(6,:)*mu,'--')%mu*Fn
    grid on
    legend('F_t','F_n','\mu*Fn')
    hold off
    xlabel('Time (s)')
    ylabel('GRF (N)')

    
figure(10)
    subplot(2,1,1)
    plot(Tc,p_tib2(1,:))
    hold all
    plot(Tc,p_tib2(2,:))
    grid on
    legend('p^x_t_2','p^y_t_2')
    hold off
    xlabel('Time (s)')
    ylabel('p_t_i_b_2 (m)')
    subplot(2,1,2)
    plot(Tc,v_hip(1,:))
    hold all
    plot(Tc,v_hip(2,:))
    grid on
    legend('v^x_h','v^y_h')
    hold off
    xlabel('Time (s)')
    ylabel('v_h_i_p (m/s)')
