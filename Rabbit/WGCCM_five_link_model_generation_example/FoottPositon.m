function [p_tib2 ] = FoottPositon(Q,L_fem, L_tib)

q1=Q(1);
q2=Q(2);
q3=Q(3);
q4=Q(4);
q5=Q(5);


p_tib2 = [ + L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5),
           - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + L_tib*cos(q2 + q4 + q5)]; 


end

