function [KE,PE ] = KineticPotentialEnergy(Q,DQ,L_fem, L_tib, L_torso, Lc_fem, Lc_tib, Lc_torso, M_fem, M_tib, M_torso, XX_fem, XX_tib, XX_torso, g )

q1=Q(1);
q2=Q(2);
q3=Q(3);
q4=Q(4);
q5=Q(5);
y1=0;

dq1=DQ(1);
dq2=DQ(2);
dq3=DQ(3);
dq4=DQ(4);
dq5=DQ(5);
dy1=0;
dx1=0;

KE=(XX_fem*(dq1 + dq5)^2)/2 + (XX_fem*(dq2 + dq5)^2)/2 + (M_tib*((dy1 + L_tib*dq1*sin(q1 + q3 + q5) + L_tib*dq3*sin(q1 + q3 + q5) + L_tib*dq5*sin(q1 + q3 + q5) - Lc_tib*dq2*sin(q2 + q4 + q5) - Lc_tib*dq4*sin(q2 + q4 + q5) - Lc_tib*dq5*sin(q2 + q4 + q5) + L_fem*dq1*sin(q1 + q5) - L_fem*dq2*sin(q2 + q5) + L_fem*dq5*sin(q1 + q5) - L_fem*dq5*sin(q2 + q5))^2 + (dx1 + L_tib*dq1*cos(q1 + q3 + q5) + L_tib*dq3*cos(q1 + q3 + q5) + L_tib*dq5*cos(q1 + q3 + q5) - Lc_tib*dq2*cos(q2 + q4 + q5) - Lc_tib*dq4*cos(q2 + q4 + q5) - Lc_tib*dq5*cos(q2 + q4 + q5) + L_fem*dq1*cos(q1 + q5) - L_fem*dq2*cos(q2 + q5) + L_fem*dq5*cos(q1 + q5) - L_fem*dq5*cos(q2 + q5))^2))/2 + (XX_torso*dq5^2)/2 + (XX_tib*(dq1 + dq3 + dq5)^2)/2 + (XX_tib*(dq2 + dq4 + dq5)^2)/2 + (M_tib*((dx1 + dq1*cos(q1 + q3 + q5)*(L_tib - Lc_tib) + dq3*cos(q1 + q3 + q5)*(L_tib - Lc_tib) + dq5*cos(q1 + q3 + q5)*(L_tib - Lc_tib))^2 + (dy1 + dq1*sin(q1 + q3 + q5)*(L_tib - Lc_tib) + dq3*sin(q1 + q3 + q5)*(L_tib - Lc_tib) + dq5*sin(q1 + q3 + q5)*(L_tib - Lc_tib))^2))/2 + (M_torso*((dx1 + dq1*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) + dq5*(L_fem*cos(q1 + q5) - cos(q5)*(L_torso - Lc_torso) + L_tib*cos(q1 + q3 + q5)) + L_tib*dq3*cos(q1 + q3 + q5))^2 + (dy1 + dq1*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) + dq5*(L_fem*sin(q1 + q5) - sin(q5)*(L_torso - Lc_torso) + L_tib*sin(q1 + q3 + q5)) + L_tib*dq3*sin(q1 + q3 + q5))^2))/2 + (M_fem*((dy1 + dq1*(L_fem*sin(q1 + q5) + L_tib*sin(q1 + q3 + q5)) + dq5*(L_fem*sin(q1 + q5) - Lc_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5)) + L_tib*dq3*sin(q1 + q3 + q5) - Lc_fem*dq2*sin(q2 + q5))^2 + (dx1 + dq1*(L_fem*cos(q1 + q5) + L_tib*cos(q1 + q3 + q5)) + dq5*(L_fem*cos(q1 + q5) - Lc_fem*cos(q2 + q5) + L_tib*cos(q1 + q3 + q5)) + L_tib*dq3*cos(q1 + q3 + q5) - Lc_fem*dq2*cos(q2 + q5))^2))/2 + (M_fem*((dy1 + dq1*(sin(q1 + q5)*(L_fem - Lc_fem) + L_tib*sin(q1 + q3 + q5)) + dq5*(sin(q1 + q5)*(L_fem - Lc_fem) + L_tib*sin(q1 + q3 + q5)) + L_tib*dq3*sin(q1 + q3 + q5))^2 + (dx1 + dq1*(cos(q1 + q5)*(L_fem - Lc_fem) + L_tib*cos(q1 + q3 + q5)) + dq5*(cos(q1 + q5)*(L_fem - Lc_fem) + L_tib*cos(q1 + q3 + q5)) + L_tib*dq3*cos(q1 + q3 + q5))^2))/2;


PE =g*(M_tib*(2*y1 - L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + Lc_tib*cos(q2 + q4 + q5)) - M_fem*(cos(q1 + q5)*(L_fem - Lc_fem) - 2*y1 + L_tib*cos(q1 + q3 + q5)) + M_fem*(2*y1 - L_fem*cos(q1 + q5) + Lc_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5)) + M_torso*(2*y1 - L_fem*cos(q1 + q5) + cos(q5)*(L_torso - Lc_torso) - L_tib*cos(q1 + q3 + q5)) + M_tib*(2*y1 - cos(q1 + q3 + q5)*(L_tib - Lc_tib)));
 


end

