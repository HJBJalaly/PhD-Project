function [value,isterminal,direction]=EventTouchDown(t,x,L_fem, L_tib)

% Locate the time when height passes through zero in a
% decreasing direction and stop integration.
q1=x(1);
q2=x(2);
q3=x(3);
q4=x(4);
q5=x(5);

Ytib=- L_fem*cos(q1 + q5) + L_fem*cos(q2 + q5) - L_tib*cos(q1 + q3 + q5) + L_tib*cos(q2 + q4 + q5);
Xtib=+ L_fem*sin(q1 + q5) - L_fem*sin(q2 + q5) + L_tib*sin(q1 + q3 + q5) - L_tib*sin(q2 + q4 + q5);
 
value = [Ytib];     % Detect height = 0
isterminal= [ Xtib>0 ];   % Stop the integration
direction = [-1 ];   % ( - --> + )
% ( + --> - )
end

