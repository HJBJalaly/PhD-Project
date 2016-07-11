function [value,isterminal,direction] = StopCond(t,theta)
% Locate the time when height passes through zero in a
% decreasing direction and stop integration.
value = [2.9*pi/2-theta(1)];     % Detect height = 0
isterminal= [+1 ];   % Stop the integration
direction = [-1 ]; 
% ( + --> - )
end