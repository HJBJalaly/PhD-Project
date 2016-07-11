function Dq=differential(q,Time,tstep)
% q is vector or matrix of trajectories, which each row is trajectory of a
% joint and each colomn is value of trajectories per step.


Dq= zeros(size(q));

for i=1:size(q,1)
    curve = fit( Time(1:4)', q(i,1:4)', 'poly2');
    P=coeffvalues(curve);
    Dq(i,1)=2*P(1)*Time(1)+P(2);
end

Dq(:,2)=(q(:,3)-q(:,1))/( diff(Time([1,3])) );              % 3 point differentiation

for i=1+2:size(q,2)-2
    
    Dq(:,i)=(-q(:,i+2)+8*q(:,i+1)-8*q(:,i-1)+q(:,i-2))/ (12*tstep);
                                                % 5 point differentiation
    
end
Dq(:,end-1)=(q(:,end)-q(:,end-2))/(diff(Time([end-2,end])));      % 3 point differentiation
% Dq(:,end-1)=(q(:,end-1)-q(:,end-2))/(diff(Time([end-2,end-1])));      % 3 point differentiation


for i=1:size(q,1)
    curve = fit( Time(end-3:end)', q(i,end-3:end)', 'poly2');
    P=coeffvalues(curve);
    Dq(i,end)=2*P(1)*Time(end)+P(2);
end
1;