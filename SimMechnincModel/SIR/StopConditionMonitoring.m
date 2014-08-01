function STOP=StopConditionMonitoring(input)

pos2X=input(1);
pos2Y=input(2);
vel2Y=input(5);
pos1X=input(7);
STOP=0;

global ReachLine


if(pos2Y>-0.01 && vel2Y>0 && pos2X>pos1X )
    ReachLine=1;
end

% if(ReachLine &&  pos2Y < -0.01 <) 
%     1;
% end

if(ReachLine &&  pos2Y <= -0.01) 
    STOP=1;
end