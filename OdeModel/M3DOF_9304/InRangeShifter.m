function Q=InRangeShifter(Q)

while(mean(Q)>pi)
    Q=Q-pi;
end

while(mean(Q)<-pi)
    Q=Q+pi;
end