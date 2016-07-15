function Q=InRangeShifter(Q)

while(mean(Q)>pi)
    Q=Q-2*pi;
end

while(mean(Q)<-pi)
    Q=Q+2*pi;
end