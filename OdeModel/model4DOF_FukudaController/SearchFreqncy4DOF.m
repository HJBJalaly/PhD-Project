function SearchFreqncy4DOF()

L = 0.25;
d0= 0.05;
SpaceRange=[0.06:0.02:0.58];
length(SpaceRange)
for i=1:length(SpaceRange)
    i
    WW(i)=OdeSirFukudaControllerFindNaFr4DOF(SpaceRange(i),0);
end

save('SearchFreqncy4DOFData.mat');
end