function SearchFreqncy2DOF()

L=0.3;
SpaceRange=[0.06:0.02:0.58];
length(SpaceRange)
for i=1:length(SpaceRange)
    i
    WW(i)=OdeSirFukudaControllerFindNaFr2DOF(SpaceRange(i),0);
end

save('SearchFreqncy2DOFData.mat');
end