function [x,fval,exitflag,output,population,score] = ...
    GA_FminCost(CostParam,nvars,lb,PopInitRange_Data,PopulationSize_Data,InitialPopulation_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = gaoptimset;
%% Modify options setting
options = gaoptimset(options,'PopInitRange', PopInitRange_Data);
options = gaoptimset(options,'PopulationSize', PopulationSize_Data);
options = gaoptimset(options,'InitialPopulation', InitialPopulation_Data);
options = gaoptimset(options,'CreationFcn', @gacreationlinearfeasible);
options = gaoptimset(options,'MutationFcn', {  @mutationuniform 0.25 });
options = gaoptimset(options,'Display', 'off');
options = gaoptimset(options,'PlotFcns', {  @gaplotbestf @gaplotbestindiv });
[x,fval,exitflag,output,population,score] = ...
            ga(CostParam,nvars,[],[],[],[],lb,[],[],[],options);
