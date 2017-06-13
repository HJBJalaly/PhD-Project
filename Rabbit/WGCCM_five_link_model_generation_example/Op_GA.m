function [x,fval,exitflag,output,population,score] =...
    Op_GA(CostFun,ConstraintFun,nvars,InitialPopulation,PopInitRange,PopulationSize,GenerationsLimit,StallGenLimit)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = gaoptimset;
%% Modify options setting
options = gaoptimset(options,'PopInitRange', PopInitRange);
options = gaoptimset(options,'PopulationSize', PopulationSize);
options = gaoptimset(options,'Generations', GenerationsLimit);
options = gaoptimset(options,'StallGenLimit', StallGenLimit);
options = gaoptimset(options,'InitialPopulation', InitialPopulation);
options = gaoptimset(options,'Display', 'off');
options = gaoptimset(options,'PlotFcns', { @gaplotbestf});
% options = optimoptions(options,'PlotFcns', {  @optimplotfval @optimplotconstrviolation });
options = gaoptimset(options,'Vectorized', 'off');
options = gaoptimset(options,'UseParallel', 'always');
options = gaoptimset(options,'Display', 'iter');
options = gaoptimset(options,'TolFun',1e-3);
options = gaoptimset(options,'TolCon',1e-3);
[x,fval,exitflag,output,population,score] = ...
    ga(CostFun,nvars,[],[],[],[],[],[],ConstraintFun,[],options);
