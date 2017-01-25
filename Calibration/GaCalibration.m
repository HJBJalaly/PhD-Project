function [x,fval,exitflag,output,population,score] = GaCalibration(lb,ub,PopulationSize_Data,DataLinkagePos)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = gaoptimset;
%% Modify options setting
options = gaoptimset(options,'PopInitRange', [lb,ub]');
options = gaoptimset(options,'PopulationSize', PopulationSize_Data);
options = gaoptimset(options,'Display', 'off');
options = gaoptimset(options,'PlotFcns', {  @gaplotbestf @gaplotbestindiv });
options = gaoptimset(options,'CrossoverFcn', {  @crossoverheuristic 1.2 });
[x,fval,exitflag,output,population,score] = ...
ga(@(x)CostError(x,DataLinkagePos),7,[],[],[],[],lb,ub,[],[],options);
