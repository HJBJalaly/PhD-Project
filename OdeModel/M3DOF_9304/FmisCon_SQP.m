function [x,fval,exitflag,output,lambda,grad,hessian] = FmisCon_SQP(CostFun,NonCons,x0,MaxFunEvals_Data,MaxIter_Data,TolFun_Data,TolX_Data,TolCon_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimoptions('fmincon');
%% Modify options setting
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'MaxFunEvals', MaxFunEvals_Data);
options = optimoptions(options,'MaxIter', MaxIter_Data);
options = optimoptions(options,'TolFun', TolFun_Data);
options = optimoptions(options,'TolX', TolX_Data);
options = optimoptions(options,'PlotFcns', {  @optimplotfval @optimplotconstrviolation });
options = optimoptions(options,'Algorithm', 'sqp');
options = optimoptions(options,'TolCon', TolCon_Data);
options = optimoptions(options,'UseParallel', 'always');
[x,fval,exitflag,output,lambda,grad,hessian] = ...
fmincon(CostFun,x0,[],[],[],[],[],[],NonCons,options);
