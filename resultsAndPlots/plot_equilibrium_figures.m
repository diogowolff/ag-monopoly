clear;
close all;
addpath('../classes');
addpath('plotFunctions');
addpath('../plotFunctions/export_fig_updated');
addpath('../plotFunctions');
rng(1);

%% Individual models
for modelNameString = { ...
        'interval_censnorm_server'};
    
    CalculationData = load(modelNameString{1});
    plotEquilibriumAndOptimum( ...
        CalculationData.Model, CalculationData.pEquilibrium, CalculationData.DEquilibrium, ...
        CalculationData.pEfficient, CalculationData.DEfficient, modelNameString{1});
end;