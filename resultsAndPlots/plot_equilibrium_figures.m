clear;
close all;
addpath('../classes');
addpath('plotFunctions');
addpath('../plotFunctions/export_fig_updated');
addpath('../plotFunctions');
rng(1);

%% Individual models
for modelNameString = { ...
        'interval',                  ...
        'interval_high_mh_variance'};
    
    CalculationData = load(modelNameString{1});
    plotEquilibriumAndOptimum( ...
        CalculationData.Model, CalculationData.pEquilibrium, CalculationData.DEquilibrium, ...
        CalculationData.pEfficient, CalculationData.DEfficient, modelNameString{1});
end;