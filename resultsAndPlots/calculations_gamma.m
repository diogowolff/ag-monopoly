clear;
addpath('../classes');
rng(1);

% Input model parameters

typeDistributionMean = ...
    [1*10^(-5), 1330, 0.03, 150000]; % Original A was 1.9*10^-3
typeDistributionLogCovariance = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
     -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
     -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
      0     0     0    0.25]; % ???
  
costOfPublicFunds = 0;

% Calculation parameters
populationSize = 5*10^6;

CalculationParametersEquilibrium.behavioralAgents = 0.01;
CalculationParametersEquilibrium.fudge            = 1e-6;
CalculationParametersEquilibrium.maxIterations    = 1e4;
CalculationParametersEquilibrium.tolerance        = 1;

CalculationParametersOptimum.maxIterations        = 1e3;
CalculationParametersOptimum.tolerance            = 0.01;
CalculationParametersOptimum.knitro               = 'true';
CalculationParametersOptimum.knitroMultistartN    = 300;

% List of models
modelName{1}              = 'interval_gamma_v9';
slopeVector{1}            = 0:.04:1;
moralHazardLogVariance{1} = 0.28;

% Loop
nSimulations = length(modelName);
i=1;
for i = 1 : nSimulations
    innerTypeDistributionLogCovariance = typeDistributionLogCovariance;

    innerTypeDistributionLogCovariance(2, 2) = moralHazardLogVariance{i};
    Model = healthcaralognormalmodel_gamma(slopeVector{i}, typeDistributionMean, innerTypeDistributionLogCovariance);
    
    Population = population(Model, populationSize);
    populationSize = length(Population.typeList);
    [pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
            Population.findequilibrium(CalculationParametersEquilibrium);
    WEquilibrium = Population.welfare(pEquilibrium, ...
                                          costOfPublicFunds);
    
    [pEfficient, PiEfficient, ComputationOutputEfficient] = ...
            findefficient(Population, CalculationParametersOptimum);
    DEfficient = Population.demand(pEfficient);

    [pWelfare, WWelfare, ComputationOutputWelfare] = ...
            findwelfaremax(Population, costOfPublicFunds, CalculationParametersOptimum);
    DWelfare = Population.demand(pWelfare);
    
    i_name = modelName{i};
    
    display(i_name);
    display(ComputationOutputEfficient);
    display(ComputationOutputEquilibrium);   
    display(ComputationOutputWelfare);
    parsave(i_name, Model, pEquilibrium, DEquilibrium, ...
                ACEquilibrium, ComputationOutputEquilibrium, ...
            WEquilibrium, pEfficient, PiEfficient, ...
            ComputationOutputEfficient, DEfficient, ...
            pWelfare, WWelfare, ComputationOutputWelfare, ...
            DWelfare, ...
            CalculationParametersEquilibrium, CalculationParametersOptimum, ...
            populationSize)


end