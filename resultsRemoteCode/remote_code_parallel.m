clear;
addpath('classes');
addpath('resultsAndPlots');
rng(1);

% Input model parameters
meanS = sqrt(25000^2 - 5100^2);
typeDistributionMean = [1*10^(-5), 1330, 4340, meanS]; % Original A was 1.9*10^-3
typeDistributionLogCovariance = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
     -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
     -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
      0     0     0    0.25]; % ???
  
costOfPublicFunds = 0;

% Calculation parameters
populationSize = 10^5;

CalculationParametersEquilibrium.behavioralAgents = 0.01;
CalculationParametersEquilibrium.fudge            = 1e-6;
CalculationParametersEquilibrium.maxIterations    = 1e4;
CalculationParametersEquilibrium.tolerance        = 1;

CalculationParametersOptimum.maxIterations        = 1e3;
CalculationParametersOptimum.tolerance            = 0.01;
CalculationParametersOptimum.knitro               = 'true';
CalculationParametersOptimum.knitroMultistartN    = 300;

% List of models
modelName{1}              = 'interval';
slopeVector{1}            = 0:0.04:1;
moralHazardLogVariance{1} = 0.28;

modelName{2}              = 'interval_high_mh_variance';
slopeVector{2}            = 0:0.04:1;
moralHazardLogVariance{2} = 0.98;

% Loop
nSimulations = length(modelName);

poolObject = parpool(10);

for i = 1 : nSimulations
    innerTypeDistributionLogCovariance = typeDistributionLogCovariance;

    innerTypeDistributionLogCovariance(2, 2) = moralHazardLogVariance{i};
    Model(i) = healthcaralognormalmodel(slopeVector{i}, typeDistributionMean, innerTypeDistributionLogCovariance);
    
    Population(i) = population(Model(i), populationSize, 10);
end

delete(poolObject);
poolObject = parpool(2);


parfor i = 1 : nSimulations
    
    [pEquilibrium, DEquilibrium, ACEquilibrium, ComputationOutputEquilibrium] = ...
            Population(i).findequilibrium(CalculationParametersEquilibrium);
    WEquilibrium = Population(i).welfare(pEquilibrium, ...
                                          costOfPublicFunds);
    
    [pEfficient, PiEfficient, ComputationOutputEfficient] = ...
            findefficient(Population(i), CalculationParametersOptimum);
    DEfficient = Population(i).demand(pEfficient);
    
    i_name = modelName{i};
    Model_i = Model(i);
    
    display(i_name);
    display(ComputationOutputEfficient);
    display(ComputationOutputEquilibrium);    
    parsave(i_name, Model_i, pEquilibrium, DEquilibrium, ...
                ACEquilibrium, ComputationOutputEquilibrium, ...
            WEquilibrium, pEfficient, PiEfficient, ...
            ComputationOutputEfficient, DEfficient, ...
            CalculationParametersEquilibrium, CalculationParametersOptimum, ...
            populationSize)


end

delete(poolObject);