clear;
addpath('../classes');
rng(1);

% Input model parameters

typeDistributionMean = ...
    [1*10^(-5), 1330, 4340, 24400]; % Original A was 1.9*10^-3
typeDistributionLogCovariance = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
     -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
     -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
      0     0     0    0.25]; % ???
  
costOfPublicFunds = 0;

% Calculation parameters
populationSize = 5*10^4;

CalculationParametersEquilibrium.behavioralAgents = 0.01;
CalculationParametersEquilibrium.fudge            = 1e-6;
CalculationParametersEquilibrium.maxIterations    = 1e4;
CalculationParametersEquilibrium.tolerance        = 1;

CalculationParametersOptimum.maxIterations        = 1e3;
CalculationParametersOptimum.tolerance            = 0.01;
CalculationParametersOptimum.knitro               = 'true';
CalculationParametersOptimum.knitroMultistartN    = 300;

% List of models
modelName{1}              = 'interval_censnorm_test';
slopeVector{1}            = 0:0.04:1;
moralHazardLogVariance{1} = 0.28;

% Loop
nSimulations = length(modelName);

for i = 1 : nSimulations
    innerTypeDistributionLogCovariance = typeDistributionLogCovariance;

    innerTypeDistributionLogCovariance(2, 2) = moralHazardLogVariance{i};
    Model = healthcaralognormalmodel_censnorm(slopeVector{i}, typeDistributionMean, innerTypeDistributionLogCovariance, -5000000);
    
    Population = population(Model, populationSize);

        Types = zeros(populationSize, 2);

    for j=1:populationSize
       Types(j,:) = [Population.typeList{j}.M, Population.typeList{j}.S];
    end

    simu_mean = mean(Types(:, 1) - normpdf((-5000000-Types(:, 1))./Types(:, 2))./...
        (1-normcdf((-5000000-Types(:, 1))./Types(:, 1))).*Types(:, 1));

    simu_var = mean(Types(:, 2).^2.*(1 + (-5000000-Types(:, 1))./Types(:, 2).*...
        normpdf((-5000000-Types(:, 1))./Types(:, 2))./(1-normcdf((-5000000-Types(:, 1))./Types(:, 2))) - ...
        (normpdf((-5000000-Types(:, 1))./Types(:, 2))./(1-normcdf((-5000000-Types(:, 1))./Types(:, 2)))).^2));

    sqrtsim = sqrt(simu_var);
    
    i_name = modelName{i};


end