clear;
addpath('../classes');
rng(1);

% Input model parameters
meanS = sqrt(25000^2 - 5100^2);

F = @(x) [4340 - x(1) - normpdf(-x(1)/x(2))/(1-normcdf(-x(1)/x(2)))*x(2);
    meanS - x(2)^2*(1 - x(1)/x(2)*normpdf(-x(1)/x(2))/(1-normcdf(-x(1)/x(2))) - ...
    (normpdf(-x(1)/x(2))/(1-normcdf(-x(1)/x(2))))^2)];

x0 = [2000; 20000];

[x,fval] = fsolve(F,x0);

typeDistributionMean = ...
    [1*10^(-5), 1330, x(1), x(2)]; % Original A was 1.9*10^-3
typeDistributionLogCovariance = ...
    [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
     -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
     -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
      0     0     0    0.25]; % ???
  
costOfPublicFunds = 0;

% Calculation parameters
populationSize = 5*10^3;
slopeVector{1}            = 0:0.01:1;

innerTypeDistributionLogCovariance = typeDistributionLogCovariance;

Model = healthcaralognormalmodel_censnorm(slopeVector{1}, typeDistributionMean, innerTypeDistributionLogCovariance);

Population = population(Model, populationSize);

Types = zeros(populationSize, 2);

for i=1:populationSize
   Types(i,:) = [Population.typeList{i}.M, Population.typeList{i}.S];
end

expected_value = mean(Types(:, 1) - normpdf(-Types(:, 1)./Types(:, 2))./...
    (1-normcdf(-Types(:, 1)./Types(:, 1))).*Types(:, 1));

variance = mean(Types(:, 2).^2.*(1 - Types(:, 1)./Types(:, 2).*...
    normpdf(-Types(:, 1)./Types(:, 2))./(1-normcdf(-Types(:, 1)./Types(:, 2))) - ...
    (normpdf(-Types(:, 1)./Types(:, 2))./(1-normcdf(-Types(:, 1)./Types(:, 2)))).^2));