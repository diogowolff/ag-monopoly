function [val] = SMM(guess)
    typeDistributionMean = ...
    [1*10^(-5), 1330, guess(1), guess(2)]; % Original A was 1.9*10^-3
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
    
    meanS = sqrt(25000^2 - 5100^2);
    vec = [expected_value, variance];
    true = [4340, meanS];
    
    val = (vec-true)*(vec-true)';
    
end

