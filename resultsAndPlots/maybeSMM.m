clear;
addpath('../Classes');
%chuta um par, calcula, pega as vars, faz o match, etc etc


%very wide range to find parameters

[teste2_a, teste2_b] = meshgrid(100:100:5000, 1000:1000:50000);

teste2_mean = zeros(50,50);
teste2_var = zeros(50,50);

for k=1:50
    for l=1:50
         typeDistributionMean = ...
        [1*10^(-5), 1330, teste2_a(1,k), teste2_b(l,1)]; % Original A was 1.9*10^-3
        typeDistributionLogCovariance = ...
            [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
             -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
             -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
              0     0     0    0.25]; % ???

        costOfPublicFunds = 0;

        % Calculation parameters
        populationSize = 5*10^1;
        slopeVector{1}            = 0:0.01:1;

        innerTypeDistributionLogCovariance = typeDistributionLogCovariance;

        Model = healthcaralognormalmodel_censnorm(slopeVector{1}, typeDistributionMean, innerTypeDistributionLogCovariance);

        Population = population(Model, populationSize);

        Types = zeros(populationSize, 2);

        for i=1:populationSize
           Types(i,:) = [Population.typeList{i}.M, Population.typeList{i}.S];
        end

        teste2_mean(k,l) = mean(Types(:, 1) - normpdf(-Types(:, 1)./Types(:, 2))./...
            (1-normcdf(-Types(:, 1)./Types(:, 1))).*Types(:, 1));

        teste2_var(k,l) = mean(Types(:, 2).^2.*(1 - Types(:, 1)./Types(:, 2).*...
            normpdf(-Types(:, 1)./Types(:, 2))./(1-normcdf(-Types(:, 1)./Types(:, 2))) - ...
            (normpdf(-Types(:, 1)./Types(:, 2))./(1-normcdf(-Types(:, 1)./Types(:, 2)))).^2));
    end
end

teste2_vardif = (teste2_var - meanS^2);
teste2_meandif = (teste2_mean - 4340);

teste2_sumsq = teste2_vardif.^2 + teste2_meandif.^2;
[row, col] = find(teste2_sumsq==min(teste2_sumsq(:)));

sqrt(teste2_var(row, col))
teste2_mean(row, col)

teste2_a(37,5)
teste2_b(37,5)

% the variance moment dominates, so i identify the sigma and vary the mu
% to find the best one, results were 8k and 37k approximately
%below, i show that the results are somewhat smooth (depends a lot on
%sample size)

typeDistributionMean = ...
    [1*10^(-5), 1330, 8000, 37000]; % Original A was 1.9*10^-3
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

simu_mean = mean(Types(:, 1) - normpdf(-Types(:, 1)./Types(:, 2))./...
    (1-normcdf(-Types(:, 1)./Types(:, 1))).*Types(:, 1));

simu_var = mean(Types(:, 2).^2.*(1 - Types(:, 1)./Types(:, 2).*...
    normpdf(-Types(:, 1)./Types(:, 2))./(1-normcdf(-Types(:, 1)./Types(:, 2))) - ...
    (normpdf(-Types(:, 1)./Types(:, 2))./(1-normcdf(-Types(:, 1)./Types(:, 2)))).^2));

sqrt(simu_var)
simu_mean



















[teste3_a, teste3_b] = meshgrid((8000-24*20):20:(8000+25*20), ...
    (teste2_b(37,5)-24*50):50:(teste2_b(37,5)+25*50));

teste3_mean = zeros(50,50);
teste3_var = zeros(50,50);

for k=1:50
    for l=1:50
         typeDistributionMean = ...
        [1*10^(-5), 1330, teste3_a(1,k), teste3_b(l,1)]; % Original A was 1.9*10^-3
        typeDistributionLogCovariance = ...
            [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
             -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
             -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
              0     0     0    0.25]; % ???

        costOfPublicFunds = 0;

        % Calculation parameters
        populationSize = 5*10^2;
        slopeVector{1}            = 0:0.01:1;

        innerTypeDistributionLogCovariance = typeDistributionLogCovariance;

        Model = healthcaralognormalmodel_censnorm(slopeVector{1}, typeDistributionMean, innerTypeDistributionLogCovariance);

        Population = population(Model, populationSize);

        Types = zeros(populationSize, 2);

        for i=1:populationSize
           Types(i,:) = [Population.typeList{i}.M, Population.typeList{i}.S];
        end

        teste3_mean(k,l) = mean(Types(:, 1) - normpdf(-Types(:, 1)./Types(:, 2))./...
            (1-normcdf(-Types(:, 1)./Types(:, 1))).*Types(:, 1));

        teste3_var(k,l) = mean(Types(:, 2).^2.*(1 - Types(:, 1)./Types(:, 2).*...
            normpdf(-Types(:, 1)./Types(:, 2))./(1-normcdf(-Types(:, 1)./Types(:, 2))) - ...
            (normpdf(-Types(:, 1)./Types(:, 2))./(1-normcdf(-Types(:, 1)./Types(:, 2)))).^2));
    end
    disp(k);
end

surf(sqrt(teste3_var))
surf(teste3_mean)