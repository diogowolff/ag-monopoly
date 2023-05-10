clear;
addpath('../Classes');
%chuta um par, calcula, pega as vars, faz o match, etc etc


%very wide range to find parameters

[teste2_a, teste2_b] = meshgrid(-370000:200:-340000, 45000:250:60000);

teste2_mean = zeros(151,61);
teste2_var = zeros(151,61);

parfor k=1:151
    for l=1:61
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
        slopeVector            = 0:0.01:1;

        innerTypeDistributionLogCovariance = typeDistributionLogCovariance;

        Model = healthcaralognormalmodel_censnorm(slopeVector, typeDistributionMean, innerTypeDistributionLogCovariance, 0);

        Population = population(Model, populationSize);

        Types = zeros(populationSize, 2);

        for i=1:populationSize
           Types(i,:) = [Population.typeList{i}.M, Population.typeList{i}.S];
        end

        teste2_mean(k,l) = mean(Types(:, 1) + normpdf(-Types(:, 1)./Types(:, 2))./...
            (1-normcdf(-Types(:, 1)./Types(:, 2))).*Types(:, 2));

        teste2_var(k,l) = mean(Types(:, 2).^2.*(1 - Types(:, 1)./Types(:, 2).*...
            normpdf(-Types(:, 1)./Types(:, 2))./(1-normcdf(-Types(:, 1)./Types(:, 2))) - ...
            (normpdf(-Types(:, 1)./Types(:, 2))./(1-normcdf(-Types(:, 1)./Types(:, 2)))).^2));

        teste2_aa(k, l) = mean(Types(:, 1));
        teste2_bb(k, l) = mean(Types(:, 2));
    end
end

meanS = sqrt(25000^2 - 5100^2);

teste2_vardif = (sqrt(teste2_var) - meanS);
teste2_meandif = (teste2_mean - 4340);

teste2 = sqrt(teste2_var);

teste2_sumsq = teste2_vardif.^2 + teste2_meandif.^2;
[row, col] = find(teste2_sumsq==min(teste2_sumsq(:)));

figure;
subplot(1,2,1);
surf(teste2_mean);
subplot(1,2,2);
surf(teste2);

sqrt(teste2_var(row, col))
teste2_mean(row, col)

teste2_a(77,34)
teste2_b(77,34)

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







[teste2_a, teste2_b] = meshgrid(-380000:500:-350000, 40000:50:70000);

teste2_mean = teste2_a + normpdf(-teste2_a./teste2_b)./...
    (1-normcdf(-teste2_a./teste2_b)).*teste2_b;

teste2_sd = sqrt(teste2_b.^2.*(1 - teste2_a./teste2_b.*...
            normpdf(-teste2_a./teste2_b)./(1-normcdf(-teste2_a./teste2_b)) - ...
            (normpdf(-teste2_a./teste2_b)./(1-normcdf(-teste2_a./teste2_b))).^2));
teste2_sd(imag(teste2_sd)~=0) = nan;

meanS = sqrt(25000^2 - 5100^2);
teste2_sddif = (teste2_sd - meanS);
teste2_meandif = (teste2_mean - 4340);

a = abs(teste2_sddif) < 250;
b = abs(teste2_meandif) < 750;
c = logical(a.*b);
teste = teste2_a(c);
teste2 = teste2_b(c);
testeidk = [teste, teste2];

% Extract the coordinates from testeidk
x = testeidk(:,1);
y = testeidk(:,2);

% Generate the surf plots
figure;
subplot(1,2,1);
surf(teste2_a(1,:), teste2_b(:,1), teste2_mean);
title('Mean');
xlabel('a');
ylabel('b');
zlabel('Mean');
hold on;

% Add the scatter plot to the first plot
scatter3(x, y, z_mean, 'r', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.8);

subplot(1,2,2);
surf(teste2_a(1,:), teste2_b(:,1),teste2_sd);
title('Standard Deviation');
xlabel('a');
ylabel('b');
zlabel('Standard Deviation');

% Add the scatter plot to the second plot
hold on;
scatter3(x, y, z_sd, 'r', 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.8);

