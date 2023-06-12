[teste2_a, teste2_b] = meshgrid(0:0.001:1, 130000:50:160000);

teste2_mean = teste2_a.*teste2_b;
teste2_sd = sqrt(teste2_a.*teste2_b.^2);


meanS = sqrt(25000^2 - 5100^2);
teste2_sddif = (teste2_sd - meanS);
teste2_meandif = (teste2_mean - 4340);

a = abs(teste2_sddif) < 25;
b = abs(teste2_meandif) < 50;
c = logical(a.*b);
teste = teste2_a(c);
teste2 = teste2_b(c);
testeidk = [teste, teste2];

% Extract the coordinates from testeidk
x = testeidk(:,1);
y = testeidk(:,2);

z_mean = interp2(teste2_a, teste2_b, teste2_mean, x, y);
z_sd = interp2(teste2_a, teste2_b, teste2_sd, x, y);


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
scatter3(x, y, z_mean, 'r', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);

subplot(1,2,2);
surf(teste2_a(1,:), teste2_b(:,1),teste2_sd);
title('Standard Deviation');
xlabel('a');
ylabel('b');
zlabel('Standard Deviation');

% Add the scatter plot to the second plot
hold on;
scatter3(x, y, z_sd, 'r', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);




%fixing params
k = 0.031;
theta = 139000;



% matching risk aversion

originalCE = 4340 + 10^-5*meanS^2/2;

ub = 1/theta;
Avec = 0:0.0000001:ub;

CE = -k.*log(1-Avec.*theta)./Avec;

newA = fzero(@(x) -k.*log(1-x.*theta)./x - originalCE, [.1*ub/2, 9*ub/10]);

1-newA*theta





%--------------------------------- ex-post version
addpath('../Classes');
[teste2_a, teste2_b] = meshgrid(.015:0.0005:.4, 150000:250:210000);

teste2_mean = zeros(771,241);
teste2_sd = zeros(771,241);

parfor k=1:771
    for l=1:241
         typeDistributionMean = ...
        [1*10^(-5), 1330, teste2_a(1,k), teste2_b(l,1)]; % Original A was 1.9*10^-3
        typeDistributionLogCovariance = ...
            [ 0.25 -0.01 -0.12 0    ; % c11 = 0.25 originally
             -0.01  0.28 -0.03 0    ; % c22 = 0.98 originally
             -0.12 -0.03  0.20 0    ; % c33 = 0.20 originally
              0     0     0    0.25]; % ???

        costOfPublicFunds = 0;

        % Calculation parameters
        populationSize = 1*10^3;
        slopeVector            = 0:0.04:1;

        innerTypeDistributionLogCovariance = typeDistributionLogCovariance;

        Model = healthcaralognormalmodel_gamma(slopeVector, typeDistributionMean, innerTypeDistributionLogCovariance);
        
        Population = population(Model, populationSize);
        populationSize = length(Population.typeList);
        Types = zeros(populationSize, 2);

        for i=1:populationSize
           Types(i,:) = [Population.typeList{i}.k, Population.typeList{i}.theta];
        end

        teste2_mean(k,l) = mean(Types(:, 1).*Types(:, 2));

        teste2_sd(k,l) = mean(sqrt(Types(:, 1).*Types(:, 2).^2));

    end
end

meanS = sqrt(25000^2 - 5100^2);

idx_mean = find(abs(teste2_mean - 4340) <= 500);
idx_sd = find(abs(teste2_sd - meanS) <= 250);

idx_intersection = intersect(idx_mean, idx_sd);

z_mean = teste2_mean(idx_intersection);
z_sd = teste2_sd(idx_intersection);

% Extract the coordinates from testeidk
t2t = teste2_a';
t2t2 = teste2_b';
x = t2t(idx_intersection);
y = t2t2(idx_intersection);

% Generate the surf plots
figure;
subplot(1,2,1);
surf(teste2_a, teste2_b, teste2_mean');
title('Mean');
xlabel('a');
ylabel('b');
zlabel('Mean');
hold on;

% Add the scatter plot to the first plot
scatter3(x, y, z_mean, 'r', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);

subplot(1,2,2);
surf(teste2_a, teste2_b,teste2_sd');
title('Standard Deviation');
xlabel('a');
ylabel('b');
zlabel('Standard Deviation');

% Add the scatter plot to the second plot
hold on;
scatter3(x, y, z_sd, 'r', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);


% 0.03 k and 209000 theta seems good ex-post, could maybe look later

k=0.03;
theta = 209000;

% matching risk aversion

originalCE = 4340 + 10^-5*meanS^2/2;

ub = 1/theta;
Avec = 0:0.0000001:ub;

CE = -k.*log(1-Avec.*theta)./Avec;

newA = fzero(@(x) -k.*log(1-x.*theta)./x - originalCE, [.1*ub/2, 9*ub/10]);

