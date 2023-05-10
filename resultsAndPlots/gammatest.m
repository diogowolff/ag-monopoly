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
