clear;
close all;
addpath('../classes');
addpath('plotFunctions');
addpath('../plotFunctions/export_fig_updated');
addpath('../plotFunctions');
rng(1);


%% Competitive equilibrium vs monopoly histogram

% Start
    modelNameString = 'interval_more_contracts';
    
    Interval = load(modelNameString); %%%%%% This has to be changed too
    
% Calculate necessary series
    nContracts = Interval.Model.nContracts;
    xGrid = zeros(1, Interval.Model.nContracts);
    for j = 1 : nContracts
        xGrid(j) = Interval.Model.contracts{j}.slope;
    end;

% Histograms
    figMonopQuantities = figure;
    set(figMonopQuantities, 'Position', get(0, 'Screensize'));
    set(figMonopQuantities, 'name', 'Quantities with and without a Monopoly', 'numbertitle', 'on');
    % Plot quantities
        bar(xGrid, ...
            nContracts .* ...
            [Interval.DEquilibrium; ...
            Interval.DEfficient]', 3);
        hold on;
        % axis([0, 1, 0, .1 .* nContracts]);
    % Labels
        legend('Competitive Equilibrium', 'Monopoly', 'Location', 'NorthEast');
        legend boxoff
        xlabel('Coverage');
        ylabel('(Density)');
   % Other options
        box off;
        set(gca,'FontSize',27);
        set(findall(gcf,'type','text'),'FontSize',27);
        
   % Save
       fileName = '../figures40m/competition_vs_monopoly_quantities_more_contracts.pdf';        
       export_fig(fileName, '-transparent');
       fileName = [fileName(1: length(fileName)-4), '.eps'];
       print(fileName, '-depsc2');
    
%% Surplus-maximizing and monopoly histogram

% Calculate necessary series
    nContracts = Interval.Model.nContracts;
    xGrid = zeros(1, Interval.Model.nContracts);
    for j = 1 : nContracts
        xGrid(j) = Interval.Model.contracts{j}.slope;
    end;

% Histograms
    figSurplusQuantities = figure;
    set(figSurplusQuantities, 'Position', get(0, 'Screensize'));
    set(figSurplusQuantities, 'name', 'Surplus-maximizing and monopoly quantities', 'numbertitle', 'on');
    % Plot quantities
        bar(xGrid, ...
            nContracts .* ...
            [Interval.DWelfare; ...
            Interval.DEfficient]', 3);
        hold on;
        % axis([0, 1, 0, .1 .* nContracts]);
    % Labels
        legend('Maximum Surplus', 'Monopoly', 'Location', 'NorthEast');
        legend boxoff
        xlabel('Coverage');
        ylabel('(Density)');
   % Other options
        box off;
        set(gca,'FontSize',27);
        set(findall(gcf,'type','text'),'FontSize',27);
        
   % Save
       fileName = '../figures40m/welfare_vs_monopoly_quantities_more_contracts.pdf';        
       export_fig(fileName, '-transparent');
       fileName = [fileName(1: length(fileName)-4), '.eps'];
       print(fileName, '-depsc2');
    
%% Monopoloy and equilibrium prices

% Parameters
populationSize           = 10^5;
nPopulationDemandProfile = min(7000, populationSize);

% Color Scheme
colorRed  = [230,77,79]   / 255;
colorBlue = [208,226,241]   / 255;
colorBlueLine = colorBlue     / 2;

% Prepare
xGrid = zeros(1, Interval.Model.nContracts);
for j = 1 : Interval.Model.nContracts
    xGrid(j) = Interval.Model.contracts{j}.slope;
end;

Population = population(Interval.Model, populationSize);

% Calculate mean loss parameter of agents purchasing each contract
meanLossTypeVector = zeros(Population.size, 1);
riskAversionVector = zeros(Population.size, 1);
for i = 1 : Population.size
    meanLossTypeVector(i) = Population.typeList{i}.M;
    riskAversionVector(i) = Population.typeList{i}.A;
end;

[~, ~, ~, choiceVectorEquilibrium] = Population.demand(Interval.pEquilibrium);
[~, ~, ~, choiceVectorEfficient]   = Population.demand(Interval.pEfficient);
meanLossEquilibrium     = zeros(1, Interval.Model.nContracts);
meanLossEfficient       = zeros(1, Interval.Model.nContracts);
for j = 1 : Interval.Model.nContracts
    buyersEquilibrium = (choiceVectorEquilibrium == j);
    meanLossVectorBuyersEquilibrium = meanLossTypeVector(buyersEquilibrium);
    meanLossEquilibrium(j)  = mean(meanLossVectorBuyersEquilibrium);
    
    buyersEfficient = choiceVectorEfficient == j;
    meanLossVectorBuyersEfficient = meanLossTypeVector(buyersEfficient);
    meanLossEfficient(j)  = mean(meanLossVectorBuyersEfficient);
end;

figEffnEqPrices = figure;
    set(figEffnEqPrices, 'Position', get(0, 'Screensize'));
    set(figEffnEqPrices, 'name', 'Monopoly and Equilibrium Price and Mean Loss', 'numbertitle', 'on');
    % Equilibrium
    % Plot prices
        plot(xGrid, Interval.pEquilibrium, 'k', 'linewidth', 3);
        hold on;
    % Plot average loss parameter
        scatter(xGrid, meanLossEquilibrium, max(1,Interval.DEquilibrium.*3600), ...
            'MarkerFaceColor', colorRed, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.5); 
    % Efficient
    % Plot prices
        plot(xGrid, Interval.pEfficient, 'color', colorBlueLine, 'linewidth', 3);
        hold on;
    % Plot average loss parameter
        scatter(xGrid, meanLossEfficient, max(1,Interval.DEfficient.*3600), ...
            'MarkerFaceColor', colorBlue, ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 1.5);        
    % Labels
        legend('Equilibrium Prices', 'Equilibrium Losses', ...
               'Monopoly Prices', 'Monopoly Losses', ...
               'Location', 'NorthWest');
        legend boxoff;
        xlabel('Contract');
        ylabel('($)');
   % Other options
        box off;
        set(gca,'FontSize',27);
        set(findall(gcf,'type','text'),'FontSize',27);
   % Financial tick label
        axis([0,1,0,max(Interval.pEfficient)])
        set(gca, 'YTickLabel', num2bank(get(gca, 'YTick')));
        
fileName = ['../figures40m/', modelNameString, '_', 'monopoly_and_equilibrium_prices.pdf'];        
export_fig(fileName, '-transparent');
fileName = [fileName(1: length(fileName)-4), '.eps'];
print(fileName, '-depsc2');


%% Exclusion region

[~, ~, ~, choiceVectorMonopolist]   = Population.demand(Interval.pEfficient);
meanLossMonopolist       = zeros(1, Interval.Model.nContracts);

for j = 1 : Interval.Model.nContracts
    buyersMonopolist = choiceVectorMonopolist == j;
    meanLossVectorBuyersMonopolist = meanLossTypeVector(buyersMonopolist);
    meanLossMonopolist(j)  = mean(meanLossVectorBuyersMonopolist);
end;

slopeVectorMonopolist = zeros(populationSize, 1);
for j = 1 : Interval.Model.nContracts
    I = (choiceVectorMonopolist == j);
    slopeVectorMonopolist(I) = Interval.Model.contracts{j}.slope;
end;

minContract = min(slopeVectorMonopolist);

exclusionVector = slopeVectorMonopolist;
exclusionVector(exclusionVector>0) = [1];

figExclusion = figure;
    set(figExclusion, 'Position', get(0, 'Screensize'));
    set(figExclusion, 'name', 'Profit-maximizing Exclusion Region', 'numbertitle', 'on');
    % Plot scatter
        gscatter(meanLossTypeVector(1:nPopulationDemandProfile), ...
            riskAversionVector(    1:nPopulationDemandProfile), ...
            exclusionVector(1:nPopulationDemandProfile), ...
        [0 0.4470 0.7410; 0.9290 0.6940 0.1250]);
    % Labels
        xlabel('Average Loss,  M_{\theta}');
        ylabel('Risk Aversion, A_{\theta}');
        legend('No coverage', 'Some coverage');
    % Properties
        set(gca, 'xscale', 'log', ...
                 'yscale', 'log');
        grid off;
   % Other options
        box off;
        set(gca,'FontSize',27);
        set(findall(gcf,'type','text'),'FontSize',27);        
    % Financial tick label
                axis([500,100000,10^(-6),10^(-4)])
        set(gca, 'XTickLabel', num2bank(get(gca, 'xtick')));     

fileName = ['../figures40m/', modelNameString, '_', 'exclusion.pdf'];        
export_fig(fileName, '-transparent');
fileName = [fileName(1: length(fileName)-4), '.eps'];
print(fileName, '-depsc2');


%% Table of exclusion and percentiles

% clear;
% modelNameString = 'interval_more_contracts';
% load(modelNameString);
% 
% table = [DEquilibrium(1), sum(DEquilibrium(19:26)), sum(DEquilibrium(21:26)), ...
%     sum(DEquilibrium(24:26)); DEfficient(1), sum(DEfficient(19:26)), ...
%     sum(DEfficient(21:26)), sum(DEfficient(24:26)); DWelfare(1), ...
%     sum(DWelfare(19:26)), sum(DWelfare(21:26)), sum(DWelfare(24:26))];
% 
% save('tails_of_insurance_levels_more.mat', 'table');