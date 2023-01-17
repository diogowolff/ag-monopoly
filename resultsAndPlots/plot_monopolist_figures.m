clear;
close all;
addpath('../classes');
addpath('plotFunctions');
addpath('../plotFunctions/export_fig_updated');
addpath('./plotFunctions');
rng(1);

for modelNameString = { ...
        'interval',                  ...
        'interval_high_mh_variance'};
    
    CalculationData = load(modelNameString{1});
    plotMonopolist( ...
        CalculationData.Model, ...
        CalculationData.pEfficient, modelNameString{1});
    
end; 


% Start
    close all;
    clear;
    Interval = load('interval');
    
% Calculate necessary series
    nContracts = Interval.Model.nContracts;
    xGrid = zeros(1, Interval.Model.nContracts);
    for j = 1 : nContracts
        xGrid(j) = Interval.Model.contracts{j}.slope;
    end;

% Histograms
    figMonopQuantities = figure;
    set(figMonopQuantities, 'Position', [0 0 1.618.*500 500]);
    set(figMonopQuantities, 'name', 'Quantities with and without a Monopoly', 'numbertitle', 'on');
    % Plot quantities
        bar(xGrid, ...
            nContracts .* ...
            [Interval.DEquilibrium; ...
            Interval.DEfficient]', 3);
        hold on;
        axis([0, 1, 0, .1 .* nContracts]);
    % Labels
        legend('Competitive Equilibrium', 'Monopoly', 'Location', 'NorthWest');
        legend boxoff
        xlabel('Coverage');
        ylabel('(Density)');
   % Other options
        box off;
        set(gca,'FontSize',27);
        set(findall(gcf,'type','text'),'FontSize',27);
        
   % Save
       fileName = '../figuresEquilibrium/competition_vs_monopoly_quantities.pdf';        
       export_fig(fileName, '-transparent');
       fileName = [fileName(1: length(fileName)-4), '.eps'];
       print(fileName, '-depsc2');
    