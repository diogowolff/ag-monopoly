function [ ] = plotMonopolist(...
    Model, pMonopolist, modelNameString)

rng(1);

% Parameters
populationSize           = 10^5;
nPopulationDemandProfile = min(7000, populationSize);

% Color Scheme

% Prepare
xGrid = zeros(1, Model.nContracts);
for j = 1 : Model.nContracts
    xGrid(j) = Model.contracts{j}.slope;
end;

Population = population(Model, populationSize);

% Calculate mean loss parameter of agents purchasing each contract
meanLossTypeVector = zeros(Population.size, 1);
riskAversionVector = zeros(Population.size, 1);
for i = 1 : Population.size
    meanLossTypeVector(i) = Population.typeList{i}.M;
    riskAversionVector(i) = Population.typeList{i}.A;
end;

[~, ~, ~, choiceVectorMonopolist]   = Population.demand(pMonopolist);
meanLossMonopolist       = zeros(1, Model.nContracts);

for j = 1 : Model.nContracts
    buyersMonopolist = choiceVectorMonopolist == j;
    meanLossVectorBuyersMonopolist = meanLossTypeVector(buyersMonopolist);
    meanLossMonopolist(j)  = mean(meanLossVectorBuyersMonopolist);
end;

slopeVectorMonopolist = zeros(populationSize, 1);
for j = 1 : Model.nContracts
    I = (choiceVectorMonopolist == j);
    slopeVectorMonopolist(I) = Model.contracts{j}.slope;
end;

minContract = min(slopeVectorMonopolist);

figOptDemandProfile = figure;
    set(figOptDemandProfile, 'Position', [0 0 1.618.*500 500]);
    set(figOptDemandProfile, 'name', 'Profit-maximizing Demand Profile', 'numbertitle', 'on');
    % Plot scatter
        scatter(meanLossTypeVector(1:nPopulationDemandProfile), ...
            riskAversionVector(    1:nPopulationDemandProfile), 20, ...
            slopeVectorMonopolist(1:nPopulationDemandProfile));
        colorbar;
            caxis([minContract, 1]);
    % Labels
        xlabel('Average Loss,  M_{\theta}');
        ylabel('Risk Aversion, A_{\theta}');
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
        
fileName = ['../figuresMonopolist/', modelNameString, '_', 'demand_profile.pdf'];        
export_fig(fileName, '-transparent');
fileName = [fileName(1: length(fileName)-4), '.eps'];
print(fileName, '-depsc2');

exclusionVector = slopeVectorMonopolist;
exclusionVector(exclusionVector>0) = [1];

figExclusion = figure;
    set(figExclusion, 'Position', [0 0 1.618.*500 500]);
    set(figExclusion, 'name', 'Profit-maximizing Exclusion Region', 'numbertitle', 'on');
    % Plot scatter
        gscatter(meanLossTypeVector(1:nPopulationDemandProfile), ...
            riskAversionVector(    1:nPopulationDemandProfile), ...
            exclusionVector(1:nPopulationDemandProfile));

    % Labels
        xlabel('Average Loss,  M_{\theta}');
        ylabel('Risk Aversion, A_{\theta}');
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

fileName = ['../figuresMonopolist/', modelNameString, '_', 'exclusion.pdf'];        
export_fig(fileName, '-transparent');
fileName = [fileName(1: length(fileName)-4), '.eps'];
print(fileName, '-depsc2');

end

