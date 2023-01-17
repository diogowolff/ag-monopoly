clear;
addpath('../classes');
rng(1);

table  = zeros(3, 2);


%Table
%                                   Interval |  Interval high MH        
%Consumer surplus in equilibrium  |          |
%Consumer surplus with monopoly   |          |
%Monopolist profit                |          |
%

Interval = load('interval');
IntervalHighMH = load('interval_high_mh_variance');

Population = population(Interval.Model, 500);
[~, ~, CS_eq, ~] = Population.demand(Interval.pEquilibrium);
[~, ~, CS_monop, ~] = Population.demand(Interval.pEfficient);
Profit = Population.profit(Interval.pEfficient);

table(:, 1) = [ ...
    CS_eq, ...
    CS_monop, ...
    Profit];


Population = population(IntervalHighMH.Model, 500);
[~, ~, CS_eq, ~] = Population.demand(IntervalHighMH.pEquilibrium);
[~, ~, CS_monop, ~] = Population.demand(IntervalHighMH.pEfficient);
Profit = Population.profit(IntervalHighMH.pEfficient);

table(:, 2) = [ ...
    CS_eq, ...
    CS_monop, ...
    Profit];

save('table');