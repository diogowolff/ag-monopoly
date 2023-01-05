# Extension of the code for "Perfect Competition in Markets with Adverse Selection" by Eduardo Azevedo and Daniel Gottlieb.

Generates in MATLAB the case of a monopolist in the economy setting of the paper. 


The scripts are divided in five folders:

## `./classes`

Contains classes with the most important functions.

- `model.m`: abstract class that defines all properties and methods that should be implemented by a model of perfect competition with adverse selection.
	- `healthcaralognormalmodel.m`: subclass for the health insurance model in the paper with linear contracts and normally distributed losses.
	- `healthcaralognormalmodel_nl`: subclass for the health insurance model in the paper with nonlinear contracts and lognormally distributed losses.
- `population.m`: A population object describes the preferences of a finite number of consumers. This class has methods for finding competitive equilibria and optima.

## `./resultsAndPlots`

Contains the code that runs the functions contained in /classes given the parameters specified, and plots results.

- `calculations_parallel.m`: creates the populations and runs the main calculations in parallel
- `plot_equilibrium_figures.m`: plots the figures. Uses `plotFunctions/plotEquilibriumAndOptimum.m` as an auxiliary function for some figures.
- `plot_monopolist_figures.m`: plots the figures specific to monopoly. Uses `plotFunctions/plotMonopolist.m` as an auxiliary function for some figures.

## `./resultsRemoteCode`

Contains the code from resultsAndPlots adapted to run in a machine with more cores; this isn't fully tested yet and i'm not sure if it runs faster than the former code.
The files are all the same except for remote_code_parallel instead of calculations_parallel.


## `./plotFunctions`

Auxiliary functions for plotting:

- `num2bank.m`: Downloaded from MATLAB central and slightly modified function for formatting numbers.
- `./export_fig`: https://github.com/altmany/export_fig

## License

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.


