classdef healthcaralognormalmodel_censnorm < model
    %healthcaralognormalmodel_censnorm Creates a health insurance model as in
    %the Azevedo and Gottlieb paper. Subclass of the model class.
    %   This models have CARA preferences, truncated normal losses, linear contracts,
    %   and lognormally distributed types. Inputs to the constructor:
    %   slopeVector is a vector of contract slopes, parameterMean is a 4
    %   dimensional vector of means of parameters and parameterLogVariance
    %   is a 4x4 matrix of log covariances. Parameters are ordered as A, H,
    %   M, S as in the Azevedo and Gottlieb paper (absolute risk aversion A
    %   , moral hazard H, mean loss M, and standard deviation of losses S).
    
    properties
        typeDistributionMean
        typeDistributionLogCovariance
    end
    
    methods
        % Constructor
        function Model = healthcaralognormalmodel_lognorm(slopeVector, ...
                typeDistributionMean, typeDistributionLogCovariance)
            
            Model.typeDistributionMean = typeDistributionMean;
            Model.typeDistributionLogCovariance = ...
                typeDistributionLogCovariance;
            
            n = length(slopeVector);
            for i = 1:n
                x.slope =         slopeVector(i) ;
                x.name  = num2str(slopeVector(i));
                Model.contracts{i} = x;
            end;
        end
        
        function x = lossPDF(~, type, l)
            % Probability distribution function for the losses
            qsi = (l-type.M)./type.S;
            alpha = -type.M./type.S;
            
            x = normpdf(qsi)./(type.S.*(1 - normcdf(alpha)));
        end
        
        function x = lossCDF(~, type, l)
            % Cumulative distribution function for the losses
            qsi = (l-type.M)./type.S;
            alpha = -type.M./type.S;
            
            x = (normcdf(qsi) - normcdf(alpha))./(1 - normcdf(alpha));
        end
        
        function u = uFunction(obj, contract, type)
            % Returns the willingness to pay of an agent for a contract,
            % No outside option in this case
                
                u = expectedValue( obj, @(l) ...
                    exPostUtility(obj, contract, type, l),  contract, type );
            
        end
            
       
        function c = cFunction(obj, contract, type)
            
                c = expectedValue( obj, @(l) ...
                    exPostCost(obj, contract, type, l), contract, type );
            
        end
        
        function e = eFunction(~, ~, ~)
                e = 0;
        end
            
        function type = typeDistribution(obj)
            % Returns a rando structure with the individual caracteristics
            % distributed according to the lognormal parameters
            
            v = lognrndfrommoments(...
                obj.typeDistributionMean, obj.typeDistributionLogCovariance, 1);
            
            type.A = v(1);
            type.H = v(2);
            type.M = v(3);
            type.S = v(4);
            
            function v = ...
                    lognrndfrommoments(meanVector, logCovMatrix, varargin)
                %   estimate_lognormal_from_moments Estimate lognormal parameters
                %   Inputs: a line vector of means, and a log covariance matrix.
                %   Outputs: A line vector of draws. Optional argument for number of lines.
                %   Warning: for some reason MATLAB chokes if one of the variables has 0
                %   variance, so always put at least a tiny bit of bogus variance.
                
                nDraws = 1;
                if (length(varargin) == 1)
                    nDraws = varargin{1};
                end;
                
                % Find lognormal parameters
                sigmaMatrix = cholcov(logCovMatrix);
                b  = sum(sigmaMatrix.^2, 1);
                mu = log(meanVector) - b/2;
                
                % Draw
                x = randn(nDraws, length(meanVector));
                v = exp(repmat(mu,[nDraws, 1]) + x * sigmaMatrix);
            end
            
        end
        
        function expenditure = exPostExpenditure(obj, contract, type, losses)
            
                [~, expenditure] = exPostUtility(obj, contract, type, losses);
            
        end
        
        function cost = exPostCost(obj, contract, type, losses)
            
                [~, expenditure, payment] = exPostUtility(obj, contract, type, losses);
                cost = expenditure - payment;
            
        end
        
        function [u, expenditure, payment, bounds] = exPostUtility(obj, contract, type, losses)
            
            % Initializing variables
            u = zeros(1, length(losses));
            expenditure = zeros(1, length(losses));
            payment = zeros(1, length(losses));
            bounds = zeros(length(losses), 3);
            
            for i = 1:length(losses)
                
                % Loss is bounded below by zero
                l = max(losses(i), 0);
                
                    
                    % Calculating the loss boundaries for each interval of expenses
                    bounds(i, 1) = max(min(contract.deductible-(1-contract.coinsurance)*type.H/2,contract.oopMax-type.H/2),0);
                    bounds(i, 2) = max(contract.deductible-(1-contract.coinsurance)*type.H/2,0);
                    bounds(i, 3) = max((contract.oopMax-(1-contract.coinsurance)*contract.deductible)/contract.coinsurance ...
                        - (2 - contract.coinsurance) * type.H / 2, 0);
                    
                    if l < bounds(i, 1)
                        
                        expenditure(i) = l;
                        u(i) = -l;
                        payment(i) = l;
                        
                    elseif (l >= bounds(i, 2)) && (l < bounds(i, 3))
                        
                        u(i) = (1-contract.coinsurance)^2*type.H/2 - (1-contract.coinsurance)*contract.deductible - contract.coinsurance*l;
                        expenditure(i) = (1-contract.coinsurance)*type.H + l;
                        payment(i) = contract.deductible + contract.coinsurance*(expenditure(i)-contract.deductible);
                        
                    else
                        
                        u(i) = type.H/2 - contract.oopMax;
                        expenditure(i) = type.H + l;
                        payment(i) = contract.oopMax;
                        
                    end
            end
        end
        
        function [populationSize, CalculationParametersEquilibrium, CalculationParametersOptimum] = ...
                suggestComputationParameters(Model, percentError)
            
            nContracts = Model.nContracts;
            Population = population(Model, 100);
            priceOrderOfMagnitude = mean(Population.cMatrix(:));
            
            populationSize = ...
                floor(nContracts / percentError^2 / 2);
            CalculationParametersEquilibrium.behavioralAgents = ...
                percentError;
            CalculationParametersEquilibrium.fudge = ...
                percentError / nContracts / 100;
            CalculationParametersEquilibrium.tolerance = ... 
                percentError * priceOrderOfMagnitude;
            CalculationParametersEquilibrium.maxIterations = ...
                floor(10 / CalculationParametersEquilibrium.fudge);
            
            CalculationParametersOptimum.tolerance = CalculationParametersEquilibrium.tolerance;
            CalculationParametersOptimum.maxIterations = 10^4;
        end;
    end
    
    
    methods ( Access = private, Hidden = true ) % Auxiliary methods
        
        
        function x = expectedValue(obj, function_handle, contract, type )
            
            
            [limits, oopMaxLoss] = integralLimits(obj, contract, type);
            
            x = leftIntegral(obj, function_handle, contract, type, limits );
            
            x = x + innerIntegral (obj, @(l) lossPDF(obj, type, l)...
                .* function_handle(l), contract, type, limits, oopMaxLoss );
            
            x = x + rightIntegral(obj, @(l) function_handle(l), contract,...
                type, limits, oopMaxLoss );
            
        end
        
        function x = leftIntegral(obj, function_handle, ~, type, limits )
            
            x = 0;
            
            if ( limits(1) == 0 )
                x = lossCDF(obj, type, 0) * function_handle(0);
            end
            
        end
        
        function x = innerIntegral(obj, function_handle, contract, type, limits, oopMaxLoss )
            
            x = 0;
            
            [~, ~, ~, bounds] = exPostUtility(obj, contract, type, 0);
            
            if ( limits(2) > 0 || limits(1) <  oopMaxLoss )
                x  = integral(@(l) function_handle(l), limits(1), limits(2),...
                    'AbsTol', 1e-15,'RelTol',1e-12,'WayPoints',...
                    [bounds(isfinite(bounds)),linspace(limits(1), limits(2),1e3)] );
            end
            
        end
        
        function x = rightIntegral(obj, function_handle, contract, type, limits, oopMaxLoss )
            
            x = 0;
            
            if ( limits(2) == oopMaxLoss )
                
                inclination = function_handle( oopMaxLoss + 1 )...
                    - function_handle( oopMaxLoss );
                
                x = (1 - lossCDF(obj, type, oopMaxLoss)) ...
                    * function_handle( oopMaxLoss );
                
                if (inclination > 0)
                    
                    x = x + inclination...
                        * ( type.M - innerIntegral(obj,@(l) l.*lossPDF(obj, type, l),...
                        contract, type, limits, oopMaxLoss) -  (1 - lossCDF(obj, type, oopMaxLoss)) ...
                        * oopMaxLoss  );
                    
                end
                
            end
            
        end
        
        function [limits, oopMaxLoss] = integralLimits(obj, contract, type)
            
            [~, ~, ~, bounds] = exPostUtility(obj, contract, type, 0);
            oopMaxLoss = max(bounds(1),bounds(3));
            
            if (lossPDF(obj,type,0) > 0);
                limits(1) = 0;
            elseif ( (type.M > oopMaxLoss ) && ...
                    lossPDF(obj, type, oopMaxLoss) == 0)
                limits(1) = oopMaxLoss;
            else
                limits(1) = max(0,findCloserNonZero(min(type.M,oopMaxLoss),type.M/100,1e-10));
            end
            
            if (lossPDF(obj,type,oopMaxLoss) > 0);
                limits(2) = oopMaxLoss;
            elseif (type.M > oopMaxLoss )
                limits(2) = oopMaxLoss;
            else
                limits(2) = max(0, findCloserNonZero(type.M,-type.M/100,1e-10));
            end
            
            function b = findCloserNonZero(b,d_init,tol)
                f_b = 0;
                d = d_init;
                while ( abs(d) > tol || f_b == 0 )
                    b = b - d;
                    f_b = obj.lossPDF(type,b);
                    if (f_b > 0)
                        if (sign(d) ~= sign(d_init))
                            d = - d / 10;
                        end
                    else
                        if (sign(d) == sign(d_init))
                            d = - d / 10;
                        end
                    end
                end
            end
        end
    end
end

