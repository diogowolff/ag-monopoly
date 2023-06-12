[~,~,~,eq_demand]=Population.demand(pEquilibrium);
meanLossVec = zeros(Population.size,1);
stdLossVec = zeros(Population.size,1);
riskAvVec = zeros(Population.size,1);	
mhVec = zeros(Population.size,1);

for i = 1 : Population.size
	if eq_demand(i) > 0
		meanLossVec(i) = Population.typeList{i}.M;
		stdLossVec(i) =  Population.typeList{i}.S;
    		riskAvVec(i) = Population.typeList{i}.A;
		mhVec(i) = Population.typeList{i}.H;
	end;	
end;



mean(meanLossVec + normpdf(-meanLossVec./stdLossVec)./...
            (1-normcdf(-meanLossVec./stdLossVec)).*stdLossVec)


mean(sqrt(stdLossVec.^2.*(1 - meanLossVec./stdLossVec.*...
    normpdf(-meanLossVec./stdLossVec)./(1-normcdf(-meanLossVec./stdLossVec)) - ...
    (normpdf(-meanLossVec./stdLossVec)./(1-normcdf(-meanLossVec./stdLossVec))).^2)))


mean(meanLossVec + riskAvVec.*stdLossVec.^2/2 + ...
    (log((1 - normcdf(-meanLossVec./stdLossVec - riskAvVec.*stdLossVec))./...
    (1 - normcdf(-meanLossVec./stdLossVec))))./riskAvVec)



Avec = 10^-7:10^-7:10^-3;

calculateMeanRisk = @(riskAvVec) ...
    mean(meanLossVec + riskAvVec .* stdLossVec.^2/2 + ...
    (log((1 - normcdf(-meanLossVec./stdLossVec - riskAvVec .* stdLossVec))./...
    (1 - normcdf(-meanLossVec./stdLossVec))))./riskAvVec);

result = calculateMeanRisk(Avec);

[~, index] = min(abs(result - originalCE));
closestValue = Avec(index);

disp(closestValue);

calculateMeanRisk(riskAvVec)





mean(meanLossVec + normpdf(-meanLossVec./stdLossVec)./...
            (normcdf(1000000-meanLossVec./stdLossVec)-normcdf(-meanLossVec./stdLossVec)).*stdLossVec)


mean(meanLossVec + riskAvVec.*stdLossVec.^2/2 + ...
    (log((normcdf(1000000-meanLossVec./stdLossVec) - normcdf(-meanLossVec./stdLossVec - riskAvVec.*stdLossVec))./...
    (normcdf(1000000-meanLossVec./stdLossVec) - normcdf(-meanLossVec./stdLossVec))))./riskAvVec)








[~,~,Equilibrium_CS,~] = demand(Population, pEquilibrium)

Eq_profit = profit(Population, pEquilibrium)

[D, TC] = Population.demand(pEquilibrium);

[~,~,MaxWelf_CS,~] = demand(Population, pWelfare)

Welf_profit = profit(Population, pWelfare)

[~,~,Monop_CS,~] = demand(Population, pEfficient)

Monop_profit = profit(Population, pEfficient)






























for i = 1 : Population.size
		meanLossVec(i) = Population.typeList{i}.k;
		stdLossVec(i) =  Population.typeList{i}.theta;
    		riskAvVec(i) = Population.typeList{i}.A;
		mhVec(i) = Population.typeList{i}.H;	
end;
aux = sum(imag(Population.uMatrix) ~= 0, 2) == 0;
kVec = meanLossVec(aux);
thetaVec = stdLossVec(aux);
mean(kVec.*thetaVec)
mean(sqrt(kVec.*thetaVec.^2))
riskAvVec = riskAvVec(aux);

mean(riskAvVec)

mean(meanLossVec.*stdLossVec)
mean(sqrt(meanLossVec.*stdLossVec.^2))
mean(-meanLossVec.*log(1 - riskAvVec.*stdLossVec)./riskAvVec)

figure
pdSix = fitdist(riskAvVec','Kernel','Width',4);
x = 0:.1:45;
ySix = pdf(pdSix,x);
plot(x,ySix,'k-','LineWidth',2)




meanS = sqrt(25000^2 - 5100^2);
originalCE = 4340 + 10^-5*meanS^2/2;

x = 0:50:30000;
y1 = gampdf(x,.03,209000);
plot(x,y1)