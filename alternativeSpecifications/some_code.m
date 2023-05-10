[~,~,~,eq_demand]=Population.demand(pEquilibrium);
meanLossVec = zeros(Population.size,1);
stdLossVec = zeros(Population.size,1);
riskAvVec = zeros(Population.size,1);	
mhVec = zeros(Population.size,1);

for i = 1 : Population.size
	if eq_demand(i) > 0
		meanLossVec(i) = Population.typeList{i}.k;
		stdLossVec(i) =  Population.typeList{i}.theta;
    		riskAvVec(i) = Population.typeList{i}.A;
		mhVec(i) = Population.typeList{i}.H;
	end;	
end;

mean(meanLossVec.*stdLossVec)
mean(sqrt(meanLossVec.*stdLossVec.^2))
CE = -k.*log(1-Avec.*theta)./Avec;
mean(-meanLossVec.*log(1 - riskAvVec.*stdLossVec)./riskAvVec)

sum(1 - riskAvVec.*stdLossVec <= 0)

mean(normpdf(-meanLossVec./stdLossVec)./...
            (1-normcdf(-meanLossVec./stdLossVec)))

4340

mean(sqrt(stdLossVec.^2.*(1 - meanLossVec./stdLossVec.*...
            normpdf(-meanLossVec./stdLossVec)./(1-normcdf(-meanLossVec./stdLossVec)) - ...
            (normpdf(-meanLossVec./stdLossVec)./(1-normcdf(-meanLossVec./stdLossVec))).^2)))