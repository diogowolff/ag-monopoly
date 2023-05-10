alpha = (0-meanLossVec)./riskAvVec;

left_expec = (meanLossVec.*riskAvVec + (riskAvVec .* stdLossVec).^2 ./ 2) + ... 
                    log((1 - normcdf(alpha - stdLossVec .* riskAvVec))./(1 - normcdf(alpha)));

left_term = exp(meanLossVec.*riskAvVec + (riskAvVec .* stdLossVec).^2 ./ 2);
right_term = (1 - normcdf(alpha - stdLossVec .* riskAvVec))./(1 - normcdf(alpha));

inside_left = meanLossVec.*riskAvVec + (riskAvVec .* stdLossVec).^2 ./ 2;
inside_rt = (riskAvVec .* stdLossVec).^2 ./ 2;
inside_lt = meanLossVec.*riskAvVec;