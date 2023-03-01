function [ oout ] = modelPerfectOutcomes(pmodel, param)
% tmpt = find(isnan(pmodel.dataMuD)) - 1;
% lastDay = tmpt(1);

oout = NaN([16, 1]);

oout(1) = 100 * max(pmodel.muD);
oout(2) = sum(pmodel.muI > 1.0e-6) / 4;

oout(3) = max(pmodel.B - pmodel.cfB) ./ 52 * 100;
oout(4) = 0.0;
oout(5) = 0.0;

oout(6) = -sum(pmodel.Y - 1) ./ 52 * 100;
oout(7) = -sum((1 - param.thetaY * pmodel.L).^0.67 - 1) ./ 52 * 100;
oout(8) = 0.0;
oout(9) = 0.0;

oout(10:14) = 100 * (pmodel.C(1) - pmodel.cfC(1));

oout(15) = -100 * (param.ce - param.ceNoCOVID) / (1 - param.beta) / 52;
oout(16) = 0.0;
end

