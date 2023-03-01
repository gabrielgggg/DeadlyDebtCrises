function [ oout ] = modelVersusModel(mmodel, rref, param)
cWindow = mmodel.yyear >= 2020 & mmodel.yyear <= 2024;

ce = @(vvv) (1 + (1 - param.sigma) * (1 - param.bet) * vvv).^(1/(1 - param.sigma));
ceGain = @(cNew, cOld) (cNew - cOld) / (1 - param.bet) / 52;
qFromSp = @(sss) 1 ./ ( sss / param.pay + 1 );

oout = NaN([6, 1]);
oout(1) = 100 * ceGain(ce(mmodel.Vpath(param.before)), ce(rref.Vpath(param.before)));
oout(2) = sum(rref.Dpath > 0) / 4 - sum(mmodel.Dpath > 0) / 4;
oout(3) = max(rref.spPath) - max(mmodel.spPath);
oout(4) = 100 * (1 - max(mmodel.muDpath) ./ max(rref.muDpath));
oout(5) = 100 * ( sum(mmodel.Cpath(cWindow)) ./ sum(rref.Cpath(cWindow)) - 1);

spOnImpact = (1 + mmodel.spPath(param.before) / 100).^(1/52) - 1;
spOnImpactRef = (1 + rref.spPath(param.before) / 100).^(1/52) - 1;
oout(6) = (qFromSp(spOnImpact) - qFromSp(spOnImpactRef)) * mmodel.Bpath(param.before) ./ 52 * 100;
end

