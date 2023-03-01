function [ oout ] = modelOutcomes(mmodel, param)
window2020 = mmodel.yyear == 2020;
window2021 = mmodel.yyear == 2021;
window2022 = mmodel.yyear == 2022;
window2023 = mmodel.yyear == 2023;
window2024 = mmodel.yyear == 2024;
count2020 = sum(mmodel.yyear == 2020);
count2021 = sum(mmodel.yyear == 2021);
count2022 = sum(mmodel.yyear == 2022);
count2023 = sum(mmodel.yyear == 2023);
count2024 = sum(mmodel.yyear == 2024);

ce = @(vvv) (1 + (1 - param.sigma) * (1 - param.bet) * vvv).^(1/(1 - param.sigma));
ceLoss = @(cOld, cNew) (cNew - cOld) / (1 - param.bet) / 52;
qFromSp = @(sss) 1 ./ ( sss / param.pay + 1 );

oout = NaN([16, 1]);

oout(1) = 100 * max(mmodel.muDpath);
oout(2) = sum(mmodel.muIpath > 1.0e-6) / 4;

oout(3) = (max(mmodel.Bpath) - mmodel.Bpath(param.before)) ./ 52 * 100;
oout(4) = max(mmodel.spPath);
oout(5) = sum(mmodel.Dpath > 0) / 4;

ySample = param.before:size(mmodel.Ypath, 1);
phiDpath = mmodel.Ypath ./ (1 - param.thetaY * mmodel.Lpath).^0.67 ./ mmodel.Zpath;

oout(6) = -sum(mmodel.Ypath(ySample) - mmodel.Ypath(param.before-1)) ./ 52 * 100;
oout(7) = -sum((1 - param.thetaY * mmodel.Lpath(ySample)).^0.67 - 1) ./ 52 * 100;
oout(8) = -sum(phiDpath - 1) ./ 52 * 100;
oout(9) = -sum(mmodel.Zpath - 1) ./ 52 * 100;

cstart=mmodel.Cpath(param.before-1); 
oout(10) = 100 * (sum(mmodel.Cpath(window2020)) / (cstart*count2020) - 1);
oout(11) = 100 * (sum(mmodel.Cpath(window2021)) / (cstart*count2021) - 1);
oout(12) = 100 * (sum(mmodel.Cpath(window2022)) / (cstart*count2022) - 1);
oout(13) = 100 * (sum(mmodel.Cpath(window2023)) / (cstart*count2023) - 1);
oout(14) = 100 * (sum(mmodel.Cpath(window2024)) / (cstart*count2024) - 1);

oout(15) = -100 * ceLoss(ce(mmodel.Vpath(param.before-1)), ce(mmodel.Vpath(param.before)));
spOnImpact = (1 + mmodel.spPath(param.before) / 100).^(1/52) - 1;
oout(16) = -(qFromSp(spOnImpact) - 1) * mmodel.Bpath(param.before) ./ 52 * 100;
end

