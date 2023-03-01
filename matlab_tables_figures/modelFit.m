function [ ffit ] = modelFit(mmodel, param, lastDay)
window2020 = mmodel.yyear == 2020;
window2021 = mmodel.yyear == 2021;
count2020 = sum(mmodel.yyear == 2020);

ffit = NaN([12, 1]);
ffit(1) = 100 * mean(param.piD + param.piDsq * mmodel.muIpath(param.before:(param.before+lastDay-1)));

ffit(2) = 100 * mean(mmodel.muDpath(window2020));
ffit(3) = 100 * mean(mmodel.muDpath(window2021));

ffit(4) = max(mmodel.Lpath(window2020));
ffit(5) = mean(mmodel.Lpath(window2020));

ffit(6) = max(mmodel.Lpath(window2021));
ffit(7) = mean(mmodel.Lpath(window2021));

ffit(8) = mmodel.Bpath(param.before-1) ./ 52 * 100;
ffit(9) = max(mmodel.Bpath(window2020)) ./ 52 * 100 - ffit(8);

cstart=mmodel.Cpath(param.before-1); 
ffit(10) = 100 * (sum(mmodel.Cpath(mmodel.yyear == 2020)) / (cstart*count2020) - 1);

ffit(11) = max(mmodel.spPath);
ffit(12) = mean(mmodel.spPath(mmodel.yyear == 2020));
end

