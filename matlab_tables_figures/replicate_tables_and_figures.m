%%
%
% Load results
%
clear; clc;
load inMatlab.mat;

normalizedWinSize = [0.1 0.3 0.6 0.6];
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesFontSize', 26);
% get(groot,'factory')
% get(groot,'default')

%%
%
% Figure 1: A path
%
figure('Name', 'Figure 1: The Time Path of A_t', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(theData.yyww, theData.wedgePath ./ param.pp, '-');
hold on;
selectWindow = (param.before+size(theData.wedgePath, 1)):(param.before + param.H);
plot(model_baseline.yyww(selectWindow), model_baseline.wedgePath(selectWindow) ./ param.pp, '--');
ylabel('A_t / (1 - \pi_I)');
xlim([ datetime(2020, 3, 1), datetime(2022, 12, 31) ]);

%%
%
% Figure 2: baseline model time paths
%
timeLimits = [model_baseline.yyww(param.before-8), model_baseline.yyww(end-6*52)];

figure('Name', 'Figure 2: Fatalities', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(model_baseline.yyww, 100 * model_baseline.muDpath);
xlim(timeLimits);

figure('Name', 'Figure 2: Social Distancing', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(model_baseline.yyww, model_baseline.Lpath);
xlim(timeLimits);

figure('Name', 'Figure 2: Infected', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(model_baseline.yyww, model_baseline.muIpath);
xlim(timeLimits);

figure('Name', 'Figure 2: Susceptible', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(model_baseline.yyww((param.before-8):(param.before+param.H-2)), model_baseline.muSpath((param.before-8):(param.before+param.H-2))); 
xline(model_baseline.yyww(param.before+param.H-1));
xlim(timeLimits);

figure('Name', 'Figure 2: Debt', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(model_baseline.yyww, model_baseline.Bpath ./ 52);
xlim(timeLimits);

figure('Name', 'Figure 2: Spread', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(model_baseline.yyww, model_baseline.spPath);
xlim(timeLimits);

figure('Name', 'Figure 2: Output', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(model_baseline.yyww, model_baseline.Ypath ./ model_baseline.Ypath(param.before-1));
xlim(timeLimits);

figure('Name', 'Figure 2: Consumption', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(model_baseline.yyww, model_baseline.Cpath ./ model_baseline.Cpath(param.before-1));
xlim(timeLimits);

%%
%
% Table 2: Model fit
%
lastDay = size(theData.muDpath, 1);
theCaption = 'Table 2: Model Fit';
stats = ["Case fatality rate average", "Cumulative fatalities 2020 (per 100)", ...
  "Cumulative fatalities 2021 (per 100)", "Social Distancing: Peak 2020", ...
  "Social Distancing: Average 2020", "Social Distancing: Peak 2021", ...
  "Social Distancing: Average 2021", "Govt debt to output 2019", ...
  "Govt debt increase 2020", "Consumption decline 2020", "Spread increase peak 2020", ...
  "Spread increase average 2020" ]';
theTable = table;
theTable.moments = stats;
theTable.dataValues = round(dataMoments, 2);
theTable.baseline = round(modelFit(model_baseline, param, lastDay), 2);
theTable.Properties.VariableNames = ["Moments", "Data", "Baseline"];
disp(theCaption);
disp(theTable);
disp(' ');
uif = uifigure('Name', theCaption);
uitable(uif, 'ColumnWidth', 'auto', 'units', 'normalized', 'Position', [0 0 1 1], 'Data', theTable);

%%
%
% Table 3: Epidemic outcomes and financial markets
%
theCaption = 'Table 3: Epidemic Outcomes and Financial Markets';
stats = [ "Fatalities (per 100)", "Health crisis length (months)", ...
  "Peak debt increase", "Peak spread", "Debt crisis length (months)", ...
  "Cumulative output loss", "...from social distancing", "...from default cost", ...
  "...from low productivity", "C 2020", "C 2021", "C 2022", "C 2023", "C 2024", ...
  "Welfare loss country", "Welfare loss lenders" ]';
theTable = table;
theTable.moments = stats;
theTable.baseline = round(modelOutcomes(model_baseline, param), 2);
theTable.perfect = round(modelPerfectOutcomes(model_perfect, paramPerfect), 2);
theTable.persistent = round(modelOutcomes(model_persistent, param), 2);

theTable.Properties.VariableNames = ["Moments", "Baseline", "Perfect Financial", "Persistent Recession"];
disp(theCaption);
disp(theTable);
disp(' ');
uif = uifigure('Name', theCaption);
uitable(uif, 'ColumnWidth', 'auto', 'units', 'normalized', 'Position', [0 0 1 1], 'Data', theTable);

%%
%
% Figure 3: Epidemic dynamics
%
timeLimits = [model_baseline.yyww(param.before-8), model_baseline.yyww(param.before+5*52)];

figure('Name', 'Figure 3: Fatalities', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(model_baseline.yyww, 100 * model_baseline.muDpath, '-'); hold on;
plot(model_perfect.dates, 100 * model_perfect.muD, '--');
plot(model_persistent.yyww, 100 * model_persistent.muDpath, '-.');
xlim(timeLimits); legend("Baseline", "Perfect", "Persistent");

figure('Name', 'Figure 3: Debt', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(model_baseline.yyww, model_baseline.Bpath ./ 52, '-'); hold on;
plot(model_perfect.dates, (model_perfect.B - model_perfect.cfB + model_baseline.Bpath(param.before-1))  ./ 52, '--');
plot(model_persistent.yyww, model_persistent.Bpath ./ 52, '-.');
xlim(timeLimits); legend("Baseline", "Perfect", "Persistent");

figure('Name', 'Figure 3: Consumption', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(model_baseline.yyww, 1 + model_baseline.Cpath - model_baseline.Cpath(param.before-1), '-'); hold on;
plot(model_perfect.dates, 1 + model_perfect.C - model_perfect.cfC, '--');
plot(model_persistent.yyww, 1 + model_persistent.Cpath - model_persistent.Cpath(param.before-1), '-.');
xlim(timeLimits); ylim([0.825, 1.025]); yline(1); legend("Baseline", "Perfect", "Persistent");

%%
%
% Table 4: Debt relief counterfactuals
%
theCaption = 'Table 4: Debt Relief Counterfactuals';
stats = [ "Country welfare gain", "Debt crisis: length reduction", ...
  "Debt crisis: reduction in peak spread", "Health crisis: deaths prevented", ...
  "Consumption, avg % increase", "Lenders gain" ]';
theTable = table;
theTable.moments = stats;
theTable.loan50 = round(modelVersusModel(model_baseline_50_loan, model_baseline_50, param), 2);
theTable.loan = round(modelVersusModel(model_baseline_loan, model_baseline, param), 2);
theTable.loan70 = round(modelVersusModel(model_baseline_70_loan, model_baseline_70, param), 2);
theTable.dssik1 = round(modelVersusModel(model_dssi_kappa1, model_baseline, param), 2);
theTable.dssikk = round(modelVersusModel(model_dssi, model_baseline, param), 2);
theTable.voluntary = round(modelVersusModel(model_baseline_vol, model_baseline, param), 2);
theTable.voluntary(end) = 0.0;

theTable.Properties.VariableNames = [ "Moments", "Loan 51", "Loan 61", ...
  "Loan 71", "DSSI kD=1", "DSSI kD=k", "Voluntary" ];
disp(theCaption);
disp(theTable);
disp(' ');
uif = uifigure('Name', theCaption);
uitable(uif, 'ColumnWidth', 'auto', 'units', 'normalized', 'Position', [0 0 1 1], 'Data', theTable);

%%
%
% Figure 4: Scope for voluntary restructuring of the outstanding debt
%
ce = @(vvv) (1 + (1 - param.sigma) * (1 - param.bet) * vvv).^(1/(1 - param.sigma));

[nnmuSgrid, nnmuIgrid, nnbGrid] = ndgrid(param.muSgrid, param.muIgrid, param.bGrid);
iq = griddedInterpolant(nnmuSgrid, nnmuIgrid, nnbGrid, policies.q001);
iqStat = griddedInterpolant(param.bGrid, stationary.q);
iV = griddedInterpolant(nnmuSgrid, nnmuIgrid, nnbGrid, policies.V001);
ibPr = griddedInterpolant(nnmuSgrid, nnmuIgrid, nnbGrid, policies.B001);
iL = griddedInterpolant(nnmuSgrid, nnmuIgrid, nnbGrid, policies.L001);
id = griddedInterpolant(nnmuSgrid, nnmuIgrid, nnbGrid, policies.d001);

baseB = model_baseline.Bpath(param.before);

fineSz = 1000;
fineB = linspace(param.bGrid(1), param.bGrid(end), fineSz);

initI = model_baseline.muIpath(param.before);
initR = model_baseline.muRpath(param.before);
initS = 1.0 - initI - initR;
RzeroInit = model_baseline.wedgePath(param.before);

baselineBprime = ibPr( initS, initI,  baseB );
baselineL = iL( initS, initI,  baseB );
baselineD = id( initS, initI,  baseB );
baselineV = iV( initS, initI,  baseB );

varyV = iV( initS * ones([1, fineSz]), initI * ones([1, fineSz]),  fineB );
varyL = iL( initS * ones([1, fineSz]), initI * ones([1, fineSz]),  fineB );
varyD = id( initS * ones([1, fineSz]), initI * ones([1, fineSz]),  fineB );
varyBprime = ibPr( initS * ones([1, fineSz]), initI * ones([1, fineSz]),  fineB );

newInfect = RzeroInit * param.pp * (1 - param.theta * varyL).^2 * initS * initI;
varyIprime = (1 - param.pp) * initI + newInfect;
varySprime = initS - newInfect;
varyQ = iq(varySprime, varyIprime, varyBprime);
varyQschedule = iqStat(fineB);

baselineQ = interp1(fineB, varyQ, baseB);
varyQstat = interp1(param.bGrid, stationary.q, fineB);
varyDstat = interp1(param.bGrid, stationary.Ed, fineB);

valueToLenders = fineB .* ( ...
  varyQ .* ( 1 - param.delta + param.kappa * (param.delta + param.rf - 1) .* varyD ) + ...
  (param.delta + param.rf - 1) * (1 - varyD) );

valueToLendersSchedule = fineB .* varyQschedule;

baselineValueLenders = baseB .* ( ...
  baselineQ .* ( 1 - param.delta + param.kappa * (param.delta + param.rf - 1) .* baselineD ) + ...
  (param.delta + param.rf - 1) * (1 - baselineD) );

valueToLendersStat = fineB .* ( ...
  varyQstat .* ( 1 - param.delta + param.kappa * (param.delta + param.rf - 1) .* varyDstat ) + ...
  (param.delta + param.rf - 1) * (1 - varyDstat) );

figure('Name', 'Figure 4: Value to Country', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
plot(fineB ./ 52, (ce(varyV) - ce(baselineV)) ./ (1-param.bet) ./ 52 * 100 ); hold on;
yline(0, '--k');
scatter(baseB ./52, 0.0, 'd', 'filled');
xlim([0.55 0.625]);
xline(0.5635);
xlabel('Level of Debt');
ylabel('Gain to Country');

figure('Name', 'Figure 4: Value to Lenders', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);orient(gcf, 'landscape');
plot(fineB ./ 52, valueToLenders ./ 52.0); hold on;
scatter(baseB ./52, baselineValueLenders ./ 52.0, 'd', 'filled');
yline(baselineValueLenders ./ 52.0, '--k');
xlim([0.55 0.625]);
xline(0.5635);
xlabel('Level of Debt');
ylabel('Market Value of Debt');

%%
%
% Table 5: Sensitivity and alternative scenarios
%
theCaption = 'Table 3: Epidemic Outcomes and Financial Markets';
stats = [ "Fatalities (per 100)", "Health crisis length (months)", ...
  "Peak debt increase", "Peak spread", "Debt crisis length (months)", ...
  "Cumulative output loss", "...from social distancing", "...from default cost", ...
  "...from low productivity", "C 2020", "C 2021", "C 2022", "C 2023", "C 2024", ...
  "Welfare loss country", "Welfare loss lenders" ]';
theTable = table;
theTable.moments = stats;
theTable.baseline = round(modelOutcomes(model_baseline, param), 2);
theTable.chi = round(modelOutcomes(model_baseline_chi, param), 2);
theTable.piD1 = round(modelOutcomes(model_baseline_piD1, param), 2);
theTable.H2 = round(modelOutcomes(model_baseline_H2, param), 2);
theTable.H4 = round(modelOutcomes(model_baseline_H4, param), 2);
theTable.w2 = round(modelOutcomes(model_baseline_2wave, param), 2);

theTable.Properties.VariableNames = [ "Moments", "Baseline", ...
  "Chi", "piD1", "H2", "H4", "2wave" ];
disp(theCaption);
disp(theTable);
disp(' ');
uif = uifigure('Name', theCaption);
uitable(uif, 'ColumnWidth', 'auto', 'units', 'normalized', 'Position', [0 0 1 1], 'Data', theTable);
%%
%
% Figure A1: Exo component of disease transmission and masking
%
figure('Name', 'Figure 1: The Time Path of A_t', 'NumberTitle', 'off');
set(gcf,'units', 'normalized', 'position', normalizedWinSize);
yyaxis left;
plot(theData.yyww, theData.wedgePath ./ param.pp, '-');
ylabel('A_t');
yyaxis right;
plot(theData.yyww, dataMasks(8:end-6), '--');
ylabel('Masking');
xlim([ datetime(2020, 3, 1), datetime(2021, 11, 31) ]);

%%
%
% Table A1: Exogenous social distancing (Appendix)
%
theCaption = 'Table A1: Exogenous Social Distancing (Appendix)';
stats = [ "Fatalities (per 100)", "Health crisis length (months)", ...
  "Peak debt increase", "Peak spread", "Debt crisis length (months)", ...
  "Cumulative output loss", "...from social distancing", "...from default cost", ...
  "...from low productivity", "C 2020", "C 2021", "C 2022", "C 2023", "C 2024", ...
  "Welfare loss country", "Welfare loss lenders" ]';
theTable = table;
theTable.moments = stats;
theTable.baseline = round(modelOutcomes(model_baseline, param), 2);

tmp = modelPerfectOutcomes(model_perfect, paramPerfect);
theTable.perfect = round(tmp, 2);

theTable.fixedL = round(modelOutcomes(model_baseline_fixedL, param), 2);
theTable.perfectFixed = round(tmp, 2);

theTable.Properties.VariableNames = [ ...
  "Moments", "Endo L: Baseline", "Endo L: Perfect", ...
  "Exo L: Baseline", "Exo L: Perfect" ];
disp(theCaption);
disp(theTable);
disp(' ');
uif = uifigure('Name', theCaption);
uitable(uif, 'ColumnWidth', 'auto', 'units', 'normalized', 'Position', [0 0 1 1], 'Data', theTable);
