clear;
%
%
% Load all results from out files (fortran output) to matlab .mat files
%
%

myReadMatrix = @(fff) readmatrix(fff, 'FileType','text', 'Delimiter','\t');

baselinePath = '../baseline/results/'; % [baselinePath, '
perfectPath = '../perfect_financial/results/'; % [perfectPath, '
persistentPath = '../persistent_recession/results/'; % [persistentPath, '
dssiPath = '../debt_suspension/results/'; % [dssiPath, '
fixedLPath = '../fixed_L/results/'; % [fixedLPath, '

ddata = myReadMatrix([baselinePath, 'parameters.tab']);
ix = 1;
param.muSsz = ddata(ix); ix = ix + 1;
param.muIsz = ddata(ix); ix = ix + 1;
param.bSz = ddata(ix); ix = ix + 1;
param.Tsz = ddata(ix); ix = ix + 1;
param.H = ddata(ix); ix = ix + 1;
param.before = ddata(ix); ix = ix + 1;
param.after = ddata(ix); ix = ix + 1;
param.chi = ddata(ix); ix = ix + 1;
param.sigma = ddata(ix); ix = ix + 1;
param.bet = ddata(ix); ix = ix + 1;
param.rf = ddata(ix); ix = ix + 1;
param.duration = ddata(ix); ix = ix + 1;
param.delta = ddata(ix); ix = ix + 1;
param.pay = ddata(ix); ix = ix + 1;
param.kappa = ddata(ix); ix = ix + 1;
param.gamm0 = ddata(ix); ix = ix + 1;
param.gamm1 = ddata(ix); ix = ix + 1;
param.gamm2 = ddata(ix); ix = ix + 1;
param.pp = ddata(ix); ix = ix + 1;
param.rho_psi = ddata(ix); ix = ix + 1;
param.psi_wave1_end = ddata(ix); ix = ix + 1;
param.psi_wave2_end = ddata(ix); ix = ix + 1;
param.piD = ddata(ix); ix = ix + 1;
param.piDsq = ddata(ix); ix = ix + 1;
param.theta = ddata(ix); ix = ix + 1;
param.thetaY = ddata(ix);

muSsz = param.muSsz;
muIsz = param.muIsz;
bSz = param.bSz;

param.muSgrid = myReadMatrix([baselinePath, 'muS.tab']);
param.muIgrid = myReadMatrix([baselinePath, 'muI.tab']);
param.bGrid = myReadMatrix([baselinePath, 'bGrid.tab']);

stationary.EC = reshape(myReadMatrix([baselinePath, 'reg_EC.tab']), [bSz, 1]);
stationary.Ed = reshape(myReadMatrix([baselinePath, 'reg_Ed.tab']), [bSz, 1]);
stationary.EbPr = reshape(myReadMatrix([baselinePath, 'reg_EbPr.tab']), [bSz, 1]);
stationary.V = reshape(myReadMatrix([baselinePath, 'reg_V.tab']), [bSz, 1]);
stationary.q = reshape(myReadMatrix([baselinePath, 'reg_q.tab']), [bSz, 1]);
stationary.polPr = reshape(myReadMatrix([baselinePath, 'reg_polPr.tab']), [bSz, bSz]);
stationary.C = reshape(myReadMatrix([baselinePath, 'reg_C.tab']), [bSz, bSz]);
stationary.d = reshape(myReadMatrix([baselinePath, 'reg_d.tab']), [bSz, bSz]);
stationary.dStar = reshape(myReadMatrix([baselinePath, 'reg_dStar.tab']), [bSz, bSz]);
stationary.spBond = param.pay * (1 ./ stationary.q - 1);

%
% Data
%
ddata = myReadMatrix([baselinePath, 'data.tab']);
ix = 1;
theData.yyear = ddata(:, ix); ix = ix + 1;
theData.mmonth = ddata(:, ix); ix = ix + 1;
theData.dday = ddata(:, ix); ix = ix + 1;
theData.yyww = datetime(theData.yyear, theData.mmonth, theData.dday);
theData.Lpath = ddata(:, ix); ix = ix + 1;
theData.muDpath = ddata(:, ix); ix = ix + 1;
theData.wedgePath = ddata(:, ix); 

dataMoments = readmatrix('../data/Output.xls', ...
  'Sheet', 'Table 1', ...
  'Range', 'C2:C13');

dataMasks = readmatrix('../data/data_mask.txt');
dataDates = readmatrix('../data/data_dates.txt');

%
% Baseline
%
ddata = myReadMatrix([baselinePath, 'sir_baseline.tab']);
model_baseline = modelFromMatrix(ddata);

ddata = myReadMatrix([baselinePath, 'sir_baseline_50.tab']);
model_baseline_50 = modelFromMatrix(ddata);

ddata = myReadMatrix([baselinePath, 'sir_baseline_70.tab']);
model_baseline_70 = modelFromMatrix(ddata);

%Loans
ddata = myReadMatrix([baselinePath, 'sir_baseline_loan.tab']);
model_baseline_loan = modelFromMatrix(ddata);

ddata = myReadMatrix([baselinePath, 'sir_baseline_50_loan.tab']);
model_baseline_50_loan = modelFromMatrix(ddata);

ddata = myReadMatrix([baselinePath, 'sir_baseline_70_loan.tab']);
model_baseline_70_loan = modelFromMatrix(ddata);

% Voluntary restructuring
ddata = myReadMatrix([baselinePath, 'sir_baseline_vol.tab']);
model_baseline_vol = modelFromMatrix(ddata);

% Sensitivity
ddata = myReadMatrix([baselinePath, 'sir_baseline_chi.tab']);
model_baseline_chi = modelFromMatrix(ddata);

ddata = myReadMatrix([baselinePath, 'sir_baseline_piD1.tab']);
model_baseline_piD1 = modelFromMatrix(ddata);

ddata = myReadMatrix([baselinePath, 'sir_baseline_H2.tab']);
model_baseline_H2 = modelFromMatrix(ddata);

ddata = myReadMatrix([baselinePath, 'sir_baseline_H4.tab']);
model_baseline_H4 = modelFromMatrix(ddata);

ddata = myReadMatrix([baselinePath, 'sir_baseline_2wave.tab']);
model_baseline_2wave = modelFromMatrix(ddata);

%
% Persistent recession
%
ddata = myReadMatrix([persistentPath, 'sir_persistent.tab']);
model_persistent = modelFromMatrix(ddata);

%
% Debt Suspension
%
ddata = myReadMatrix([dssiPath, 'sir_dssi.tab']);
model_dssi = modelFromMatrix(ddata);

ddata = myReadMatrix([dssiPath, 'sir_dssi_k1.tab']);
model_dssi_kappa1 = modelFromMatrix(ddata);

% Fixed L (appendix)
ddata = myReadMatrix([fixedLPath, 'sir_fixedL.tab']);
model_baseline_fixedL = modelFromMatrix(ddata);

%
% Perfect Financial Markets
%
ddata = readmatrix([perfectPath, 'perfect_param.tab'], 'FileType', 'text', 'Delimiter', '\t');
ix = 1;
paramPerfect.initS = ddata(ix); ix = ix + 1;
paramPerfect.initI = ddata(ix); ix = ix + 1;
paramPerfect.initR = ddata(ix); ix = ix + 1;
paramPerfect.pp = ddata(ix); ix = ix + 1;
paramPerfect.piD0 = ddata(ix); ix = ix + 1;
paramPerfect.piD1 = ddata(ix); ix = ix + 1;
paramPerfect.theta = ddata(ix); ix = ix + 1;
paramPerfect.thetaY = ddata(ix); ix = ix + 1;
paramPerfect.rfRate = ddata(ix); ix = ix + 1;
paramPerfect.beta = ddata(ix); ix = ix + 1;
paramPerfect.chi = ddata(ix); ix = ix + 1;
paramPerfect.crra = ddata(ix); ix = ix + 1;
paramPerfect.drs = ddata(ix); ix = ix + 1;
paramPerfect.days = ddata(ix); ix = ix + 1;
paramPerfect.maxT = ddata(ix); ix = ix + 1;
paramPerfect.freeVaxAt = ddata(ix); ix = ix + 1;
paramPerfect.waveTwoAt = ddata(ix); ix = ix + 1;
paramPerfect.rho_psi = ddata(ix); ix = ix + 1;
paramPerfect.psi_wave1_end = ddata(ix); ix = ix + 1;
paramPerfect.psi_wave2_end = ddata(ix); ix = ix + 1;
paramPerfect.welfare = ddata(ix); ix = ix + 1;
paramPerfect.welfareNoCOVID = ddata(ix); ix = ix + 1;
paramPerfect.ce = ddata(ix); ix = ix + 1;
paramPerfect.ceNoCOVID = ddata(ix); ix = ix + 1;
paramPerfect.outputLoss = ddata(ix); ix = ix + 1;
paramPerfect.totalDeaths = ddata(ix); ix = ix + 1;
paramPerfect.lossValL = ddata(ix); ix = ix + 1;
paramPerfect.lossValD = ddata(ix);

ddata = readmatrix([perfectPath, 'perfect_sim.tab'], 'FileType', 'text', 'Delimiter', '\t');
model_perfect = modelPerfectFromMatrix(ddata);
clear ddata;

%
% Baseline Policies
%
policies.V001 = loadBinary([baselinePath, 'sir_V_001.bin'], [muSsz, muIsz, bSz]);
policies.q001 = loadBinary([baselinePath, 'sir_q_001.bin'], [muSsz, muIsz, bSz]);
policies.B001 = loadBinary([baselinePath, 'sirEb_001.bin'], [muSsz, muIsz, bSz]);
policies.L001 = loadBinary([baselinePath, 'sirEL_001.bin'], [muSsz, muIsz, bSz]);
policies.C001 = loadBinary([baselinePath, 'sirEC_001.bin'], [muSsz, muIsz, bSz]);
policies.d001 = loadBinary([baselinePath, 'sirEd_001.bin'], [muSsz, muIsz, bSz]);

clear ix myReadMatrix *Path;
save inMatlab.mat;