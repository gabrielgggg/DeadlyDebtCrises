function [ mmodel ] = modelPerfectFromMatrix(ddata)
ix = 1;
mmodel.timeIx = ddata(:, ix); ix = ix + 1;
mmodel.R0 = ddata(:, ix); ix = ix + 1;
mmodel.wedge = ddata(:, ix); ix = ix + 1;
mmodel.muS = ddata(:, ix); ix = ix + 1;
mmodel.muI = ddata(:, ix); ix = ix + 1;
mmodel.muD = ddata(:, ix); ix = ix + 1;
mmodel.muR = 1.0 - mmodel.muS - mmodel.muI - mmodel.muD;
mmodel.Y = ddata(:, ix); ix = ix + 1;
mmodel.C = ddata(:, ix); ix = ix + 1;
mmodel.cfC = ddata(:, ix); ix = ix + 1;
mmodel.cfB = ddata(:, ix); ix = ix + 1;
mmodel.L = ddata(:, ix); ix = ix + 1;
mmodel.B = ddata(:, ix); ix = ix + 1;
mmodel.pushCash = ddata(:, ix); ix = ix + 1;

mmodel.dataMuD = ddata(:, ix); ix = ix + 1;
mmodel.dataL = ddata(:, ix); ix = ix + 1;

mmodel.dataMuD(mmodel.dataMuD < -1) = NaN;
mmodel.dataL(mmodel.dataL < -1) = NaN;

mmodel.beta = ddata(:, ix); ix = ix + 1;
mmodel.dates = datetime(ddata(:, ix), 'ConvertFrom', 'epochtime', 'TicksPerSecond', 1);
end

