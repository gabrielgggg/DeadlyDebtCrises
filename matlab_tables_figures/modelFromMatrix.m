function [ mmodel ] = modelFromMatrix(ddata)
ix = 1;
mmodel.xxVals = ddata(:, ix); ix = ix + 1;
mmodel.muSpath = ddata(:, ix); ix = ix + 1;
mmodel.muIpath = ddata(:, ix); ix = ix + 1;
mmodel.muRpath = ddata(:, ix); ix = ix + 1;
mmodel.muDpath = ddata(:, ix); ix = ix + 1;
mmodel.Bpath = ddata(:, ix); ix = ix + 1;
mmodel.Cpath = ddata(:, ix); ix = ix + 1;
mmodel.Dpath = ddata(:, ix); ix = ix + 1;
mmodel.spPath = 100.0 * ddata(:, ix); ix = ix + 1;
mmodel.Ypath = ddata(:, ix); ix = ix + 1;
mmodel.Vpath = ddata(:, ix); ix = ix + 1;
mmodel.Zpath = ddata(:, ix); ix = ix + 1;
mmodel.Lpath = ddata(:, ix); ix = ix + 1;
mmodel.cashPath = ddata(:, ix); ix = ix + 1;
mmodel.wedgePath = ddata(:, ix); ix = ix + 1;
mmodel.yyear = ddata(:, ix); ix = ix + 1;
mmodel.mmonth = ddata(:, ix); ix = ix + 1;
mmodel.dday = ddata(:, ix);
mmodel.yyww = datetime(mmodel.yyear, mmodel.mmonth, mmodel.dday);
end

