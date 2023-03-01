function [ out ] = loadBinary(fname, sz)
    fileid = fopen(fname, 'r');
    out = fread(fileid, [prod(sz), 1], 'float64');
    out = reshape(out, sz);
    fclose(fileid);
end
