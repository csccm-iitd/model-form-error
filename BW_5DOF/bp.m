function [bp_data] = bp(data, fsample, lfc, hfc)
fnyq = fsample/2;
wl = lfc/fnyq;
wh = hfc/fnyq;
P = 2*ceil(8*4/(wh-wl));
Coeff = fir1(P,[wl, wh],'noscale');
bp_data = filter(Coeff, 1, data);
end