function [phiID] = PhiID_BOLD(my_BOLD)
%PHIID_BOLD 此处显示有关此函数的摘要
%   此处显示详细说明

n = size(my_BOLD,1);

rtr = zeros(n);
rtx = zeros(n);
rty = zeros(n);
rts = zeros(n);
xtr = zeros(n);
xtx = zeros(n);
xty = zeros(n);
xts = zeros(n);
ytr = zeros(n);
ytx = zeros(n);
yty = zeros(n);
yts = zeros(n);
str = zeros(n);
stx = zeros(n);
sty = zeros(n);
sts = zeros(n);

for i = 1:n-1
    for j = i+1:n
        atoms = PhiIDFull([my_BOLD(i,:); my_BOLD(j, :)]);
        rtr(i,j) = atoms.rtr;
        rtx(i,j) = atoms.rtx;
        rty(i,j) = atoms.rty;
        rts(i,j) = atoms.rts;
        xtr(i,j) = atoms.xtr;
        xtx(i,j) = atoms.xtx;
        xty(i,j) = atoms.xty;
        xts(i,j) = atoms.xts;
        ytr(i,j) = atoms.ytr;
        ytx(i,j) = atoms.ytx;
        yty(i,j) = atoms.yty;
        yts(i,j) = atoms.yts;
        str(i,j) = atoms.str;
        stx(i,j) = atoms.stx;
        sty(i,j) = atoms.sty;
        sts(i,j) = atoms.sts;

        rtr(j,i) = atoms.rtr;
        rtx(j,i) = atoms.rty;
        rty(j,i) = atoms.rtx;
        rts(j,i) = atoms.rts;
        xtr(j,i) = atoms.ytr;
        xtx(j,i) = atoms.yty;
        xty(j,i) = atoms.ytx;
        xts(j,i) = atoms.yts;
        ytr(j,i) = atoms.xtr;
        ytx(j,i) = atoms.xty;
        yty(j,i) = atoms.xtx;
        yts(j,i) = atoms.xts;
        str(j,i) = atoms.str;
        stx(j,i) = atoms.sty;
        sty(j,i) = atoms.stx;
        sts(j,i) = atoms.sts;
    end
end

phiID = struct('rtr', rtr, 'rtx', rtx, 'rty', rty, 'rts', rts, ...
              'xtr', xtr, 'xtx', xtx, 'xty', xty, 'xts', xts, ...
              'ytr', ytr, 'ytx', ytx, 'yty', yty, 'yts', yts, ...
              'str', str, 'stx', stx, 'sty', sty, 'sts', sts);

end

