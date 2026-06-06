function [phiID] = PhiID_BOLD(my_BOLD, varargin)
%PHIID_BOLD MATLAB R2018£©

% ============ inputParser arguments ==========
p = inputParser;
addRequired(p, 'my_BOLD');
addParameter(p, 'infodefi', 'mmi', @(x) any(validatestring(x, {'mmi', 'ccs'})));
addParameter(p, 'datatype', 'continuous', @(x) any(validatestring(x, {'continuous', 'discrete'})));
parse(p, my_BOLD, varargin{:});

infodefi = p.Results.infodefi;
datatype = p.Results.datatype;
% ==============================================

key = string(infodefi) + "_" + string(datatype);

switch key
    case "mmi_continuous"
        disp('Branch 1: infodefi = mmi, datatype = continuous');
        phiFunc = @PhiIDFull;
        measure = 'mmi';

    case "mmi_discrete"
        disp('Branch 2: infodefi = mmi, datatype = discrete');
        phiFunc = @PhiIDFullDiscrete;
        measure = 'mmi';

    case "ccs_continuous"
        disp('Branch 3: infodefi = ccs, datatype = continuous');
        phiFunc = @PhiIDFull;
        measure = 'ccs';

    case "ccs_discrete"
        disp('Branch 4: infodefi = ccs, datatype = discrete');
        phiFunc = @PhiIDFullDiscrete;
        measure = 'ccs';

    otherwise
        error('Unknown key: %s', key);
end

tic;

n = size(my_BOLD,1);

rtr = zeros(n); rtx = zeros(n); rty = zeros(n); rts = zeros(n);
xtr = zeros(n); xtx = zeros(n); xty = zeros(n); xts = zeros(n);
ytr = zeros(n); ytx = zeros(n); yty = zeros(n); yts = zeros(n);
str = zeros(n); stx = zeros(n); sty = zeros(n); sts = zeros(n);

% smooth
my_BOLD = smoothdata(my_BOLD, 2, "movmean", 8);

% filter
fs = 1; % Sampling frequency, measured in Hz
fpass = [0.008, 0.09]; % Passband frequency range of the bandpass filter, measured in Hz
order = 5; % Order of the filter
% Set up the bandpass filter
[b, a] = butter(order, fpass / (fs / 2), 'bandpass');
for i = 1:n
    my_BOLD(i,:) = filtfilt(b, a, my_BOLD(i,:));
end

for i = 1:n-1
    for j = i+1:n
        
        atoms = phiFunc([my_BOLD(i,:); my_BOLD(j,:)], 1, measure);

        rtr(i,j) = atoms.rtr; rtx(i,j) = atoms.rtx; rty(i,j) = atoms.rty; rts(i,j) = atoms.rts;
        xtr(i,j) = atoms.xtr; xtx(i,j) = atoms.xtx; xty(i,j) = atoms.xty; xts(i,j) = atoms.xts;
        ytr(i,j) = atoms.ytr; ytx(i,j) = atoms.ytx; yty(i,j) = atoms.yty; yts(i,j) = atoms.yts;
        str(i,j) = atoms.str; stx(i,j) = atoms.stx; sty(i,j) = atoms.sty; sts(i,j) = atoms.sts;

        rtr(j,i) = atoms.rtr; rtx(j,i) = atoms.rty; rty(j,i) = atoms.rtx; rts(j,i) = atoms.rts;
        xtr(j,i) = atoms.ytr; xtx(j,i) = atoms.yty; xty(j,i) = atoms.ytx; xts(j,i) = atoms.yts;
        ytr(j,i) = atoms.xtr; ytx(j,i) = atoms.xty; yty(j,i) = atoms.xtx; yts(j,i) = atoms.xts;
        str(j,i) = atoms.str; stx(j,i) = atoms.sty; sty(j,i) = atoms.stx; sts(j,i) = atoms.sts;
    end
end

phiID = struct('rtr', rtr, 'rtx', rtx, 'rty', rty, 'rts', rts, ...
              'xtr', xtr, 'xtx', xtx, 'xty', xty, 'xts', xts, ...
              'ytr', ytr, 'ytx', ytx, 'yty', yty, 'yts', yts, ...
              'str', str, 'stx', stx, 'sty', sty, 'sts', sts);

elapsedTime = toc;
fprintf('Total computation time: %.4f seconds.\n', elapsedTime);

end