function satPos = getPos(rxTime, RAW)
% GETPOS 每颗卫星保留一个信号：优先最小 signalId，同信号取最高 C/N0。
% 1:3 位置，4:6 速度，7 修正伪距，8 C/N0，9:11 Doppler/变化率/相位。
% 12:14 星座(0 GPS/1 Galileo)、PRN、signalId，供诊断保留身份。
satPos = zeros(0,14);
if isempty(RAW), return; end
rows = RAW(abs([RAW.rxTime]-rxTime)<1e-6);
if isempty(rows), return; end
sig = [rows.signalId]';
system = double(sig>=2);
data = [reshape([rows.ECEF],3,[])',reshape([rows.Vel],3,[])', ...
    [rows.prRaw]'+[rows.clockBias]',[rows.cn0]',[rows.fdRaw]', ...
    [rows.fdRateRaw]',[rows.carrierPhase]',system,[rows.svId]',sig];
valid = all(isfinite(data(:,[1:3,7])),2) & data(:,7)>0 & ismember(sig,0:4);
data = data(valid,:);
if isempty(data), return; end
cn0 = data(:,8); cn0(~isfinite(cn0)) = -Inf;
[~,order] = sortrows([data(:,12:14),-cn0],[1 2 3 4]);
data = data(order,:);
[~,keep] = unique(data(:,12:13),'rows','stable');
satPos = data(keep,:);
end
