function data = selectEpochMeasurements(rxTime, raw)
% SELECTEPOCHMEASUREMENTS 取出一个接收历元的测量，每颗卫星保留一个信号：
% 优先最小 signalId，同信号取最高 C/N0。返回 n×14 矩阵，列定义：
% 1:3 卫星位置，4:6 卫星速度，7 卫星钟差修正后的伪距，8 C/N0，9:11 Doppler/变化率/相位，
% 12:14 星座(0 GPS/1 Galileo)、PRN、signalId。
data = zeros(0,14);
if isempty(raw), return; end
rows = raw(abs([raw.rxTime]-rxTime)<1e-6);
if isempty(rows), return; end
signal = [rows.signalId]';
system = double(signal>=2);
candidates = [reshape([rows.satPosition],3,[])',reshape([rows.satVelocity],3,[])', ...
    [rows.prRaw]'+[rows.satClockM]',[rows.cn0]',[rows.fdRaw]', ...
    [rows.fdRateRaw]',[rows.carrierPhase]',system,[rows.svId]',signal];
valid = all(isfinite(candidates(:,[1:3,7])),2) & candidates(:,7)>0 & ismember(signal,0:4);
candidates = candidates(valid,:);
if isempty(candidates), return; end
cn0 = candidates(:,8); cn0(~isfinite(cn0)) = -Inf;
[~,order] = sortrows([candidates(:,12:14),-cn0],[1 2 3 4]);
candidates = candidates(order,:);
[~,keep] = unique(candidates(:,12:13),'rows','stable');
data = candidates(keep,:);
end
