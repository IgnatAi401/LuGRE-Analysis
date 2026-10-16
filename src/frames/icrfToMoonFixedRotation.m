function rotation = icrfToMoonFixedRotation(gpsSeconds)
% ICRFTOMOONFIXEDROTATION ICRF -> 月固系（DE405 月球主轴系 PA）旋转矩阵；时间输入 GPST。
% 姿态来自 Aerospace Toolbox moonLibration（DE405，本机自带，无需下载）：
% 欧拉角 [phi theta psi] 为 3-1-3 顺序，R = Rz(psi)*Rx(theta)*Rz(phi)，把 ICRF 分量转到月固系分量。
% moonLibration 的时间为 TDB 儒略日；TDB-TT < 2 ms，取 TT = GPST + 19 s + 32.184 s。
% gpsSeconds 可为向量，返回 3x3xN。
if exist('moonLibration', 'file') ~= 2
    error('LuGRE:AerospaceRequired', '月球姿态需要 Aerospace Toolbox 的 moonLibration。');
end
gpsSeconds = gpsSeconds(:);
% 两列 [J0 d] 输入保持亚毫秒时间精度。
jdTt = [2444244.5*ones(numel(gpsSeconds),1), (gpsSeconds + 19 + 32.184)/86400];
angles = moonLibration(jdTt);
phi = angles(:,1); theta = angles(:,2); psi = angles(:,3);
n = numel(gpsSeconds);
rotation = zeros(3,3,n);
for k = 1:n
    rotation(:,:,k) = rotZ(psi(k))*rotX(theta(k))*rotZ(phi(k));
end
end

function r = rotZ(a)
r = [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];
end

function r = rotX(a)
r = [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
end
