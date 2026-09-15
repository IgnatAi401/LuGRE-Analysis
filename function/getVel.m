function satVel = getVel(rxTime, RAW)
% GETPOS 根据 rxTime 匹配数据，输出卫星坐标和观测信息
%
% 输出 satPos (Nx8):
% [X Y Z 0 0 0 lambda cn0]

    c = 299792458; % 光速

    % 1. 提取所有时间并查找匹配索引
    allTimes = [RAW.rxTime];
    tol = 1e-12;
    idx = abs(allTimes - rxTime) < tol;

    if ~any(idx)
        satVel = [];
        return
    end

    %% --- ECEF Vel ---
    ecefTemp = [RAW(idx).ECEF];
    ecefMat = reshape(ecefTemp,3,[]).';
    velTemp = [RAW(idx).Vel];
    velMat = reshape(velTemp,3,[]).';

    %% --- signalId ---
    signalIdVec = [RAW(idx).signalId]';

    %% --- cn0 ---
    cn0vec = [RAW(idx).cn0]';

    %% --- 计算 lambda ---
    lambda = zeros(length(signalIdVec),1);

    for i = 1:length(signalIdVec)

        switch signalIdVec(i)
            case 0 % GPS L1CA
                f = 1575.42e6;

            case 1 % GPS L5
                f = 1176.45e6;

            case 2 % GAL E1BC
                f = 1575.42e6;

            case 3 % GAL E5A
                f = 1176.45e6;

            case 4 % GAL E5B
                f = 1207.14e6;

            otherwise
                f = NaN;
        end

        lambda(i) = c / f;

    end


    %% --- 拼接 ---
    satPos_Raw = [ecefMat, velMat, lambda, cn0vec];

    %% --- 去重 (按卫星位置) ---
    [~, idu] = unique(satPos_Raw(:,1:3),'rows','stable');
    satVel = satPos_Raw(idu,:);

end