function result = kalmanGnssPositioning(NAV, RAW_Pos, cfg)
%KALMANGNSSPOSITIONING 伪距 + Doppler 的 GNSS 扩展卡尔曼滤波定位。
%
% result = kalmanGnssPositioning(NAV, RAW_Pos)
% result = kalmanGnssPositioning(NAV, RAW_Pos, cfg)
%
% 每个历元通过下面的已有函数取得观测：
%   satPos = getPos(NAV(i).rxTime, RAW_Pos)
%
% satPos 各列：
%   1:3   卫星 ECEF 位置 (m)
%   4:6   卫星 ECEF 速度 (m/s)
%   7     修正后的伪距 (m)
%   8     信噪比/CN0 (dB 或 dB-Hz)
%   9     Doppler (Hz)
%   10    Doppler rate (Hz/s)，当前模型不使用
%
% 状态定义：
%   x = [x y z vx vy vz cb cd]'
%   cb 为接收机钟差乘以光速后的距离量 (m)
%   cd 为接收机钟漂乘以光速后的距离速率量 (m/s)
%
% 重要假设：
%   1. 第 7 列已完成卫星钟差等必要改正；
%   2. 第 1:6 列卫星状态与该观测对应，且处于同一 ECEF 坐标系；
%   3. 所有卫星使用同一时间系统（只有一个接收机钟差状态）；
%   4. 默认 Doppler 定义满足 rangeRate = -lambda * Doppler。

    if nargin < 3 || isempty(cfg)
        cfg = struct();
    end
    cfg = fillDefaultConfig(cfg);

    nEpoch = numel(NAV);
    if nEpoch == 0
        error('NAV 为空，无法进行滤波。');
    end

    result.state = nan(nEpoch, 8);
    result.posECEF = nan(nEpoch, 3);
    result.velECEF = nan(nEpoch, 3);
    result.clockBiasM = nan(nEpoch, 1);
    result.clockDriftMps = nan(nEpoch, 1);
    result.nPseudoUsed = zeros(nEpoch, 1);
    result.nDopplerUsed = zeros(nEpoch, 1);
    result.innovationRms = nan(nEpoch, 1);
    result.initialized = false(nEpoch, 1);
    result.config = cfg;

    x = zeros(8, 1);
    P = initialCovariance(cfg);
    isInitialized = false;

    for i = 1:nEpoch
        satPos = getPos(NAV(i).rxTime, RAW_Pos);
        validateSatPos(satPos, i);

        if ~isInitialized
            [ok, x0] = initializeState(satPos, cfg);
            if ~ok
                continue;
            end
            x = x0;
            P = initialCovariance(cfg);
            isInitialized = true;
        elseif i > 1
            [x, P] = predictState(x, P, cfg.dt, cfg);
        end

        [z, h, H, R, type] = buildMeasurements(x, satPos, cfg);

        if ~isempty(z)
            innovation = z - h;

            % 按单个观测的归一化新息剔除粗差。
            innovationVar = sum((H * P) .* H, 2) + diag(R);
            normalizedInnovation = abs(innovation) ./ sqrt(max(innovationVar, eps));
            keep = normalizedInnovation <= cfg.innovationGateSigma;

            z = z(keep);
            h = h(keep);
            H = H(keep, :);
            R = R(keep, keep);
            type = type(keep);

            if ~isempty(z)
                innovation = z - h;
                S = H * P * H' + R;
                K = (P * H') / S;
                x = x + K * innovation;

                % Joseph 形式比直接使用 (I-KH)P 更能保持协方差正定。
                I8 = eye(8);
                P = (I8 - K * H) * P * (I8 - K * H)' + K * R * K';
                P = (P + P') / 2;

                result.nPseudoUsed(i) = sum(type == 1);
                result.nDopplerUsed(i) = sum(type == 2);
                result.innovationRms(i) = sqrt(mean(innovation.^2));
            end
        end

        result.state(i, :) = x';
        result.posECEF(i, :) = x(1:3)';
        result.velECEF(i, :) = x(4:6)';
        result.clockBiasM(i) = x(7);
        result.clockDriftMps(i) = x(8);
        result.initialized(i) = true;
    end

    result.posLLH = ecefToLlhWgs84(result.posECEF);
end


function cfg = fillDefaultConfig(cfg)
    defaults.dt = 1.0;                    % 历元间隔 (s)
    defaults.c = 299792458.0;             % 光速 (m/s)
    defaults.carrierFrequencyHz = 1575.42e6; % GPS L1 (Hz)
    defaults.dopplerSign = -1.0;          % rangeRate = sign*lambda*Doppler

    % 连续白噪声强度。动态较剧烈时可增大 accelPSD。
    defaults.accelPSD = 0.5;              % m^2/s^3
    defaults.clockAccelPSD = 1.0;         % m^2/s^3
    defaults.clockBiasRandomWalkPSD = 0.1;% m^2/s

    % 观测标准差会再根据信噪比缩放。
    defaults.pseudorangeSigmaM = 5.0;
    defaults.rangeRateSigmaMps = 0.5;
    defaults.referenceSnrDb = 40.0;
    defaults.minPseudorangeSigmaM = 0.8;
    defaults.maxPseudorangeSigmaM = 50.0;
    defaults.minRangeRateSigmaMps = 0.03;
    defaults.maxRangeRateSigmaMps = 5.0;
    defaults.useSnrWeighting = true;

    defaults.innovationGateSigma = 6.0;
    defaults.initialPositionStdM = 50.0;
    defaults.initialVelocityStdMps = 5.0;
    defaults.initialClockBiasStdM = 100.0;
    defaults.initialClockDriftStdMps = 10.0;
    defaults.wlsMaxIterations = 20;
    defaults.wlsPositionToleranceM = 1e-3;

    names = fieldnames(defaults);
    for k = 1:numel(names)
        name = names{k};
        if ~isfield(cfg, name) || isempty(cfg.(name))
            cfg.(name) = defaults.(name);
        end
    end
end


function validateSatPos(satPos, epoch)
    if ~isnumeric(satPos) || ndims(satPos) ~= 2 || size(satPos, 2) < 10
        error('第 %d 历元的 satPos 必须是至少 10 列的数值矩阵。', epoch);
    end
end


function P = initialCovariance(cfg)
    std0 = [repmat(cfg.initialPositionStdM, 1, 3), ...
            repmat(cfg.initialVelocityStdMps, 1, 3), ...
            cfg.initialClockBiasStdM, cfg.initialClockDriftStdMps];
    P = diag(std0.^2);
end


function [x, P] = predictState(x, P, dt, cfg)
    F = eye(8);
    F(1:3, 4:6) = dt * eye(3);
    F(7, 8) = dt;

    q11 = dt^3 / 3;
    q12 = dt^2 / 2;
    q22 = dt;

    Q = zeros(8);
    Q(1:3, 1:3) = cfg.accelPSD * q11 * eye(3);
    Q(1:3, 4:6) = cfg.accelPSD * q12 * eye(3);
    Q(4:6, 1:3) = cfg.accelPSD * q12 * eye(3);
    Q(4:6, 4:6) = cfg.accelPSD * q22 * eye(3);

    Q(7, 7) = cfg.clockAccelPSD * q11 + ...
              cfg.clockBiasRandomWalkPSD * dt;
    Q(7, 8) = cfg.clockAccelPSD * q12;
    Q(8, 7) = cfg.clockAccelPSD * q12;
    Q(8, 8) = cfg.clockAccelPSD * q22;

    x = F * x;
    P = F * P * F' + Q;
    P = (P + P') / 2;
end


function [ok, x] = initializeState(satPos, cfg)
    validPr = all(isfinite(satPos(:, 1:3)), 2) & isfinite(satPos(:, 7));
    if sum(validPr) < 4
        ok = false;
        x = nan(8, 1);
        return;
    end

    satR = satPos(validPr, 1:3);
    pseudo = satPos(validPr, 7);
    snr = satPos(validPr, 8);
    sigmaPr = observationSigma(snr, cfg.pseudorangeSigmaM, ...
        cfg.minPseudorangeSigmaM, cfg.maxPseudorangeSigmaM, cfg);

    xb = zeros(4, 1); % [接收机位置(3); 钟差(m)]
    ok = true;

    for iter = 1:cfg.wlsMaxIterations
        delta = satR - xb(1:3)';
        range = sqrt(sum(delta.^2, 2));
        if any(range < 1)
            ok = false;
            break;
        end

        u = delta ./ range;
        A = [-u, ones(size(u, 1), 1)];
        residual = pseudo - (range + xb(4));
        W = diag(1 ./ sigmaPr.^2);
        normal = A' * W * A;

        if rcond(normal) < 1e-12
            ok = false;
            break;
        end

        correction = normal \ (A' * W * residual);
        xb = xb + correction;
        if norm(correction(1:3)) < cfg.wlsPositionToleranceM && ...
                abs(correction(4)) < cfg.wlsPositionToleranceM
            break;
        end
    end

    if ~ok || any(~isfinite(xb))
        ok = false;
        x = nan(8, 1);
        return;
    end

    x = zeros(8, 1);
    x(1:3) = xb(1:3);
    x(7) = xb(4);

    % 用首历元 Doppler 对速度和钟漂做加权最小二乘初始化。
    validDop = validPr & all(isfinite(satPos(:, 4:6)), 2) & ...
               isfinite(satPos(:, 9));
    if sum(validDop) >= 4
        satR = satPos(validDop, 1:3);
        satV = satPos(validDop, 4:6);
        snr = satPos(validDop, 8);
        delta = satR - x(1:3)';
        range = sqrt(sum(delta.^2, 2));
        u = delta ./ range;

        lambda = cfg.c / cfg.carrierFrequencyHz;
        measuredRate = cfg.dopplerSign * lambda * satPos(validDop, 9);
        y = measuredRate - sum(u .* satV, 2);
        A = [-u, ones(size(u, 1), 1)];
        sigmaRate = observationSigma(snr, cfg.rangeRateSigmaMps, ...
            cfg.minRangeRateSigmaMps, cfg.maxRangeRateSigmaMps, cfg);
        W = diag(1 ./ sigmaRate.^2);
        normal = A' * W * A;

        if rcond(normal) >= 1e-12
            velocityClockDrift = normal \ (A' * W * y);
            x(4:6) = velocityClockDrift(1:3);
            x(8) = velocityClockDrift(4);
        end
    end
end


function [z, h, H, R, type] = buildMeasurements(x, satPos, cfg)
    z = zeros(0, 1);
    h = zeros(0, 1);
    H = zeros(0, 8);
    variances = zeros(0, 1);
    type = zeros(0, 1); % 1=伪距，2=Doppler 伪距率

    lambda = cfg.c / cfg.carrierFrequencyHz;
    nSat = size(satPos, 1);

    for k = 1:nSat
        if ~all(isfinite(satPos(k, 1:3)))
            continue;
        end

        satR = satPos(k, 1:3)';
        delta = satR - x(1:3);
        geometricRange = norm(delta);
        if geometricRange < 1
            continue;
        end
        u = delta / geometricRange;
        snr = satPos(k, 8);

        if isfinite(satPos(k, 7))
            sigmaPr = observationSigma(snr, cfg.pseudorangeSigmaM, ...
                cfg.minPseudorangeSigmaM, cfg.maxPseudorangeSigmaM, cfg);
            row = zeros(1, 8);
            row(1:3) = -u';
            row(7) = 1;

            z(end + 1, 1) = satPos(k, 7); %#ok<AGROW>
            h(end + 1, 1) = geometricRange + x(7); %#ok<AGROW>
            H(end + 1, :) = row; %#ok<AGROW>
            variances(end + 1, 1) = sigmaPr^2; %#ok<AGROW>
            type(end + 1, 1) = 1; %#ok<AGROW>
        end

        if all(isfinite(satPos(k, 4:6))) && isfinite(satPos(k, 9))
            satV = satPos(k, 4:6)';
            relativeV = satV - x(4:6);
            predictedRate = u' * relativeV + x(8);
            measuredRate = cfg.dopplerSign * lambda * satPos(k, 9);

            % 对位置的导数也保留；远距离卫星时这一项很小，但模型更完整。
            dRateDPosition = -(relativeV' * (eye(3) - u * u')) / geometricRange;
            row = zeros(1, 8);
            row(1:3) = dRateDPosition;
            row(4:6) = -u';
            row(8) = 1;

            sigmaRate = observationSigma(snr, cfg.rangeRateSigmaMps, ...
                cfg.minRangeRateSigmaMps, cfg.maxRangeRateSigmaMps, cfg);
            z(end + 1, 1) = measuredRate; %#ok<AGROW>
            h(end + 1, 1) = predictedRate; %#ok<AGROW>
            H(end + 1, :) = row; %#ok<AGROW>
            variances(end + 1, 1) = sigmaRate^2; %#ok<AGROW>
            type(end + 1, 1) = 2; %#ok<AGROW>
        end
    end

    R = diag(variances);
end


function sigma = observationSigma(snr, baseSigma, minSigma, maxSigma, cfg)
    sigma = baseSigma * ones(size(snr));
    if ~cfg.useSnrWeighting
        return;
    end

    valid = isfinite(snr);
    scale = 10.^((cfg.referenceSnrDb - snr(valid)) / 20);
    sigma(valid) = min(max(baseSigma .* scale, minSigma), maxSigma);
end


function llh = ecefToLlhWgs84(ecef)
% 输出 [纬度(deg), 经度(deg), 椭球高(m)]，不依赖工具箱。
    llh = nan(size(ecef));
    a = 6378137.0;
    f = 1 / 298.257223563;
    e2 = f * (2 - f);

    for k = 1:size(ecef, 1)
        if ~all(isfinite(ecef(k, :)))
            continue;
        end

        x = ecef(k, 1);
        y = ecef(k, 2);
        z = ecef(k, 3);
        p = hypot(x, y);
        lon = atan2(y, x);

        if p < 1e-8
            lat = sign(z) * pi / 2;
            height = abs(z) - a * sqrt(1 - e2);
        else
            lat = atan2(z, p * (1 - e2));
            for iter = 1:10
                sinLat = sin(lat);
                N = a / sqrt(1 - e2 * sinLat^2);
                height = p / cos(lat) - N;
                newLat = atan2(z, p * (1 - e2 * N / (N + height)));
                if abs(newLat - lat) < 1e-13
                    lat = newLat;
                    break;
                end
                lat = newLat;
            end
            sinLat = sin(lat);
            N = a / sqrt(1 - e2 * sinLat^2);
            height = p / cos(lat) - N;
        end

        llh(k, :) = [rad2deg(lat), rad2deg(lon), height];
    end
end
