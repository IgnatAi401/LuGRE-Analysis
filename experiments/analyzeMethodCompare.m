function analyzeMethodCompare(folder)
% ANALYZEMETHODCOMPARE 在 compareMethods 全部 OP 结果上做跨 OP 的诊断统计，回答“什么方法在什么条件下好、为什么”。
% 输出到 folder/analysis/：
%   opFeatures.csv     每个 OP 的观测条件（时长、距离、≥4 星占比与所在时段、平均星数）与批处理残差尺度
%   batchDrift.csv     转移段批处理误差（近似惯性系）的线性拟合 e(t)=e0+δv·(t−tc)：速度误差、拟合解释度、
%                      误差最小时刻，以及 ≥4 星历元与 1～3 星历元上的误差中位数
%   ekfGrowth.csv      EKF：启动时刻、多星时段内误差、最后一个 ≥4 星历元处与窗口末误差
%   bySatCount.csv     逐历元 LS 与 NAV 的误差中位数按当历元星数分档
%   fairCompare.csv    同一组历元（EKF 已启动 ∩ NAV 历元 ∩ 真值）上各方法误差中位数，用于和 EKF 公平比较
%   clockCompare.csv   dpeGrid 每历元钟差与 dpeIrls 样条钟差之差（去线性趋势后 RMS）及各自的短时起伏
%   epochPrior.csv     逐历元网格 DPE（先验）与轨道预测、逐历元 LS 的距离，网格边界解比例
% 只读结果文件，不重新定位。
root = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1, folder = fullfile(root,'results','method_compare_all'); end
outFolder = fullfile(folder,'analysis');
if ~exist(outFolder,'dir'), mkdir(outFolder); end
files = dir(fullfile(folder,'OP*_*.mat'));
files = files(~contains({files.name},'_epochL1'));
key = cellfun(@(n)sscanf(n,'OP%d_%d')'*[100;1], {files.name});
[~,order] = sort(key); files = files(order);

omegaEarth = 7.2921150e-5;
batchNames = {'batchDpeGrid','dpeGrid'; 'batchDpeIrls','dpeIrls'; 'batchLs','ls'};
ekfNames = {'ekfOrbit','ekfCvBlend','ekfLunarStatic'};
fairNames = {'batchDpeGrid','batchDpeIrls','batchLs','epochDpeOrbit','epochLs','ekfOrbit','ekfCvBlend', ...
    'ekfLunarStatic','nav'};
features = {}; drift = {}; growth = {}; bySat = {}; fair = {}; clocks = {}; prior = {};
for k = 1:numel(files)
    s = load(fullfile(folder,files(k).name),'series','info');
    info = s.info; se = s.series;
    if info.epochsWith4 == 0, continue; end
    t = se.gpsSeconds(:); n = numel(t);
    tRel = (t - t(1))/max(t(end) - t(1), eps);
    sat = se.satCount(:);
    has4 = sat >= 4;
    truthOk = all(isfinite(se.truth),2);
    lunar = isfield(info,'batchDynamics') && strcmp(info.batchDynamics,'lunarStatic');
    rs = nan(1,3); fs = nan(1,3);
    if isfield(info,'robustSigma')
        rs = [info.robustSigma.dpeGrid info.robustSigma.dpeIrls info.robustSigma.ls];
        fs = [info.formalStdMedian.dpeGrid info.formalStdMedian.dpeIrls info.formalStdMedian.ls];
    end
    features(end+1,:) = {info.window, lunar, (t(end)-t(1))/3600, info.meanRangeRe, n, info.epochsWith4, ...
        mean(has4), mean(sat(sat>0)), mean(sat==0), mean(sat>=1 & sat<=3), ...
        min(tRel(has4)), max(tRel(has4)), sum(sat), rs(1), rs(2), rs(3), fs(1), fs(2), fs(3)}; %#ok<AGROW>

    % 批处理误差的线性结构（转移段）。地心 30～60 地球半径处引力梯度约 μ/r³≈1e-10 s⁻²，一小时内
    % 轨道误差几乎按直线传播；ECEF 误差先绕 z 轴转回近似惯性系（只差常值旋转，不影响模长）。
    theta = omegaEarth*(t - t(1));
    for j = 1:size(batchNames,1)
        name = batchNames{j,1};
        if ~isfield(se,name), continue; end
        e = se.(name).positionEcef - se.truth;
        ok = all(isfinite(e),2);
        if nnz(ok) < 10, continue; end
        eI = [cos(theta).*e(:,1) - sin(theta).*e(:,2), sin(theta).*e(:,1) + cos(theta).*e(:,2), e(:,3)];
        tc = mean(t(ok));
        A = [ones(nnz(ok),1), t(ok) - tc];
        coef = A \ eI(ok,:);
        fitted = A*coef;
        explained = 1 - sum(vecnorm(eI(ok,:) - fitted,2,2).^2)/sum(vecnorm(eI(ok,:) - mean(eI(ok,:)),2,2).^2);
        d = se.(name).distanceM;
        [dMin, iMin] = min(d);
        converged = true;
        if isfield(info,'converged'), converged = info.converged.(batchNames{j,2}); end
        drift(end+1,:) = {info.window, name, converged, norm(coef(2,:)), norm(coef(1,:))/1000, explained, ...
            dMin/1000, tRel(iMin), median(d(ok & has4))/1000, median(d(ok & ~has4 & sat>0))/1000, ...
            median(d(ok & tRel>=0.25 & tRel<=0.75))/1000, median(d(ok & (tRel<0.1 | tRel>0.9)))/1000, ...
            rms(se.(name).radialM(ok))/1000, rms(se.(name).transverseM(ok))/1000}; %#ok<AGROW>
    end

    % EKF：只用过去观测，最后一个 ≥4 星历元之后靠动力学外推。
    for j = 1:numel(ekfNames)
        name = ekfNames{j};
        if ~isfield(se,name), continue; end
        d = se.(name).distanceM; ok = isfinite(d);
        if ~any(ok), continue; end
        first = find(ok,1); last4 = find(ok & has4,1,'last');
        okIdx = find(ok);
        % 多星时段（启动到最后一个 ≥4 星历元）内误差的线性漂移率，近似为速度误差
        span = ok & (1:n)' <= last4;
        e = se.(name).positionEcef - se.truth;
        eI = [cos(theta).*e(:,1) - sin(theta).*e(:,2), sin(theta).*e(:,1) + cos(theta).*e(:,2), e(:,3)];
        coef = [ones(nnz(span),1), t(span) - mean(t(span))] \ eI(span,:);
        growth(end+1,:) = {info.window, name, tRel(first), median(d(ok & has4))/1000, ...
            median(d(ok & ~has4))/1000, d(last4)/1000, d(okIdx(end))/1000, tRel(last4), ...
            rms(se.(name).radialM(ok))/1000, rms(se.(name).transverseM(ok))/1000, norm(coef(2,:))}; %#ok<AGROW>
    end

    % 钟差模型：每历元中位数钟差（dpeGrid）与 120 s 样条（dpeIrls）
    if isfield(info,'clockM')
        cg = info.clockM.dpeGrid(:); ci = info.clockM.dpeIrls(:);
        ok = isfinite(cg) & isfinite(ci);
        if nnz(ok) > 30
            A = [ones(nnz(ok),1), t(ok) - mean(t(ok))];
            dc = cg(ok) - ci(ok); dc = dc - A*(A\dc);
            clocks(end+1,:) = {info.window, nnz(ok), rms(dc), ...
                std(cg(ok) - movmedian(cg(ok),121,'SamplePoints',t(ok))), ...
                std(ci(ok) - movmedian(ci(ok),121,'SamplePoints',t(ok))), ...
                median(se.batchDpeIrls.distanceM(se.common))/median(se.batchDpeGrid.distanceM(se.common))}; %#ok<AGROW>
        end
    end

    % 逐历元网格 DPE（先验）是否被轨道预测主导
    if all(isfield(se,{'epochDpeOrbit','orbitPrediction','epochLs'}))
        ok = se.common & all(isfinite(se.epochDpeOrbit.positionEcef),2);
        boundary = NaN; if isfield(info,'epochDpeBoundary'), boundary = info.epochDpeBoundary; end
        prior(end+1,:) = {info.window, nnz(ok), ...
            median(vecnorm(se.epochDpeOrbit.positionEcef(ok,:) - se.orbitPrediction.positionEcef(ok,:),2,2))/1000, ...
            median(vecnorm(se.epochDpeOrbit.positionEcef(ok,:) - se.epochLs.positionEcef(ok,:),2,2))/1000, ...
            median(se.epochDpeOrbit.distanceM(ok))/1000, median(se.orbitPrediction.distanceM(ok))/1000, ...
            median(se.epochLs.distanceM(ok))/1000, boundary}; %#ok<AGROW>
    end

    % 逐历元 LS / NAV 误差按星数分档
    for name = {'epochLs','nav','epochDpeOrbit'}
        if ~isfield(se,name{1}), continue; end
        d = se.(name{1}).distanceM;
        for bin = 4:8
            if bin < 8, sel = sat == bin; else, sel = sat >= 8; end
            sel = sel & isfinite(d);
            if nnz(sel) < 5, continue; end
            bySat(end+1,:) = {info.window, name{1}, bin, nnz(sel), median(d(sel))/1000, ...
                median(abs(se.(name{1}).radialM(sel)))/1000, median(se.(name{1}).transverseM(sel))/1000}; %#ok<AGROW>
        end
    end

    % 公平比较：所有有结果的方法共同有值的 NAV 历元（EKF 存在时即 EKF 启动后的历元）
    present = fairNames(isfield(se,fairNames));
    both = se.common(:);
    for name = present, both = both & isfinite(se.(name{1}).distanceM); end
    for name = present
        d = se.(name{1}).distanceM(both);
        fair(end+1,:) = {info.window, name{1}, nnz(both), median(d)/1000, prctile(d,95)/1000}; %#ok<AGROW>
    end
end

writetable(cell2table(features,'VariableNames',{'window','lunar','durationH','rangeRe','rawEpochs', ...
    'epochsWith4','frac4','meanSatNonzero','frac0','frac1to3','first4Rel','last4Rel','observations', ...
    'robustSigmaGridM','robustSigmaIrlsM','robustSigmaLsM','formalStdGridM','formalStdIrlsM','formalStdLsM'}), ...
    fullfile(outFolder,'opFeatures.csv'));
writetable(cell2table(drift,'VariableNames',{'window','method','converged','velocityErrorMps','midErrorKm', ...
    'linearExplained','minErrorKm','minErrorRel','median4Km','median1to3Km','medianCenterKm','medianEdgeKm', ...
    'radialRmsKm','transverseRmsKm'}), fullfile(outFolder,'batchDrift.csv'));
writetable(cell2table(growth,'VariableNames',{'window','method','startRel','median4Km','median1to3Km', ...
    'atLast4Km','atEndKm','last4Rel','radialRmsKm','transverseRmsKm','velocityErrorMps'}), ...
    fullfile(outFolder,'ekfGrowth.csv'));
writetable(cell2table(clocks,'VariableNames',{'window','epochs','gridMinusIrlsRmsM','gridRoughnessM', ...
    'irlsRoughnessM','irlsOverGridError'}), fullfile(outFolder,'clockCompare.csv'));
writetable(cell2table(prior,'VariableNames',{'window','epochs','dpeToPredictionKm','dpeToLsKm', ...
    'dpeErrorKm','predictionErrorKm','lsErrorKm','boundaryFraction'}), fullfile(outFolder,'epochPrior.csv'));
writetable(cell2table(bySat,'VariableNames',{'window','method','satCount','epochs','medianKm', ...
    'medianAbsRadialKm','medianTransverseKm'}), fullfile(outFolder,'bySatCount.csv'));
writetable(cell2table(fair,'VariableNames',{'window','method','epochs','medianKm','p95Km'}), ...
    fullfile(outFolder,'fairCompare.csv'));
fprintf('分析结果已写入 %s\n', outFolder);
end
