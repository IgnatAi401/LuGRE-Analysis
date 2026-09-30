function table = gridAblation(windows, outputFolder)
% GRIDABLATION dpeIrls 最弱方向一维网格（cfg.batch.weakSearch）的消融实验：网格开/关 × 不同初值偏差。
% 偏差加在默认初值（快照 LS 拟合）上，方向按参考历元的地心径向、横向定义。
% 判据：最终 L1 代价（网格关闭时代价更高即陷入局部极小）、与真值误差、收敛与耗时。
% 例：gridAblation({'OP5_0','OP9_0','OP14_0'})
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, genpath(fullfile(root,'src')));
if nargin < 1, windows = {'OP5_0','OP9_0','OP14_0'}; end
if nargin < 2, outputFolder = fullfile(root,'results','grid_ablation'); end
% 偏差：[径向位置 m, 横向位置 m, 径向速度 m/s]
perturbations = [0 0 0; 5e3 0 0; 2e4 0 0; 5e4 0 0; 2e5 0 0; 1e6 0 0; -5e4 0 0; ...
    0 5e4 0; 0 1e6 0; 0 0 10; 0 0 50; 0 0 200; 5e4 5e4 50];
gridSettings = [0 2]; % weakSearch.outerIterations：0 关闭，2 开启
rows = {};
for w = 1:numel(windows)
    cfg = lugreConfig(); cfg.data.window = windows{w};
    cfg.run.plot = false; cfg.batch.methods = {'dpeIrls'};
    fprintf('\n===== %s：读取数据与真值 =====\n', windows{w});
    base = runBatchPositioning(cfg);
    assert(isfield(base,'truth'),'LuGRE:Experiment','%s 没有真值，无法评估。',windows{w});
    epochs = base.gpsSeconds; m = numel(epochs);
    rotations = icrfToItrfRotation(epochs, cfg.frame);
    truth = base.truth.position;
    radialUnit = truth./vecnorm(truth,2,2);
    r0 = base.init.state(1:3); v0 = base.init.state(4:6);
    radial = r0/norm(r0);
    transverse = cross(radial, cross(v0, radial)); % 惯性速度在横向平面内的方向
    if norm(transverse) < eps, transverse = null(radial'); transverse = transverse(:,1); end
    transverse = transverse/norm(transverse);
    for p = 1:size(perturbations,1)
        d = perturbations(p,:);
        start = base.init.state + [d(1)*radial + d(2)*transverse; d(3)*radial];
        for g = gridSettings
            c = cfg; c.batch.weakSearch.outerIterations = g;
            c.batch.maxOuterIterations = 20; % 大偏差需要更多外迭代才能收敛，两组相同
            lastwarn('');
            timer = tic;
            try
                s = solveBatchIrls(base.observations, epochs, rotations, start, base.moon0, base.tRef, 'dpeIrls', c);
                seconds = toc(timer);
                e = s.positionEcef - truth; good = all(isfinite(e),2);
                radialError = sum(e.*radialUnit,2);
                rows(end+1,:) = {windows{w}, d(1)/1e3, d(2)/1e3, d(3), g>0, s.converged, s.iterations, ...
                    s.cost, median(vecnorm(e(good,:),2,2)), rms(radialError(good)), seconds, ...
                    ~isempty(lastwarn), s.onBoundary, s.state'}; %#ok<AGROW>
            catch ME
                rows(end+1,:) = {windows{w}, d(1)/1e3, d(2)/1e3, d(3), g>0, false, NaN, NaN, NaN, NaN, ...
                    toc(timer), true, false, nan(1,6)}; %#ok<AGROW>
                fprintf('  失败 %s\n', ME.message);
            end
            fprintf('  %-7s dR=%6.0f km dT=%6.0f km dV=%4.0f m/s grid=%d: conv=%d iter=%2d cost=%11.1f median=%9.0f m radialRMS=%9.0f m %.1f s\n', ...
                windows{w}, rows{end,2:7}, rows{end,8:11});
        end
    end
end
names = {'window','dRadialKm','dTransverseKm','dRadialVelocityMps','grid','converged','outerIterations', ...
    'cost','medianErrorM','radialRmsM','seconds','warned','gridOnBoundary','state'};
table = cell2table(rows,'VariableNames',names);
if ~isfolder(outputFolder), mkdir(outputFolder); end
save(fullfile(outputFolder,'grid_ablation.mat'),'table');
writetable(removevars(table,'state'), fullfile(outputFolder,'grid_ablation.csv'));
summarize(table);
end

function summarize(table)
% 同一初值下网格开/关的代价差与解差：代价差 > 0 表示关网格陷入更差的极小。
fprintf('\n%-7s %8s %8s %6s | %12s %12s | %10s %10s\n','窗口','dR/km','dT/km','dV', ...
    '代价差(关-开)','位置差/m','误差(关)/m','误差(开)/m');
keys = unique(table(:,{'window','dRadialKm','dTransverseKm','dRadialVelocityMps'}),'stable');
for k = 1:height(keys)
    match = strcmp(table.window,keys.window{k}) & table.dRadialKm==keys.dRadialKm(k) & ...
        table.dTransverseKm==keys.dTransverseKm(k) & table.dRadialVelocityMps==keys.dRadialVelocityMps(k);
    off = table(match & ~table.grid,:); on = table(match & table.grid,:);
    fprintf('%-7s %8.0f %8.0f %6.0f | %12.2f %12.1f | %10.0f %10.0f\n', keys.window{k}, ...
        keys.dRadialKm(k), keys.dTransverseKm(k), keys.dRadialVelocityMps(k), off.cost-on.cost, ...
        norm(off.state(1:3)-on.state(1:3)), off.medianErrorM, on.medianErrorM);
end
end
