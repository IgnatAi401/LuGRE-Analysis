function plotLugreBatch(results)
% PLOTLUGREBATCH 批处理结果：真值误差（含形式 3σ）、径向/横向分解、DPE 残差、网格代价。
time = results.gpsSeconds - results.gpsSeconds(1);
methods = results.config.batch.methods;
if isfield(results,'truth')
    figure('Name','批处理 与真值的三维误差'); hold on; grid on; set(gca,'YScale','log');
    names = fieldnames(results.truth);
    names = names(~ismember(names,{'position','valid'}));
    for k = 1:numel(names)
        plot(time,results.truth.(names{k}).distanceM,'DisplayName',strrep(names{k},'_','\_'));
    end
    for k = 1:numel(methods)
        plot(time,3*results.(methods{k}).formalStd3dM,'--', ...
            'DisplayName',['形式 3\sigma ' upper(methods{k})]);
    end
    xlabel('相对时间 / s'); ylabel('三维误差 / m'); legend('show','Location','best');
    figure('Name','批处理 径向/横向误差'); hold on; grid on;
    for k = 1:numel(methods)
        s = results.truth.(['batch_' methods{k}]);
        plot(time,s.radialM,'DisplayName',[upper(methods{k}) ' 径向']);
        plot(time,s.transverseM,'DisplayName',[upper(methods{k}) ' 横向']);
    end
    xlabel('相对时间 / s'); ylabel('误差 / m'); legend('show','Location','best');
end
if isfield(results,'dpe')
    obs = results.observations;
    figure('Name','批处理 DPE 残差'); hold on; grid on;
    for sig = unique(obs.signalId(results.dpe.used))'
        rows = results.dpe.used & obs.signalId == sig;
        plot(time(obs.epochIndex(rows)),results.dpe.residualM(rows),'.', ...
            'DisplayName',sprintf('signal %d',sig));
    end
    xlabel('相对时间 / s'); ylabel('伪距残差 / m'); legend('show','Location','best');
    if ~isempty(results.dpe.grid)
        figure('Name','DPE 最弱方向网格代价'); grid on;
        plot(results.dpe.grid.offsetsM/1000,results.dpe.grid.cost);
        xlabel('沿最弱方向的位置偏移 / km'); ylabel('加权 L1 代价 / m');
    end
end
end
