function plotFilterResults(results)
% PLOTFILTERRESULTS 卡尔曼滤波结果：真值三维误差（含形式 3σ）、径向误差与每历元卫星数。
time = results.gpsSeconds - results.gpsSeconds(1);
methods = results.config.filter.methods;
if isfield(results,'truth')
    figure('Name','卡尔曼滤波 与真值的三维误差'); hold on; grid on; set(gca,'YScale','log');
    names = fieldnames(results.truth);
    names = names(~ismember(names,{'position','valid'}));
    for k = 1:numel(names)
        plot(time,results.truth.(names{k}).distanceM,'DisplayName',names{k});
    end
    for k = 1:numel(methods)
        plot(time,3*results.(methods{k}).formalStd3dM,'--','DisplayName',['形式 3\sigma ' methods{k}]);
    end
    xlabel('相对时间 / s'); ylabel('三维误差 / m'); legend('show','Location','best');
    figure('Name','卡尔曼滤波 径向误差与卫星数'); tiledlayout(2,1,'TileSpacing','compact');
    nexttile; hold on; grid on;
    for k = 1:numel(methods)
        plot(time,results.truth.(methods{k}).radialM/1e3,'DisplayName',methods{k});
    end
    ylabel('径向误差 / km'); legend('show','Location','best');
    nexttile; stairs(time,results.observationCount); grid on;
    xlabel('相对时间 / s'); ylabel('卫星数');
end
end
