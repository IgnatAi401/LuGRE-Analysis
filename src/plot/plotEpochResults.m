function plotEpochResults(results)
% PLOTEPOCHRESULTS 逐历元结果：只显示开启的方法，标签明确为与 NAV 的差异。
cfg = results.config;
time = results.gpsSeconds-results.gpsSeconds(1);
figure('Name','LuGRE 与 NAV 的位置差'); hold on; grid on;
for method = {'dpe','ls'}
    name = method{1};
    if cfg.enable.(name)
        plot(time,results.(name).differenceNavM,'DisplayName',upper(name));
    end
end
if cfg.enable.orbit
    plot(time,results.orbit.differenceNavM,'DisplayName','轨道预测');
end
xlabel('相对时间 / s'); ylabel('与 NAV 的三维位置差 / m'); legend('show');
figure('Name','接收机钟差'); hold on; grid on;
for method = {'dpe','ls'}
    name = method{1};
    if cfg.enable.(name)
        plot(time,results.(name).clockM,'DisplayName',upper(name));
    end
end
xlabel('相对时间 / s'); ylabel('钟差 c·dt / m'); legend('show');
figure('Name','ECEF 轨迹'); hold on; grid on; axis equal;
p = results.navPosition;
plot3(p(:,1),p(:,2),p(:,3),'DisplayName','NAV');
for method = {'dpe','ls'}
    name = method{1};
    if cfg.enable.(name)
        p = results.(name).position;
        plot3(p(:,1),p(:,2),p(:,3),'DisplayName',upper(name));
    end
end
xlabel('X / m'); ylabel('Y / m'); zlabel('Z / m'); legend('show'); view(3);
if isfield(results,'truth'), plotEpochTruth(results); end
end
