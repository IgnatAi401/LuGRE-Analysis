function figures = plotLugreTruth(results)
% PLOTLUGRETRUTH 可复用于历史 MAT 的真值误差及轨迹图。
t = results.gpsSeconds(:)-results.gpsSeconds(1);
names = {'nav'}; labels = {'NAV'};
for name = {'dpe','ls','orbit'}
    if results.config.enable.(name{1})
        names{end+1} = name{1}; %#ok<AGROW>
        labels{end+1} = upper(name{1}); %#ok<AGROW>
    end
end
colors = [0.00 0.35 0.70; 0.85 0.25 0.05; 0.10 0.55 0.25; 0.55 0.15 0.70];
figures(1) = figure('Name','Ground truth - 3D error','Color','w','Theme','light','Position',[100 100 1050 600]);
hold on; grid on;
for k = 1:numel(names)
    v = results.truth.(names{k});
    plot(t,v.distanceM,'Color',colors(k,:),'LineWidth',1.2, ...
        'DisplayName',sprintf('%s (RMSE %.1f m)',labels{k},v.rmseM));
end
xlabel('Elapsed GPS time (s)'); ylabel('3D position error (m)');
title('Position error relative to ground truth'); legend('show','Location','best');
figures(2) = figure('Name','Ground truth - XYZ errors','Color','w','Theme','light','Position',[100 100 1050 800]);
layout = tiledlayout(3,1,'TileSpacing','compact');
axesNames = {'X','Y','Z'};
for j = 1:3
    nexttile; hold on; grid on;
    for k = 1:numel(names)
        plot(t,results.truth.(names{k}).errorM(:,j),'Color',colors(k,:), ...
            'LineWidth',1.1,'DisplayName',labels{k});
    end
    ylabel(sprintf('%s error (m)',axesNames{j}));
    if j==1, legend('show','Location','best'); end
end
xlabel(layout,'Elapsed GPS time (s)'); title(layout,'ITRF position minus ground truth');
figures(3) = figure('Name','Ground truth - trajectory','Color','w','Theme','light','Position',[100 100 1050 700]);
hold on; grid on; axis equal;
p = results.truth.position/1000;
plot3(p(:,1),p(:,2),p(:,3),'k--','LineWidth',2,'DisplayName','Ground truth');
for k = 1:numel(names)
    if strcmp(names{k},'nav'), p = results.referencePosition;
    else
        p = results.(names{k}).position;
        if isfield(results.(names{k}),'valid'), p(~results.(names{k}).valid,:) = NaN; end
    end
    p = p/1000;
    plot3(p(:,1),p(:,2),p(:,3),'Color',colors(k,:),'DisplayName',labels{k});
end
xlabel('ITRF X (km)'); ylabel('ITRF Y (km)'); zlabel('ITRF Z (km)');
view(3); legend('show','Location','best'); title('Trajectory comparison');
end
