function [img, x, y] = lunarBasemap(centerPa, limitKm, nPixels, cacheFolder)
% LUNARBASEMAP 以月面点为中心的当地东-北平面底图（LRO WAC 全月拼接影像，灰度 0~1）。
% centerPa 为月固系（DE405 主轴系 PA）位置 3x1（m），底图覆盖东、北各 ±limitKm，nPixels x nPixels 像素；
% 返回的 x、y（km）对应 img 的列、行，可直接 image(x, y, img)（需 YDir normal）。
% 正射投影：东-北平面上 (E, N) 对应月面点 p = 中心·sqrt(R²-E²-N²)/R + E·east + N·north（月面点的水平误差即其
% 在当地水平面上的投影，与俯视图误差定义一致）；E²+N² > R² 的区域在月盘外，返回 NaN。
% 影像在月球平均地球系（ME）中，PA -> ME 用 DE403/DE405 的常数旋转（Rx(-0.1462")Ry(-79.0768")Rz(-63.8986")）；
% 转换后 Blue Ghost 着陆点为 18.5621°N、61.8103°E，与公开着陆坐标相差约 25 m。
% 影像瓦片来自 NASA Moon Trek WMTS（LRO_WAC_Mosaic_Global_303ppd_v02），按需下载并缓存到 cacheFolder。
% 地理栅格参考与插值用 Mapping Toolbox（georefcells、geointerp）。
if nargin < 3 || isempty(nPixels), nPixels = 800; end
if ~exist(cacheFolder,'dir'), mkdir(cacheFolder); end
centerPa = centerPa(:);
up = centerPa/norm(centerPa);
lon0 = atan2(up(2), up(1)); lat0 = asin(up(3));
east = [-sin(lon0); cos(lon0); 0];
north = [-sin(lat0)*cos(lon0); -sin(lat0)*sin(lon0); cos(lat0)];

x = linspace(-limitKm, limitKm, nPixels);
y = linspace(-limitKm, limitKm, nPixels);
[E, N] = meshgrid(x*1000, y*1000);
radius = norm(centerPa);
p = centerPa*(sqrt(max(radius^2 - E(:)'.^2 - N(:)'.^2, 0))/radius) + east*E(:)' + north*N(:)';
outside = reshape(E(:).^2 + N(:).^2 > radius^2, size(E));
p = paToMe()*p;
lat = reshape(asind(p(3,:)./vecnorm(p)), size(E));
lon = reshape(atan2d(p(2,:), p(1,:)), size(E));
lat(outside) = NaN; lon(outside) = NaN;

% 瓦片层级：第 n 层每块 256 像素、跨 180/2^n 度；取像素尺寸不大于输出像素的最低层（最高 8 层，约 83 m/像素）。
neededPpd = nPixels/(2*limitKm) * radius/1000*pi/180;
level = min(max(ceil(log2(neededPpd*180/256)), 0), 8);
tileDeg = 180/2^level;
rows = floor((90 - max(lat(:)))/tileDeg) : floor((90 - min(lat(:)))/tileDeg);
cols = floor((min(lon(:)) + 180)/tileDeg) : floor((max(lon(:)) + 180)/tileDeg);
rows = rows(rows >= 0 & rows < 2^level);
cols = cols(cols >= 0 & cols < 2^(level+1));
mosaic = zeros(256*numel(rows), 256*numel(cols));
for i = 1:numel(rows)
    for j = 1:numel(cols)
        mosaic(256*(i-1)+(1:256), 256*(j-1)+(1:256)) = readTile(level, rows(i), cols(j), cacheFolder);
    end
end
R = georefcells(90 - [rows(end)+1, rows(1)]*tileDeg, -180 + [cols(1), cols(end)+1]*tileDeg, ...
    size(mosaic), 'ColumnsStartFrom', 'north');
img = geointerp(mosaic, R, lat, lon, 'linear');
end

function tile = readTile(level, row, col, cacheFolder)
file = fullfile(cacheFolder, sprintf('wac_%d_%d_%d.jpg', level, row, col));
if ~isfile(file)
    url = sprintf(['https://trek.nasa.gov/tiles/Moon/EQ/LRO_WAC_Mosaic_Global_303ppd_v02/1.0.0/' ...
        'default/default028mm/%d/%d/%d.jpg'], level, row, col);
    websave(file, url, weboptions('Timeout', 60));
end
tile = double(imread(file))/255;
if size(tile,3) > 1, tile = mean(tile, 3); end
end

function r = paToMe()
% DE403/DE405 月球主轴系 -> 平均地球系（角度单位：角秒）。
a = [0.1462 79.0768 63.8986]/3600;
rx = [1 0 0; 0 cosd(-a(1)) sind(-a(1)); 0 -sind(-a(1)) cosd(-a(1))];
ry = [cosd(-a(2)) 0 -sind(-a(2)); 0 1 0; sind(-a(2)) 0 cosd(-a(2))];
rz = [cosd(-a(3)) sind(-a(3)) 0; -sind(-a(3)) cosd(-a(3)) 0; 0 0 1];
r = rx*ry*rz;
end
