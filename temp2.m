clc;clear;
load(fullfile(fileparts(mfilename('fullpath')), 'matlab3.mat'))

n = numel(NAV);

X2 = zeros(n,1);
Y2 = zeros(n,1);
Z2 = zeros(n,1);

for i = 1:n
    X2(i) = NAV(i).posX;
    Y2(i) = NAV(i).posY;
    Z2(i) = NAV(i).posZ;
end

figure;
hold on;
grid on;
axis equal;



scatter3(X2, Y2, Z2, ...
         40, 'b', 'filled');

xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('ECEF Coordinates Comparison');



view(3);