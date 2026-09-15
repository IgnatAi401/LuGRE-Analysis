function [pos, dt, dop, others] = leastSquarePos(data)

%% Constants
c  = 299792458;                 % speed of light (m/s)
pi = 3.141592653589793;

%% Parse input
satpos = data(:,1:3);           
obs    = data(:,7);             
CN0    = data(:,8);             
Ns     = size(data,1);

%% Initialization
% pos = zeros(4,1);
% pos = [-7.827397544975670e+07;5.440271540064030e+07;4.676309044066320e+07;-2.464356080937330e+05];
pos = [-8.827397544975670e+07;4.440271540064030e+07;3.676309044066320e+07;-140000];
numIter  = 7;

G   = zeros(Ns,4);
W   = eye(Ns);
v   = zeros(Ns,1);
el  = zeros(Ns,1);

%% Iterative least squares
for iter = 1:numIter
    for i = 1:Ns
        
        % geometric distance
        rho = norm(satpos(i,:)' - pos(1:3));
        
        % line-of-sight unit vector
        los = (satpos(i,:)' - pos(1:3)) / rho;
        
        % design matrix
        G(i,:) = [-los' 1];
        
        % % elevation angle
        % [~, el(i), ~] = topocent(pos(1:3), satpos(i,:)' - pos(1:3));
        % 
        % % ionospheric correction (simple model)
        % ion = 5e-9 * (1 + 16*(0.53 - el(i)/180)^3) * c;
        % 
        % % tropospheric correction
        % trop = tropo(sin(el(i)*pi/180), ...
        %              0.0, 1013.0, 293.0, 50.0, ...
        %              0.0, 0.0, 0.0);
        
        % observation minus computed
        % v(i) = obs(i) - rho - pos(4) - ion - trop;
        v(i) = obs(i) - rho - pos(4);
        
        % weight matrix (CN0-based)
        if iter > 2
            W(i,i) = CN0(i);
        end
    end
    
    % least squares solution
    % dX = (G' * W * G) \ (G' * W * v);
    dX = (G' * G) \ (G' * v);
    
    % update state
    pos = pos + dX;
end

%% Output
dt  = pos(4) / c;        % receiver clock bias (seconds)
pos = pos(1:3)';         % receiver position (ECEF)

%% DOP calculation (geometry only)
Q = inv(G' * G);
dop(1) = sqrt(trace(Q));                  % GDOP
dop(2) = sqrt(Q(1,1)+Q(2,2)+Q(3,3));      % PDOP
dop(3) = sqrt(Q(1,1)+Q(2,2));             % HDOP
dop(4) = sqrt(Q(3,3));                    % VDOP
dop(5) = sqrt(Q(4,4));                    % TDOP

%% Debug outputs
others.G  = G;
others.v  = v;
others.el = el;
others.dX = dX;

end
