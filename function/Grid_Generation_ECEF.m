function Grid_ECEF = Grid_Generation_ECEF(Origin_ECEF)
%==========================================================================
% Grid_Generation_Moon_ECEF
%
% 3D DPE search grid for lunar navigation
% Coordinate frame: Moon-Centered Moon-Fixed (ECEF-like)
%
% Output grid structure:
% [X; Y; Z; dT; Xdot; Ydot; Zdot; dTdot]
%
% Author: (you)
%==========================================================================

%% === 3D position & clock bias search =====================================

% --- Position search domain (meters) ---
dpos = -20000:1000:20000;              % ±20 m, step 2 m
gridnum = length(dpos);

dZ = kron(dpos, ones(1, gridnum));
dY = kron(dZ, ones(1, gridnum));
dX = kron(dY, ones(1, gridnum));

dY = repmat(dY, 1, gridnum);
dZ = repmat(dZ, 1, gridnum^2);

dXYZ = [dX; dY; dZ];

% --- Clock bias search (meters equivalent) ---
% c * dt absorbed into range domain
dT = repmat(dpos * 0.6, 1, gridnum^3);

% %% === Velocity & clock drift search =======================================
% 
% % --- Velocity search domain (m/s) ---
dvel = -20:1:20;              % ±10 m/s
gridnum = length(dvel);

dZdot = kron(dvel, ones(1, gridnum));
dYdot = kron(dZdot, ones(1, gridnum));
dXdot = kron(dYdot, ones(1, gridnum));

dYdot = repmat(dYdot, 1, gridnum);
dZdot = repmat(dZdot, 1, gridnum^2);

dXYZdot = [dXdot; dYdot; dZdot];

% --- Clock drift search ---
dTdot = repmat(dvel * 0.25, 1, gridnum^3);

%% === Assemble full grid ==================================================

Grid_ECEF = [
    Origin_ECEF + dXYZ;
    dT;
    dXYZdot;
    dTdot
];

end
