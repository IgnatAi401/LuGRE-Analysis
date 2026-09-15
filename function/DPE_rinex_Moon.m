function [pos, dt, idx, dpos] = DPE_rinex_Moon(satPos, Grid_ECEF, Origin_ECEF)
%==========================================================================
% Function: DPE for MOON Receiver (with Ionospheric Correction)
% Inputs:
%       satPos      - [n x 5] double
%                     Cols 1-3: GPS Sat positions in Earth ECEF (X, Y, Z)
%                     Col  4:   Pseudorange observations (meters)
%                     Col  5:   C/N0 (dB-Hz)
%       Grid_ECEF   - Candidate positions [4 x N_grid]
%                     Rows 1-3: Pos X, Y, Z
%                     Row  4:   Clock Bias Candidate (c * dt)
%       Origin_ECEF - Current estimated center [3 x 1]
% Outputs:
%       pos         - Receiver position [3 x 1]
%       dt          - Receiver clock error (meters)
%       idx         - Index of best grid
%       dpos        - Offset from Origin_ECEF
%==========================================================================

    % == Constants ===========================================================
    vel_c   = 299792458;          % Speed of light (m/s)
    Omega_e = 7.2921151467e-5;    % Earth rotation rate (rad/s)
    pi      = 3.141592653589793;

    % ---- Ionosphere model (same as DPE_rinex) ----
    ion_Iz  = 5e-9;               % Zenith ionospheric delay (seconds)

    n_sat   = size(satPos, 1);
    n_grid  = size(Grid_ECEF, 2);

    val_all = zeros(1, n_grid);

    % == Main Loop ===========================================================
    for i = 1 : n_sat

        % --- Satellite data ---
        curr_sat_pos = satPos(i, 1:3)';   % 3x1
        curr_obs_pr  = satPos(i, 4);      % pseudorange (m)
        cn0          = satPos(i, 5);      % C/N0 (dB-Hz)

        % --- Sagnac correction ---
        traveltime = curr_obs_pr / vel_c;
        angle = traveltime * Omega_e;

        R = [ cos(angle)  sin(angle) 0;
             -sin(angle)  cos(angle) 0;
              0           0          1 ];

        Rot_satpos = R * curr_sat_pos;

        % --- Elevation angle computation ---
        los = Rot_satpos - Origin_ECEF;
        r   = norm(los);
        el  = asin(los(3) / r) * 180 / pi;   % elevation in degrees

        % --- Ionospheric correction ---
        ion_F = 1 + 16 * (0.53 - el/180)^3;
        ion   = ion_Iz * ion_F * vel_c;      % meters

        % --- Geometric distance to all grid points ---
        ri_candidate = vecnorm(Rot_satpos - Grid_ECEF(1:3, :));

        % --- Model pseudorange ---
        model_pr = ri_candidate + Grid_ECEF(4, :) + ion;

        % --- Cost accumulation (CN0 weighted L1 norm) ---
        val_all = val_all + 10^(cn0/10) * abs(curr_obs_pr - model_pr);

    end

    % == Minimum Search ======================================================
    [~, idx] = min(val_all);
    res  = Grid_ECEF(:, idx);

    pos  = res(1:3);
    dt   = res(4);
    dpos = pos - Origin_ECEF;

end
