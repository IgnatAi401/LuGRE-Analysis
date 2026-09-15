function [pos, dt, idx, dpos] = DPE_rinex_Moon_dop(satVel, velX, velY, velZ, Grid_ECEF, Origin_ECEF)
%==========================================================================
% Function: Doppler DPE for MOON Receiver
%
% Inputs:
%   satVel      - [n x 8]
%                 Col1-3 : Satellite position (ECEF)
%                 Col4-6 : Satellite velocity (ECEF)
%                 Col7   : Carrier wavelength (lambda)
%                 Col8   : C/N0 (dB-Hz)
%
%   velX,Y,Z    - Receiver velocity (ECEF)
%
%   Grid_ECEF   - Candidate positions [3 x N_grid]
%
%   Origin_ECEF - Current estimated center [3 x 1]
%
% Outputs:
%   pos         - Receiver position [3 x 1]
%   dt          - Receiver clock error (set to 0)
%   idx         - Index of best grid
%   dpos        - Offset from Origin_ECEF
%==========================================================================

    n_sat  = size(satVel,1);
    n_grid = size(Grid_ECEF,2);

    val_all = zeros(1,n_grid);

    % Receiver velocity vector
    vr = [velX; velY; velZ];

    % ================= Main Loop =================
    for i = 1:n_sat

        % ---- Satellite state ----
        rs = satVel(i,1:3)';      % position
        vs = satVel(i,4:6)';      % velocity
        lambda = satVel(i,7);     % carrier wavelength
        cn0 = satVel(i,8);        % C/N0

        % ---- LOS vectors to all grid points ----
        los = rs - Grid_ECEF(1:3,:);           % 3 x Ngrid
        r   = vecnorm(los);                    % distance

        % ---- unit LOS ----
        u = los ./ r;

        % ---- relative velocity ----
        v_rel = vs - vr;

        % ---- modeled Doppler ----
        % fd = -(1/lambda) * dot(v_rel, u)
        model_fd = -(1/lambda) * (v_rel' * u);

        % ---- observation (true Doppler assumed from origin) ----
        los0 = rs - Origin_ECEF;
        u0   = los0 / norm(los0);
        obs_fd = -(1/lambda) * dot(v_rel, u0);

        % ---- cost accumulation ----
        val_all = val_all + 10^(cn0/10) * abs(obs_fd - model_fd);

    end

    % ================= Minimum Search =================
    [~,idx] = min(val_all);

    pos = Grid_ECEF(1:3,idx);
    dt  = 0;   % Doppler version does not estimate clock bias

    dpos = pos - Origin_ECEF;

end