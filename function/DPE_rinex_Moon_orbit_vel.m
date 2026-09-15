function [pos, dt, idx, dpos] = DPE_rinex_Moon_orbit_vel(satPos, Grid_ECEF, Origin_ECEF, ve, dt_int, k)

%==========================================================================
% satPos structure
% 1-3  sat position (ECEF)
% 4-6  sat velocity (ECEF)
% 7    pseudorange
% 8    CN0
%
% Grid_ECEF
% 1-3  position
% 4    clock bias (m)
%
% ve
% receiver velocity vector (3x1)
%
% dt_int
% time interval between epochs
%==========================================================================

vel_c   = 299792458;
Omega_e = 7.2921151467e-5;
pi      = 3.141592653589793;

ion_Iz  = 5e-9;

n_sat   = size(satPos,1);
n_grid  = size(Grid_ECEF,2);

val_all = zeros(1,n_grid);

lambda_v = 20;   % velocity constraint weight (需要调)

for i = 1:n_sat

    %% satellite data
    sat_pos = satPos(i,1:3)';
    sat_vel = satPos(i,4:6)';
    obs_pr  = satPos(i,7);
    cn0     = satPos(i,8);

    %% Sagnac correction
    traveltime = obs_pr / vel_c;
    angle = traveltime * Omega_e;

    R = [ cos(angle)  sin(angle) 0;
         -sin(angle)  cos(angle) 0;
          0           0          1 ];

    sat_pos_rot = R * sat_pos;
    sat_vel_rot = R * sat_vel;

    %% elevation (for ionosphere)
    los0 = sat_pos_rot - Origin_ECEF;
    r0   = norm(los0);
    el   = asin(los0(3)/r0) * 180/pi;

    %% ionosphere
    ion_F = 1 + 16*(0.53 - el/180)^3;
    ion   = ion_Iz * ion_F * vel_c;

    %% LOS vector for all grid points
    los = sat_pos_rot - Grid_ECEF(1:3,:);
    range = vecnorm(los);

    %% pseudorange model
    model_pr = range + Grid_ECEF(4,:) + ion;

    %% pseudorange cost
    pr_cost = 10^(cn0/10) * abs(obs_pr - model_pr);

    %% ---------- velocity constraint ----------

    % receiver velocity replicated to grid
    ve_mat = repmat(ve,1,n_grid);

    % relative velocity
    rel_vel = sat_vel_rot - ve_mat;

    % range rate
    range_rate = sum(los .* rel_vel,1) ./ range;

    % predicted receiver motion range-rate
    v_grid = vecnorm(Grid_ECEF(1:3,:) - Origin_ECEF) / dt_int;

    vel_cost = lambda_v * abs(range_rate - v_grid);

    %% accumulate cost
    if k > 300
        val_all = val_all + pr_cost + vel_cost;
    else
        val_all = val_all + pr_cost;
    end

end

%% minimum search
[~,idx] = min(val_all);

res = Grid_ECEF(:,idx);

pos = res(1:3);
dt  = res(4);

dpos = pos - Origin_ECEF;

end