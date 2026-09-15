function plot_PR_projection_DPE(satPos, Grid_ECEF, Origin_ECEF, z_slice, clk_slice, res_th)
%==========================================================================
% Function: Plot pseudorange projection on DPE grid (diagnostic)
%
% Inputs:
%   satPos      - [n x 4] satellite position + pseudorange
%   Grid_ECEF   - [4 x N] grid (X,Y,Z,c*dt)
%   Origin_ECEF - [3 x 1] reference center
%   z_slice     - fixed Z value for slicing (meters)
%   clk_slice   - fixed clock bias (meters)
%   res_th      - residual threshold (meters), e.g. 1000
%==========================================================================

    vel_c   = 299792458;
    Omega_e = 7.2921151467e-5;

    n_sat = size(satPos,1);

    % === Select grid slice ==================================================
    idx = abs(Grid_ECEF(3,:) - z_slice) < 1e-3 & ...
          abs(Grid_ECEF(4,:) - clk_slice) < 1e-3;

    Grid2D = Grid_ECEF(:,idx);
    X = Grid2D(1,:);
    Y = Grid2D(2,:);

    figure; hold on; grid on;
    colormap jet;

    % === Loop satellites ====================================================
    for i = 1:n_sat

        curr_sat_pos = satPos(i,1:3)';
        curr_obs_pr  = satPos(i,4);

        % --- Sagnac correction ---
        traveltime = curr_obs_pr / vel_c;
        angle = traveltime * Omega_e;
        R = [ cos(angle)  sin(angle) 0;
             -sin(angle)  cos(angle) 0;
              0           0          1];
        Rot_satpos = R * curr_sat_pos;

        % --- Compute residual ---
        ri = vecnorm(Rot_satpos - Grid2D(1:3,:));
        model_pr = ri + Grid2D(4,:);
        residual = abs(curr_obs_pr - model_pr);

        % --- Plot allowed region ---
        scatter(X(residual < res_th), ...
                Y(residual < res_th), ...
                15, 'filled', 'DisplayName', ['Sat ',num2str(i)]);
    end

    xlabel('X (m)');
    ylabel('Y (m)');
    title('Pseudorange Projection on DPE Grid (Slice View)');
    legend;
    axis equal;
end
