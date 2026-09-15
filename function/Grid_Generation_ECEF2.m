function Grid_ECEF = Grid_Generation_ECEF2(Origin_ECEF, k)
   
    if k < 201
        dtmpx = -20000:1000:20000;
        dtmpy = -20000:1000:20000;
        dtmpz = -20000:1000:20000;
        dtmp_dt = -25000:1000:10000;
    else
        dtmpx = -10000:500:10000;
        dtmpy = -10000:500:10000;
        dtmpz = -10000:500:10000;
        dtmp_dt = -25000:1000:10000;
    end
    
    [dX, dY, dZ, dT] = ndgrid(dtmpx, dtmpy, dtmpz, dtmp_dt);
    
    % 展平成行向量
    dX = dX(:)';
    dY = dY(:)';
    dZ = dZ(:)';
    dT = dT(:)';
    
    dXYZ = [dX; dY; dZ];
    % dXYZT = [dX; dY; dZ;dT];
    % plot3(dX,dY,dZ,'o')
    
    dtmpdotx = dtmpx;
    dtmpdoty = dtmpy;
    dtmpdotz = dtmpz;
    dtmpdot_dt = dtmp_dt;
    
    [dXdot, dYdot, dZdot, dTdot] = ndgrid(dtmpdotx, dtmpdoty, dtmpdotz, dtmpdot_dt);
    dXdot = dXdot(:)';
    dYdot = dYdot(:)';
    dZdot = dZdot(:)';
    dTdot = dTdot(:)';
    
    dXYZdot = [dXdot; dYdot; dZdot];
    
    
    %% === Assemble full grid ==================================================
    
    Grid_ECEF = [
        Origin_ECEF + dXYZ;
        dT;
        dXYZdot;
        dTdot
        ];
end