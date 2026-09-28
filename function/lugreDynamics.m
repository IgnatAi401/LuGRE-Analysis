function derivative = lugreDynamics(t, state, anchorGps, cfg)
% LUGREDYNAMICS 地心 ICRF：航天器与月球的 [r;v]，单位 m、m/s。
r = state(1:3); v = state(4:6); moon = state(7:9);
muE = cfg.constants.muEarth; muM = cfg.constants.muMoon;
relative = r-moon;
a = -muE*r/norm(r)^3 - muM*relative/norm(relative)^3 - muM*moon/norm(moon)^3;
if cfg.enable.j2
    % 将 J2 在地固系计算后转回惯性系，不能把惯性 Z 轴当作当日地轴。
    rotation = lugreFrameRotation(anchorGps+t,cfg.frame);
    fixed = rotation*r;
    radius = norm(fixed); z2 = (fixed(3)/radius)^2;
    factor = -1.5*muE*cfg.constants.j2*cfg.constants.earthRadiusM^2/radius^5;
    j2 = factor * [fixed(1)*(1-5*z2);fixed(2)*(1-5*z2);fixed(3)*(3-5*z2)];
    a = a + rotation'*j2;
end
moonAcceleration = -(muE+muM)*moon/norm(moon)^3;
derivative = [v;a;state(10:12);moonAcceleration];
end
