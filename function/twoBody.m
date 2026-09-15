function [r_pred, v_pred] = twoBody(r, v, dt)

    mu = 3.986004418e14;

    x0 = [r; v];

    % 二体动力学微分方程
    function dx = dyn(~, x)
        r = x(1:3);
        v = x(4:6);

        a = -mu / norm(r)^3 * r;

        dx = [v; a];
    end

    % 数值积分
    [~, x] = ode45(@dyn, [0 dt], x0);

    r_pred = x(end,1:3)';
    v_pred = x(end,4:6)';

end