function [X_next, V_next] = RK4(dydx, dt, X_SC, V_SC, i)
% RK4 numerical solver
    dv1 = dt * dydx(X_SC(i,:));
    dx1 = dt * V_SC(i,:);

    dv2 = dt * dydx(X_SC(i,:) + 0.5 * dx1);
    dx2 = dt * (V_SC(i,:) + 0.5 * dv1);

    dv3 = dt * dydx(X_SC(i,:) + 0.5  * dx2);
    dx3 = dt * (V_SC(i,:) +  0.5 * dv2);

    dv4 = dt * dydx(X_SC(i,:) + dx3);
    dx4 = dt * (V_SC(i,:) + dv3);

    V_SC(i + 1,:) = V_SC(i,:) +  (dv1 + 2*dv2 + 2*dv3 + dv4) / 6;
    X_SC(i + 1,:) = X_SC(i,:) + (dx1 + 2*dx2 + 2*dx3 + dx4) / 6;

    X_next = X_SC;
    V_next = V_SC;
end