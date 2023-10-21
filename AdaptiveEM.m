function [t, y] = AdaptiveEM(f, t0, tN, y0, h)
    t = t0;
    y = y0;
    tol = 1e-8;

    while(t(end) < tN)
        cur_t = t(end);
        cur_y = y(end);
        
        % Y from one Euler step of size h
        Y = cur_y + h * f(cur_t, cur_y);

        % Z from two successive Euler steps of size h/2
        Z = cur_y + h/2 * f(cur_t, cur_y);
        Z = Z + h/2 * f(cur_t + h/2, Z);
        
        % Error estimate of Y and Z
        D = Z - Y;

        if abs(D) < tol
            %then the step is successful
            t = [t, cur_t + h];
            y = [y, Z + D];
        else
            % Update the step size
            h = 0.9 * h * min(max(tol / abs(D), 0.3), 2);
        end
    end
end
