function [t, y] = DE2_wangq323(f, t0, tN, y0, y1, h)
    N = round((tN - t0) / h);
    
    % Allocate space for time vector
    t = linspace(t0, tN, N + 1);
    
    % Allocate space for y vector
    y = zeros(1, N + 1);
    y(1) = y0;
    y(2) = y0 + y1 * h;
    
    % Solving using loop
    for i = 2:N
        dy = (y(i) - y(i - 1)) / h;
        ddy = f(t(i), dy, y(i)); 
        y(i + 1) = (h^2) * ddy + 2 * y(i) - y(i - 1);
    end
end