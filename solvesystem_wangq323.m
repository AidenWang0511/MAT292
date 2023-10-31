function [t, x] = solvesystem_wangq323(f, g, t0, tN, x0, h)
    
    % Calculating number of steps
    N = round((tN - t0) / h) + 1;

    % Allocating space for t and y
    t = linspace(t0, tN, N);
    x = zeros(2, N);
    
    % Setting initial condition
    x(1, 1) = x0(1);
    x(2, 1) = x0(2);
    
    for i = 2:N
       % Predictor step (Euler method)
        x1_nxt = x(1, i-1) + h * f(t(i-1), x(1, i-1), x(2, i-1));
        x2_nxt = x(2, i-1) + h * g(t(i-1), x(1, i-1), x(2, i-1));

        % Corrector step (Euler method)
        x(1, i) = x(1, i-1) + 0.5 * h * (f(t(i-1), x(1, i-1), x(2, i-1)) + f(t(i), x1_nxt, x2_nxt));
        x(2, i) = x(2, i-1) + 0.5 * h * (g(t(i-1), x(1, i-1), x(2, i-1)) + g(t(i), x1_nxt, x2_nxt));
    end
end