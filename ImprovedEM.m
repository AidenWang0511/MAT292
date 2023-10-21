function [t, y] = ImprovedEM(f, t0, tN, y0, h)
    %Inputs:
    % f: function
    % t0: start time
    % tN: end time
    % y0: init condition
    % h: step size

    %Finding # of points:
    N = round((tN - t0) / h) + 1;

    %Allocating space
    t = zeros(1, N);
    y = zeros(1, N);

    %init value
    t(1) = t0;
    y(1) = y0;

    %For loop to do the Improved EM
    for i= 1:N-1
        temp_y = y(i) + h * f(t(i), y(i));
        t(i+1) = t(i) + h;
        y(i+1) = y(i) + 0.5 * h * ( f(t(i), y(i)) + f(t(i+1), temp_y) ) ;
    end
end
        