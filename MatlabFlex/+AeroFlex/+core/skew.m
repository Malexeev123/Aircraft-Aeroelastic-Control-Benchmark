function S = skew(V)
    % Check the size of the input vector
    if length(V) == 2
        % Create a 2D skew-symmetric matrix
        S = [0, -V(1);
             V(1), 0];
    elseif length(V) == 3
        % Create a 3D skew-symmetric matrix
        S = [0, -V(3), V(2);
             V(3), 0, -V(1);
             -V(2), V(1), 0];
    else
        error('Input vector must be 2D or 3D.');
    end
end