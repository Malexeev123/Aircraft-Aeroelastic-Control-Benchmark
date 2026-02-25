classdef StateQuat
    % Compact rigid–body state using quaternions
    %   x = [ q0 q1 q2 q3      % attitude (unit quat, q0 scalar part)
    %         u  v  w          % body‐frame transl. vel.  [m/s]
    %         p  q  r          % body‐frame ang. vel.     [rad/s]
    %         Xe Ye Ze ]'      % inertial position        [m]
    %
    properties (Constant)
        N = 13;        % length of state vector
    end
    methods (Static)
        function q = eyeQuat,  q = [1 0 0 0].'; end
        function mat = asMatrix(x)
            % helper that returns a struct formatted for Simulink mux
            mat.q  = x(1:4);     mat.Vb = x(5:7);
            mat.Om = x(8:10);    mat.Pi = x(11:13);
        end
    end
end