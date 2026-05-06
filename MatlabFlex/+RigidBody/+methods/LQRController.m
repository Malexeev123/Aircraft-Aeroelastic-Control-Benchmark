classdef LQRController < handle
    % cascaded attitude-+velocity regulator (outer loop)
    properties
        A  double;  B double;  K double
        ref  double = zeros(6,1);   % [phi theta psi u v w]
    end
    methods
        function obj = LQRController(A,B,Q,R)
            [obj.A,obj.B] = aircraft.linearise(trim); % Change to get from trimmed A and B matrix from STM
            % Q = diag([20 20 20  5*ones(1,3)]);   % tune here
            % R = diag([1 1 1 1]);
            obj.K = lqr(obj.A,obj.B,Q,R);
            obj.xTrim = trim.x;
            obj.uTrim = trim.u;
        end
        function [u, du] = command(obj,xCmd,xMeas)
            dx = xMeas - xCmd;
            du = -obj.K*dx;
            u  = obj.uTrim + du;
        end
    end
end
function [u_cmd] = outerLQR(xRB,xRef,uTrim,params)
    q   = xRB(1:4);    pqr = xRB(5:7);   U = xRB(8:10);
    qErr = quatmultiply( qRefConj , q )';      % 3-vec part used
    xLQR = [ qErr(2:4) ; pqr ; U-Uref ];        % 9-state vector
    u_cmd = -params.K * xLQR + uTrim;           % K from dlqr()
end