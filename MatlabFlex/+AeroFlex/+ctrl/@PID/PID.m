classdef PID < ControllerBase
%PID  Baseline proportional–integral–derivative controller.
%
%   cfg structure must contain:
%       Kp, Ki, Kd   : gains (vectors size nu)
%       dt           : sample time
%       refFcn(t)    : desired reference trajectory (handle)

    properties (Access = private)
        Kp Ki Kd
        dt
        e_int
        refFcn
    end

    methods
        function obj = PID(cfg)
            obj@ControllerBase(cfg);
            obj.Kp = cfg.Kp; obj.Ki = cfg.Ki; obj.Kd = cfg.Kd;
            obj.dt = cfg.dt;
            obj.refFcn = cfg.refFcn;
            obj.e_int = zeros(size(obj.prevU));
        end

        function u = computeControl(obj,xhat,t)
            r  = obj.refFcn(t);
            e  = r - xhat;
            de = (e - obj.e_int)/obj.dt;  % use previous error for derivative
            obj.e_int = obj.e_int + e*obj.dt;

            u = obj.Kp.*e + obj.Ki.*obj.e_int + obj.Kd.*de;
            u = obj.saturate(u);
            obj.prevU = u;
        end
    end
end
