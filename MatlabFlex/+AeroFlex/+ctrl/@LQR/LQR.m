classdef LQR < ControllerBase
%LQR  Gain‑scheduled Linear‑Quadratic Regulator.
%
%   cfg must provide fields:
%       A,B,Q,R           : nominal LTI matrices  *or*
%       gainScheduleFcn() : handle returning K for a given operating point
%
%   If gainScheduleFcn exists it overrides the static design.

    properties (Access = private)
        K           % feedback gain matrix
        gainFcn
    end

    methods
        function obj = LQR(cfg)
            obj@ControllerBase(cfg);
            if isfield(cfg,'gainScheduleFcn')
                obj.gainFcn = cfg.gainScheduleFcn;
            else
                obj.K = lqr(cfg.A,cfg.B,cfg.Q,cfg.R);
            end
        end

        function u = computeControl(obj,xhat,t)
            if ~isempty(obj.gainFcn)
                obj.K = obj.gainFcn(xhat,t);
            end
            u = -obj.K*xhat;
            u = obj.saturate(u);
            obj.prevU = u;
        end
    end
end
