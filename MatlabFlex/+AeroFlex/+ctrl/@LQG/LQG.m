classdef LQG < ControllerBase
%LQG  Linear‑Quadratic‑Gaussian controller = KF + LQR.
%
%   Expects cfg.estimator (an EstimatorBase) + LQR fields.
    properties (Access = private)
        K
        estimator EstimatorBase
    end

    methods
        function obj = LQG(cfg)
            obj@ControllerBase(cfg);
            obj.estimator = cfg.estimator;
            obj.K = lqr(cfg.A,cfg.B,cfg.Q,cfg.R);
        end

        function u = computeControl(obj,~,t) %#ok<INUSD>
            xhat = obj.estimator.xhat;   % get latest KF estimate
            u = -obj.K*xhat;
            u = obj.saturate(u);
            obj.prevU = u;
        end
    end
end
