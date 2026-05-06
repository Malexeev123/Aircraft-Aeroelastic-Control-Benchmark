classdef EKF < EstimatorBase
%EKF  Classic discrete Extended Kalman‑Filter implementation.

    properties (Access = private)
        f   % continuous dynamics  handle  ẋ = f(x,u,t)
        h   % measurement model    handle  z  = h(x,u,t)
        F   % Jacobian of f wrt x  handle
        H   % Jacobian of h wrt x  handle
        Q   % process‑noise cov
        R   % measurement‑noise cov
        dt  % integration step
    end

    methods
        function obj = EKF(cfg)
            obj@EstimatorBase(cfg);
            obj.f  = cfg.f;   obj.h  = cfg.h;
            obj.F  = cfg.F;   obj.H  = cfg.H;
            obj.Q  = cfg.Q;   obj.R  = cfg.R;
            obj.dt = cfg.dt;
        end

        function xhat = estimate(obj,z,u,t)
            %----- Prediction (1st‑order Euler) ------------------------
            xpred = obj.xhat + obj.f(obj.xhat,u,t)*obj.dt;
            Fk    = obj.F(obj.xhat,u,t);
            Ppred = Fk*obj.P*Fk.' + obj.Q;

            %----- Correction ------------------------------------------
            Hk    = obj.H(xpred,u,t);
            innov = z - obj.h(xpred,u,t);
            S     = Hk*Ppred*Hk.' + obj.R;
            K     = Ppred*Hk.'/S;          % optimal gain
            obj.xhat = xpred + K*innov;
            obj.P    = (eye(size(K,1))-K*Hk)*Ppred;

            xhat = obj.xhat;
        end
    end
end
