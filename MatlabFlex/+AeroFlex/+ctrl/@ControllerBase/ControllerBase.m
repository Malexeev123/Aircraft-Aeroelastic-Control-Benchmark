classdef (Abstract) ControllerBase < handle
%CONTROLLERBASE  Abstract controller superclass for AeroFlexProject.
%
%   All controllers (PID, LQR, nMPC …) inherit from this class and
%   implement `computeControl`.  The interface is intentionally
%   minimal so that high‑level code stays generic.
%
%   Properties
%   ----------
%   prevU   :  last command vector   (m×1 double)
%   limits  :  [umin umax] or []     (m×2 double)
%   cfg     :  user‑supplied struct
%
%   Methods (Abstract)
%   ------------------
%   u = computeControl(obj,xhat,t)
%       • xhat : estimator output at tk
%       • t    : current time (double)

    properties
        prevU  double
        limits double
        cfg
    end

    methods
        function obj = ControllerBase(cfg)
            arguments, cfg struct, end
            obj.cfg    = cfg;
            obj.prevU  = zeros(cfg.ctrl.n_surf*cfg.ctrl.var_per,1);
            if isfield(cfg,'limits'), obj.limits = cfg.limits; end
        end
    end

    methods (Abstract)
        u = computeControl(obj,xhat,t)
    end

    methods (Access = protected)
        function u = saturate(obj,u)
            if ~isempty(obj.limits)
                u = min(max(u,obj.limits(:,1)),obj.limits(:,2));
            end
        end
    end
end
