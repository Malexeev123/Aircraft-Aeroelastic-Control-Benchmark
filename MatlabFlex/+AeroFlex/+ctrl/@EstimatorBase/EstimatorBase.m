classdef (Abstract) EstimatorBase < handle
%ESTIMATORBASE  Abstract state–estimator superclass for AeroFlexProject.
%
%   All concrete observers *must* inherit from this class and
%   implement the `estimate` method so that SimRunner can swap them
%   simply by changing cfg.control.strategy.
%
%   Properties
%   ----------
%   xhat    :  latest state estimate  (n×1 double)
%   P       :  estimation‐error covariance (n×n double, optional)
%   cfg     :  user configuration struct passed at construction
%
%   Methods (Abstract)
%   ------------------
%   xhat = estimate(obj,z,u,t)
%       • z : measurement vector at tk
%       • u : control vector  at tk
%       • t : current time   (double)
%
%   The concrete method must update the internal state and return the
%   fresh x̂_k in *one* call.

    properties
        xhat   double   % state estimate
        P      double   % covariance (when applicable)
        yhat   double        
        what   double
        cfg             % arbitrary struct with user parameters
    end

    methods
        function obj = EstimatorBase(cfg)
            arguments, cfg struct, end
            obj.cfg = cfg;
        if isfield(cfg,'x0')
            obj.xhat = cfg.x0;
        elseif isfield(cfg,'trim')
            obj.xhat = cfg.trim.states;   % <<— new fallback
        else
            obj.xhat = zeros( cfg.nx ,1);
        end            
        if isfield(cfg,'P0'), obj.P    = cfg.P0; end
        end
    end

    methods (Abstract)
        % xhat = estimate(obj,z,u,t,varargin)

        xhat = estimate(obj,z,u,t)
    end
end
