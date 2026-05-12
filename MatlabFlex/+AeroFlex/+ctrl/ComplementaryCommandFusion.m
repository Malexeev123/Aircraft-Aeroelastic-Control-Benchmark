classdef ComplementaryCommandFusion < handle
%======================================================================
% COMPLEMENTARYCOMMANDFUSION
%======================================================================
% Fuses slow outer-loop command with fast inner-loop nMPC command:
%
%   u = uTrim + LP(uOuter-uTrim) + HP(uInner-uTrim)
%
% where:
%   LP(s) = wc/(s+wc)
%   HP(s) = s/(s+wc) = 1 - LP(s)
%
% Discrete implementation:
%   yLP[k] = a*yLP[k-1] + (1-a)*x[k],
%   a = exp(-wc*Ts)
%======================================================================

    properties
        Ts double
        wc double
        a double

        uTrim double

        yOuterLP double
        yInnerLP double

        initialized logical = false
        innerIsDeviation logical = false
    end

    methods
        function obj = ComplementaryCommandFusion(cfg, trim)
            obj.Ts = cfg.ctrl.Ts;

            if isfield(cfg,'fusion') && isfield(cfg.fusion,'omega_c')
                obj.wc = cfg.fusion.omega_c;
            else
                obj.wc = 2.5;
            end

            obj.a = exp(-obj.wc*obj.Ts);

            nSurf = cfg.ctrl.n_surf;
            varPer = cfg.ctrl.var_per;
            nu = nSurf*varPer;

            if isfield(cfg,'trim') && isfield(cfg.trim,'Uinpt')
                u0 = cfg.trim.Uinpt(:);
            elseif isfield(trim,'uTrim')
                u0 = trim.uTrim(:);
            else
                u0 = zeros(nu,1);
            end

            if numel(u0) == nSurf && varPer == 2
                u0 = [u0; zeros(nSurf,1)];
            elseif isscalar(u0)
                u0 = repmat(u0,nu,1);
            elseif numel(u0) ~= nu
                u0 = zeros(nu,1);
            end

            obj.uTrim = u0(:);

            obj.yOuterLP = zeros(nu,1);
            obj.yInnerLP = zeros(nu,1);

            if isfield(cfg,'fusion') && isfield(cfg.fusion,'innerIsDeviation')
                obj.innerIsDeviation = logical(cfg.fusion.innerIsDeviation);
            end
        end

        function [uCmd,info] = fuse(obj,uOuter,uInner)
            uOuter = uOuter(:);
            uInner = uInner(:);

            if numel(uOuter) ~= numel(obj.uTrim)
                error('ComplementaryCommandFusion:Dimension', ...
                      'uOuter length %d does not match uTrim length %d.', ...
                      numel(uOuter), numel(obj.uTrim));
            end

            if numel(uInner) ~= numel(obj.uTrim)
                error('ComplementaryCommandFusion:Dimension', ...
                      'uInner length %d does not match uTrim length %d.', ...
                      numel(uInner), numel(obj.uTrim));
            end

            dOuter = uOuter - obj.uTrim;

            if obj.innerIsDeviation
                dInner = uInner;
            else
                dInner = uInner - obj.uTrim;
            end

            obj.yOuterLP = obj.a*obj.yOuterLP + (1-obj.a)*dOuter;
            obj.yInnerLP = obj.a*obj.yInnerLP + (1-obj.a)*dInner;

            dOuterLow = obj.yOuterLP;
            dInnerHigh = dInner - obj.yInnerLP;

            uCmd = obj.uTrim + dOuterLow + dInnerHigh;

            info = struct();
            info.uOuter = uOuter;
            info.uInner = uInner;
            info.dOuterLow = dOuterLow;
            info.dInnerHigh = dInnerHigh;
            info.uCmd = uCmd;
            info.wc = obj.wc;
        end
    end
end